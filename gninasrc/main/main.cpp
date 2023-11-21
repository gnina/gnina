#include <boost/program_options.hpp>

#include <iostream>
#include <string>
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <algorithm>
#include <iterator>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include <boost/lexical_cast.hpp>
#include <boost/assign.hpp>
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/timer/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>
#include <boost/thread/thread.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/unordered_map.hpp>

#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "cache.h"
#include "cache_gpu.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "non_cache_gpu.h"
#include "non_cache_cnn.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "quasi_newton.h"
#include "tee.h"
#include "custom_terms.h"
#include "cnn_scorer.h"
#include "coords.h"
#include "obmolopener.h"
#include "gpucode.h"
#include "precalculate_gpu.h"
#include "array3d.h"
#include "grid.h"
#include "molgetter.h"
#include "result_info.h"
#include "box.h"
#include "flexinfo.h"
#include "covinfo.h"
#include "builtinscoring.h"
#include "sem.h"
#include "user_opts.h"
#include "version.h"

#include <cuda_profiler_api.h>

using namespace boost::iostreams;
using boost::filesystem::path;

void doing(int verbosity, const std::string& str, tee& log)
    {
  if (verbosity > 1)
      {
    log << str << std::string(" ... ");
    log.flush();
  }
}

void done(int verbosity, tee& log)
    {
  if (verbosity > 1)
      {
    log << "done.";
    log.endl();
  }
}
std::string default_output(const std::string& input_name)
    {
  std::string tmp = input_name;
  if (tmp.size() >= 6 && tmp.substr(tmp.size() - 6, 6) == ".pdbqt")
    tmp.resize(tmp.size() - 6); // FIXME?
  return tmp + "_out.pdbqt";
}

void write_all_output(model& m, const output_container& out, sz how_many,
    std::ostream& outstream)
    {
  if (out.size() < how_many)
    how_many = out.size();
  VINA_FOR(i, how_many)
  {
    m.set(out[i].c);
    m.write_model(outstream, i + 1); // so that model numbers start with 1
  }
}

//set m to a random conformer
fl do_randomization(model& m, const vec& corner1,
    const vec& corner2, int seed, int verbosity, tee& log)
    {
  conf init_conf = m.get_initial_conf(false);
  rng generator(static_cast<rng::result_type>(seed));
  if (verbosity > 1)
      {
    log << "Using random seed: " << seed;
    log.endl();
  }
  const sz attempts = 100;
  conf best_conf = init_conf;
  fl best_clash_penalty = 0;
  VINA_FOR(i, attempts)
  {
    conf c = init_conf;
    c.randomize(corner1, corner2, generator);
    m.set(c);
    fl penalty = m.clash_penalty();
    if (i == 0 || penalty < best_clash_penalty)
        {
      best_conf = c;
      best_clash_penalty = penalty;
      if (penalty == 0) break;
    }
  }
  m.set(best_conf);
  if (verbosity > 1)
      {
    log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
    log.endl();
  }
  return best_clash_penalty;
}

void refine_structure(model& m, const precalculate& prec, non_cache& nc,
    output_type& out, const vec& cap, const minimization_params& minparm,
    grid& user_grid, int verbosity, tee& log,non_cache& nc_new)
    {
  // std::cout << m.get_name() << " | pose " << m.get_pose_num() << " | refining structure\n";
  change g(m.get_size(), nc.move_receptor());

  nc.adjust_center(m); //for cnn, set cnn box
  if (!nc.within(m)) {
    std::cout << m.get_name() << " | pose " << m.get_pose_num()
        << " | initial pose not within box\n";
  }
  quasi_newton quasi_newton_par(minparm);
  const fl slope_orig = nc.getSlope();
  //try 5 times to get ligand into box
  //dkoes - you don't need a very strong constraint to keep ligands in the box,
  //but this factor can really bias the energy landscape
  fl slope = 10;
  VINA_FOR(p, 5)
  {
    nc.setSlope(slope);
    quasi_newton_par(m, prec, nc, out, g, cap, user_grid); //quasi_newton operator
    m.set(out.c); // just to be sure
    if (nc.within(m)) {
      break;
    }
    std::cout << m.get_name() << " | pose " << m.get_pose_num()
        << " | ligand outside box\n";
    slope *= 10;
  }
  out.coords = m.get_heavy_atom_movable_coords();
  if (!nc.within(m))
    out.e = max_fl;
  nc.setSlope(slope_orig);
  if (verbosity>1){
    //log total and empirical energy, useful for testing CNN + empirical merge 
    log.endl();
    fl final_e = nc.eval_deriv(m, cap[1], user_grid);
    log << "Total energy after refinement: " << std::fixed << std::setprecision(5) << final_e; 
    log.endl();
    //non_cache::eval for empirical energy
    fl final_emp_e = nc_new.eval(m, cap[1]);
    log << "Empirical energy after refinement: " << std::fixed << std::setprecision(5) << final_emp_e; 
    log.endl();
  }
}

std::string vina_remark(fl e, fl lb, fl ub)
    {
  std::ostringstream remark;
  remark.setf(std::ios::fixed, std::ios::floatfield);
  remark.setf(std::ios::showpoint);
  remark << "REMARK VINA RESULT: " << std::setw(9) << std::setprecision(1)
      << e << "  " << std::setw(9) << std::setprecision(3) << lb << "  "
      << std::setw(9) << std::setprecision(3) << ub << '\n';
  return remark.str();
}

output_container remove_redundant(const output_container& in, fl min_rmsd)
{
  output_container tmp;

  VINA_FOR_IN(i, in) {
    std::pair<sz, fl> closest_rmsd = find_closest(in[i].coords, tmp);
    if (closest_rmsd.first >= tmp.size() || closest_rmsd.second > min_rmsd) {
      tmp.push_back(new output_type(in[i])); //not redundant
    }
  }
  return tmp;
}

//print info to log about cnn scoring
static void get_cnn_info(model& m, CNNScorer& cnn, tee& log, float& cnnscore,
    float& cnnaffinity, float& cnnvariance) {
  float loss = 0;
  cnnscore = 0;
  cnnaffinity = 0;
  cnnscore = cnn.score(m, false, cnnaffinity, loss, cnnvariance);
  if (cnnscore >= 0 && cnn.options().moving_receptor())
  {
    //recalculate with ligand at center
    if (cnn.options().verbose) {
      log << "CNNscore1: " << std::fixed << std::setprecision(10) << cnnscore;
      log.endl();
      log << "CNNaffinity1: " << std::fixed << std::setprecision(10)
          << cnnaffinity;
      log.endl();
    }

    cnn.set_center_from_model(m);
    cnnscore = cnn.score(m, false, cnnaffinity, loss, cnnvariance);
  }

  if (cnn.options().verbose) {
    log << "CNNscore: " << std::fixed << std::setprecision(10) << cnnscore;
    log.endl();
    log << "CNNaffinity: " << std::fixed << std::setprecision(10)
        << cnnaffinity;
    log.endl();
  }
}

//dkoes - return all energies and rmsds to original conf with result
void do_search(model& m, const boost::optional<model>& ref,
    const weighted_terms& sf, const precalculate& prec, igrid& ig,
    non_cache& nc, // nc.slope is changed
    const vec& corner1, const vec& corner2,
    const parallel_mc& par, const user_settings& settings,
    bool compute_atominfo, tee& log,
    const terms *t, grid& user_grid, CNNScorer& cnn,
    std::vector<result_info>& results,szv_grid_cache& grid_cache, const grid_dims& gd, fl slope)
{
  boost::timer::cpu_timer time;
  try {

    precalculate_exact exact_prec(sf); //use exact computations for final score
    conf_size s = m.get_size();
    conf c = m.get_initial_conf(nc.move_receptor());
    fl e = max_fl;
    fl intramolecular_energy = max_fl;
    fl cnnscore = 0, cnnaffinity = 0, cnnvariance = 0;
    fl rmsd = 0;
    if (settings.gpu_docking)
      m.initialize_gpu();
    const vec authentic_v(settings.forcecap, settings.forcecap,
        settings.forcecap); //small cap restricts initial movement from clash

    cnn.set_center_from_model(m);
    non_cache nc_new = non_cache(grid_cache, gd, &prec, slope);
    if (settings.score_only)
    {
      cnn.freeze_receptor();
      intramolecular_energy = m.eval_intramolecular(exact_prec,
          authentic_v, c);
      naive_non_cache nnc(&exact_prec); // for out of grid issues
      e = m.eval_adjusted(sf, exact_prec, nnc, authentic_v, c,
          intramolecular_energy, user_grid);
      get_cnn_info(m, cnn, log, cnnscore, cnnaffinity, cnnvariance);

      log << "Affinity: " << std::fixed << std::setprecision(5) << e
          << " (kcal/mol)\n";

      log << "CNNscore: " << std::fixed << std::setprecision(5) << cnnscore << " "
          << "\nCNNaffinity: " << cnnaffinity;
      if(cnnvariance > 0) {
        log << "\nCNNvariance: " << std::fixed << std::setprecision(5) << cnnvariance;
      }
      log.endl();

      std::vector<flv> atominfo;
      flv term_values = t->evale_robust(m);
      log << "Intramolecular energy: " << std::fixed << std::setprecision(5)
          << intramolecular_energy << "\n";

      log
      << "Term values, before weighting:\n";
      log << std::setprecision(5);
      log << "## " << boost::replace_all_copy(m.get_name(), " ", "_");

      VINA_FOR_IN(i, term_values)
      {
        log << ' ' << term_values[i];
      }

      conf_independent_inputs in(m);
      const flv nonweight(1, 1.0);
      for (unsigned i = 0, n = t->conf_independent_terms.size(); i < n; i++) {
        flv::const_iterator pos = nonweight.begin();
        log << " " <<
            t->conf_independent_terms[i].eval(in, (fl) 0.0, pos);
      }
      log << '\n';

      results.push_back(result_info(e, cnnscore, cnnaffinity, cnnvariance, -1, m));

      if (compute_atominfo)
        results.back().setAtomValues(m, &sf);
    }
    else if (settings.local_only)
    {
      vecv origcoords = m.get_heavy_atom_movable_coords();
      output_type out(c, e);
      doing(settings.verbosity, "Performing local search", log);
      refine_structure(m, prec, nc, out, authentic_v, par.mc.ssd_par.minparm,
          user_grid,settings.verbosity,log,nc_new);
      done(settings.verbosity, log);
      m.set(out.c);

      //be as exact as possible for final score
      naive_non_cache nnc(&exact_prec); // for out of grid issues

      fl intramolecular_energy = m.eval_intramolecular(exact_prec,
          authentic_v,
          out.c);
      e = m.eval_adjusted(sf, exact_prec, nnc, authentic_v, out.c,
          intramolecular_energy, user_grid);

      //reset the center after last call to set
      get_cnn_info(m, cnn, log, cnnscore, cnnaffinity, cnnvariance);

      vecv newcoords = m.get_heavy_atom_movable_coords();
      assert(newcoords.size() == origcoords.size());
      for (unsigned i = 0, n = newcoords.size(); i < n; i++) {
        rmsd += (newcoords[i] - origcoords[i]).norm_sqr();
      }
      rmsd /= newcoords.size();
      rmsd = sqrt(rmsd);
      log << "Affinity: " << std::fixed << std::setprecision(5) << e << "  "
          << intramolecular_energy
          << " (kcal/mol)\nRMSD: " << rmsd << "\n";
      log << "CNNscore: " << std::fixed << std::setprecision(5) << cnnscore << " "
          << "\nCNNaffinity: " << cnnaffinity;
      if(cnnvariance > 0) {
        log << "\nCNNvariance: " << std::fixed << std::setprecision(5) << cnnvariance;
      }
      log.endl();

      if (!nc.within(m))
        log << "WARNING: not all movable atoms are within the search space\n";

      done(settings.verbosity, log);
      results.push_back(result_info(e, cnnscore, cnnaffinity, cnnvariance, rmsd, m));

      if (compute_atominfo)
        results.back().setAtomValues(m, &sf);
    }
    else //docking
    {
      rng generator(static_cast<rng::result_type>(settings.seed));
      log << "Using random seed: " << settings.seed;
      log.endl();

      output_container out_cont;
      doing(settings.verbosity, "Performing search", log);
      par(m, out_cont, prec, ig, corner1, corner2, generator, user_grid,nc);
      done(settings.verbosity, log);
      doing(settings.verbosity, "Refining results", log);

      
      VINA_FOR_IN(i, out_cont) {
        if (settings.cnnopts.cnn_scoring==CNNmetropolisrescore)//don't refine with cnn if rescoring
        {
          refine_structure(m, prec, nc_new, out_cont[i], authentic_v,
              par.mc.ssd_par.minparm, user_grid,settings.verbosity,log,nc_new);
        }
        else
        refine_structure(m, prec, nc, out_cont[i], authentic_v,
              par.mc.ssd_par.minparm, user_grid,settings.verbosity,log,nc_new);
        
        get_cnn_info(m, cnn, log, cnnscore, cnnaffinity, cnnvariance);

        out_cont[i].cnnscore = cnnscore;
        out_cont[i].cnnaffinity = cnnaffinity;
        out_cont[i].cnnvariance = cnnvariance;

        if (not_max(out_cont[i].e)) {
            intramolecular_energy = m.eval_intramolecular(exact_prec, authentic_v, out_cont[i].c);
            //we want vina energies not CNN
            out_cont[i].e = m.eval_adjusted(sf, exact_prec, nc_new, authentic_v, out_cont[i].c, intramolecular_energy, user_grid);
            out_cont[i].intramol = intramolecular_energy;
        }
      }

      auto sorter = [settings](const output_type& lhs, const output_type& rhs) {
        switch(settings.sort_order) {
        case Energy:
          return lhs.e < rhs.e;
        case CNNaffinity:
          return lhs.cnnaffinity > rhs.cnnaffinity; //reverse
        case CNNscore:
        default:
          return lhs.cnnscore > rhs.cnnscore; //reverse
        }
      };

      out_cont.sort(sorter);
      out_cont = remove_redundant(out_cont, settings.out_min_rmsd);

      done(settings.verbosity, log);

      log.setf(std::ios::fixed, std::ios::floatfield);
      log.setf(std::ios::showpoint);
      log << '\n';
      log << "mode |  affinity  |  intramol  |    CNN     |   CNN\n";
      log << "     | (kcal/mol) | (kcal/mol) | pose score | affinity\n";
      log << "-----+------------+------------+------------+----------\n";

      model best_mode_model = m;
      if (!out_cont.empty())
        best_mode_model.set(out_cont.front().c);

      sz how_many = 0;
      VINA_FOR_IN(i, out_cont)
      {
        if (!not_max(out_cont[i].e))
          continue; //may sort by something other than energy, so do not break
        if (how_many >= settings.num_modes)
          break; // check energy_range sanity FIXME
        ++how_many;
        m.set(out_cont[i].c);
        log << std::setw(5) << how_many << std::setw(12)
            << std::setprecision(2) << out_cont[i].e 
            << std::setw(12) << std::setprecision(2) << out_cont[i].intramol; 
        log << " " << std::setw(12) << std::setprecision(4) << out_cont[i].cnnscore << "  "
            << std::setw(9) << std::setprecision(3) << out_cont[i].cnnaffinity;
        log.endl();

        //dkoes - setup result_info
        results.push_back(
            result_info(out_cont[i].e, out_cont[i].cnnscore, out_cont[i].cnnaffinity, out_cont[i].cnnvariance, -1, m));

        if (compute_atominfo)
          results.back().setAtomValues(m, &sf);

      }
      done(settings.verbosity, log);

      if (how_many < 1)
          {
        log
            << "WARNING: Could not find any conformations completely within the search space.\n"
            << "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
        log.endl();
      }
    }
  } catch(numerical_error& ne) {
    log << "ERROR processing " << m.get_name() << "\n";
    log << ne.what() << "\n";
  }
  //std::cout << "Refine time " << time.elapsed().wall / 1000000000.0 << "\n";
}

void load_ent_values(const grid_dims& gd, std::istream& user_in,
    array3d<fl>& user_data)
    {
  std::string line;
  user_data = array3d<fl>(gd[0].n + 1, gd[1].n + 1, gd[2].n + 1);

  for (sz z = 0; z < gd[2].n + 1; z++)
      {
    for (sz y = 0; y < gd[1].n + 1; y++)
        {
      for (sz x = 0; x < gd[0].n + 1; x++)
          {
        std::getline(user_in, line);
        user_data(x, y, z) = ::atof(line.c_str());
      }
    }
  }
  std::cout << user_data(gd[0].n - 3, gd[1].n, gd[2].n) << "\n";
}

void main_procedure(model &m, precalculate &prec,
    const boost::optional<model> &ref, // m is non-const (FIXME?)
    const user_settings &settings,
    bool no_cache, bool compute_atominfo,
    const grid_dims &gd, minimization_params minparm,
    const weighted_terms &wt, tee &log,
    std::vector<result_info> &results, grid &user_grid, CNNScorer &cnn)
{
  doing(settings.verbosity, "Setting up the scoring function", log);

  done(settings.verbosity, log);
  log << std::fixed << std::setprecision(10);

  vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
  vec corner2(gd[0].end, gd[1].end, gd[2].end);

  parallel_mc par;
  sz heuristic = m.num_movable_atoms()
      + 10 * m.get_size().num_degrees_of_freedom();
  par.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
  if (settings.num_mc_steps > 0) {
    par.mc.num_steps = settings.num_mc_steps;
  }
  if (settings.max_mc_steps > 0 && par.mc.num_steps > settings.max_mc_steps) {
    par.mc.num_steps = settings.max_mc_steps;
  }
  if (settings.temperature > 0) {
    par.mc.temperature = settings.temperature; //expose temperature to user for cnn metropolis
  }

  par.mc.ssd_par.evals = unsigned((25 + m.num_movable_atoms()) / 3);
  if (minparm.maxiters == 0)
    minparm.maxiters = par.mc.ssd_par.evals;
  par.mc.ssd_par.minparm = minparm;
  par.mc.min_rmsd = 1.0;
  par.mc.num_saved_mins = settings.num_modes > settings.num_mc_saved ? settings.num_modes : settings.num_mc_saved;
  par.mc.hunt_cap = vec(10, 10, 10);
  par.num_tasks = settings.exhaustiveness;
  par.num_threads = settings.cpu;
  par.display_progress = true;

  szv_grid_cache gridcache(m, prec.cutoff_sqr());
  const fl slope = 1e3; // FIXME: too large? used to be 100
  if (settings.randomize_only)
  {
    for (unsigned i = 0; i < settings.num_modes; i++) {
      fl e = do_randomization(m, corner1, corner2, settings.seed + i,
          settings.verbosity, log);
      results.push_back(result_info(e, -1, 0, 0, -1, m));
    }
    return;
  }
  else
  {
    non_cache *nc = NULL;
    if (settings.cnnopts.cnn_scoring >= CNNrefinement) {
      nc = new non_cache_cnn(gridcache, gd, &prec, slope, cnn);
    } else if(settings.gpu_docking) {
      log << "WARNING: --gpu with empirical scoring is experimental and not recommended\n";
      precalculate_gpu *gprec = dynamic_cast<precalculate_gpu*>(&prec);
      if (!gprec)
        abort();
      nc = new non_cache_gpu(gridcache, gd, gprec, slope);
    } else {
      nc = new non_cache(gridcache, gd, &prec, slope);
    }

    if (no_cache || settings.cnnopts.cnn_scoring == CNNall)  {
      do_search(m, ref, wt, prec, *nc, *nc, corner1, corner2, par,
          settings, compute_atominfo, log,
          wt.unweighted_terms(), user_grid, cnn,
          results,gridcache,gd,slope);
    }
    else
    {
      bool cache_needed = !(settings.score_only || settings.randomize_only
          || settings.local_only);

      if (cache_needed)
        doing(settings.verbosity, "Analyzing the binding site", log);
      std::unique_ptr<cache> c(
          (settings.gpu_docking) ?
              new cache_gpu("scoring_function_version001",
                  gd, slope, dynamic_cast<precalculate_gpu*>(&prec)) :
              new cache("scoring_function_version001", gd, slope));
      if (cache_needed)
      {
        std::vector<smt> atom_types_needed;
        m.get_movable_atom_types(atom_types_needed);
        c->populate(m, prec, atom_types_needed, user_grid);
        done(settings.verbosity, log);
      }
      do_search(m, ref, wt, prec, *c, *nc, corner1, corner2, par,
          settings, compute_atominfo, log,
          wt.unweighted_terms(), user_grid, cnn, results,gridcache,gd,slope);
    }

    delete nc;
  }
}

struct options_occurrence
{
    bool some;
    bool all;
    options_occurrence()
        :
            some(false), all(true)
    {
    } // convenience
    options_occurrence& operator+=(const options_occurrence& x)
        {
      some = some || x.some;
      all = all && x.all;
      return *this;
    }
};

options_occurrence get_occurrence(boost::program_options::variables_map& vm,
    boost::program_options::options_description& d)
    {
  options_occurrence tmp;
  VINA_FOR_IN(i, d.options())
  {
    const std::string& str = (*d.options()[i]).long_name();
    if ((str.substr(0, 4) == "size" || str.substr(0, 6) == "center"))
    {
      if (vm.count(str))
        tmp.some = true;
      else
        tmp.all = false;
    }
  }
  return tmp;
}

void check_occurrence(boost::program_options::variables_map& vm,
    boost::program_options::options_description& d)
    {
  VINA_FOR_IN(i, d.options())
  {
    const std::string& str = (*d.options()[i]).long_name();
    if ((str.substr(0, 4) == "size" || str.substr(0, 6) == "center")
        && !vm.count(str))
      std::cerr << "Required parameter --" << str << " is missing!\n";
  }
}

template<class T>
inline void read_atomconstants_field(smina_atom_type::info& info,
    T (smina_atom_type::info::*field), unsigned int line,
    const std::string& field_name, std::istream& in)
    {
  if (!(in >> (info.*field)))
  {
    throw usage_error(
        "Error at line " + boost::lexical_cast<std::string>(line)
            + " while reading field '" + field_name
            + "' from the atom constants file.");
  }
}

void setup_atomconstants_from_file(const std::string& atomconstants_file)
    {
  std::ifstream file(atomconstants_file.c_str());
  if (file)
  {
    // create map from atom type names to indices
    boost::unordered_map<std::string, unsigned> atomindex;
    for (size_t i = 0u; i < smina_atom_type::NumTypes; ++i)
        {
      atomindex[smina_atom_type::default_data[i].smina_name] = i;
    }

    //parse each line of the file
    std::string line;
    unsigned lineno = 1;
    while (std::getline(file, line))
    {
      std::string name;
      std::stringstream ss(line);

      if (line.length() == 0 || line[0] == '#')
        continue;

      ss >> name;

      if (atomindex.count(name))
          {
        unsigned i = atomindex[name];
        smina_atom_type::info& info = smina_atom_type::data[i];

        //change this atom's parameters
#				define read_field(field) read_atomconstants_field(info, &smina_atom_type::info::field, lineno, #field, ss)
        read_field(ad_radius);
        read_field(ad_depth);
        read_field(ad_solvation);
        read_field(ad_volume);
        read_field(covalent_radius);
        read_field(xs_radius);
        read_field(xs_hydrophobe);
        read_field(xs_donor);
        read_field(xs_acceptor);
        read_field(ad_heteroatom);
#				undef read_field
      }
      else
      {
        std::cerr << "Line " << lineno << ": ommitting atom type name "
            << name << "\n";
      }
      lineno++;
    }
  }
  else
    throw usage_error(
        "Error opening atom constants file:  " + atomconstants_file);
}

void print_atom_info(std::ostream& out)
{
  out
      << "#Name radius depth solvation volume covalent_radius xs_radius xs_hydrophobe xs_donor xs_acceptr ad_heteroatom\n";
  VINA_FOR(i, smina_atom_type::NumTypes)
  {
    smina_atom_type::info& info = smina_atom_type::data[i];
    out << info.smina_name;
    out << " " << info.ad_radius;
    out << " " << info.ad_depth;
    out << " " << info.ad_solvation;
    out << " " << info.ad_volume;
    out << " " << info.covalent_radius;
    out << " " << info.xs_radius;
    out << " " << info.xs_hydrophobe;
    out << " " << info.xs_donor;
    out << " " << info.xs_acceptor;
    out << " " << info.ad_heteroatom;
    out << "\n";
  }
}

const fl box_granularity = 0.375;

//set grid dims to match center/size
static void setup_grid_dims(fl center_x, fl center_y, fl center_z, fl size_x, fl size_y, fl size_z, grid_dims& gd) {
  vec span(size_x, size_y, size_z);
  vec center(center_x, center_y, center_z);
  VINA_FOR_IN(i, gd)
  {
    gd[i].n = sz(std::ceil(span[i] / box_granularity));
    fl real_span = box_granularity * gd[i].n;
    gd[i].begin = center[i] - real_span / 2;
    gd[i].end = gd[i].begin + real_span;
  }
}
void setup_user_gd(grid_dims& gd, std::ifstream& user_in)
{
  std::string line;
  size_t pLines = 3;
  std::vector<std::string> temp;
  fl center_x = 0, center_y = 0, center_z = 0, size_x = 0, size_y = 0,
      size_z = 0;

  for (; pLines > 0; --pLines) //Eat first 3 lines
    std::getline(user_in, line);
  pLines = 3;

//Read in SPACING
  std::getline(user_in, line);
  boost::algorithm::split(temp, line, boost::algorithm::is_space());
  const fl granularity = ::atof(temp[1].c_str());
//Read in NELEMENTS
  std::getline(user_in, line);
  boost::algorithm::split(temp, line, boost::algorithm::is_space());
  size_x = (::atof(temp[1].c_str()) + 1) * granularity; // + 1 here?
  size_y = (::atof(temp[2].c_str()) + 1) * granularity;
  size_z = (::atof(temp[3].c_str()) + 1) * granularity;
//Read in CENTER
  std::getline(user_in, line);
  boost::algorithm::split(temp, line, boost::algorithm::is_space());
  center_x = ::atof(temp[1].c_str()) + 0.5 * granularity;
  center_y = ::atof(temp[2].c_str()) + 0.5 * granularity;
  center_z = ::atof(temp[3].c_str()) + 0.5 * granularity;

  vec span(size_x, size_y, size_z);
  vec center(center_x, center_y, center_z);
  VINA_FOR_IN(i, gd)
  {
    gd[i].n = sz(std::ceil(span[i] / granularity));
    fl real_span = granularity * gd[i].n;
    gd[i].begin = center[i] - real_span / 2;
    gd[i].end = gd[i].begin + real_span;
  }

}

// work queue job format
struct worker_job {
  unsigned int molid;
  model *m;
  std::vector<result_info> *results;
  grid_dims gd;

  worker_job(unsigned int molid, model *m, std::vector<result_info> *results, grid_dims gd)
      : molid(molid), m(m), results(results), gd(gd){};

  worker_job() : molid(0), m(NULL), results(NULL) {
    for (int i = 0; i < 3; i++) {
      gd[i] = grid_dim();
    }
  };
};

// writer queue job format
struct writer_job {
  unsigned int molid;
  std::vector<result_info> *results;

  writer_job(unsigned int molid, std::vector<result_info> *results) : molid(molid), results(results){};

  writer_job() : molid(0), results(NULL){};
};

template <typename T> struct job_queue {
  job_queue() : jobs(0){};

  void push(T &job) {
    jobs.push(job);
    has_work.signal();
  }

  // Returns false and doesn't modify job iff the queue has been
  // closed. Should not be called again in the same thread afterwards.
  bool wait_and_pop(T &job) {
    has_work.wait();
    return !jobs.pop(job);
  }

  // Signal that all jobs are done. num_possible_waiters will be waiting
  // on the sem, so arrange for them to wake. They'll find an empty
  // queue that distinguishes these special signals.
  void close(size_t num_possible_waiters) {
    for (size_t i = 0; i < num_possible_waiters; i++)
      has_work.signal();
  }

  sem has_work;
  boost::lockfree::queue<T> jobs;
};

//A struct of parameters that define the current run. These are packed together
//because of boost's restriction on the number of arguments you can 
//give to bind (max args is 9, but I need 10+ for the following thread
//functions) so I can reduce the number of args I pass.
struct global_state
{
    user_settings* settings;
    boost::shared_ptr<precalculate> prec;
    minimization_params* minparms;
    weighted_terms* wt;
    grid* user_grid;
    tee* log;
    std::ofstream* atomoutfile;
    cnn_options cnnopts;

    global_state(user_settings* settings, boost::shared_ptr<precalculate> prec,
        minimization_params* minparms, weighted_terms* wt,
        grid* user_grid, tee* log, std::ofstream* atomoutfile, const cnn_options& co):
        settings(settings), prec(prec), minparms(minparms), wt(wt),
            user_grid(user_grid), log(log), atomoutfile(atomoutfile),
            cnnopts(co)
    {
    }
    ;
};

//function to occupy the worker threads with individual ligands from the work queue
//TODO: see if implementing weight sharing between CNNScorer instances results
//in enough memory efficiency to avoid using a single one
void threads_at_work(job_queue<worker_job> *wrkq,
    job_queue<writer_job> *writerq, global_state *gs,
    MolGetter *mols, int *nligs, CNNScorer cnn_scorer) //copy cnn_scorer so it can maintain state
{
  if(!gs->settings->no_gpu)
    initializeCUDA(gs->settings->device);
  
  if (gs->settings->gpu_docking)
    thread_buffer.init(available_mem(gs->settings->cpu));

  worker_job j;
  while (!wrkq->wait_and_pop(j))
  {
    __sync_fetch_and_add(nligs, 1);

    main_procedure(*(j.m), *gs->prec, boost::optional<model>(),
        *gs->settings,
        false, // no_cache == false
        gs->atomoutfile->is_open()
            || gs->settings->include_atom_info, j.gd,
        *gs->minparms, *gs->wt, *gs->log, *(j.results),
        *gs->user_grid, cnn_scorer);

    writer_job k(j.molid, j.results);
    writerq->push(k);
    delete j.m;
  }
}

void write_out(std::vector<result_info> &results, ozfile &outfile,
    std::string &outext,
    user_settings &settings, const weighted_terms &wt, ozfile &outflex,
    std::string &outfext, std::ofstream &atomoutfile)  {
  if (outfile)
  {
    //write out molecular data
    for (unsigned j = 0, nr = results.size(); j < nr; j++) {
      results[j].write(outfile, outext, settings.include_atom_info, &wt,
          j + 1);
    }
  }
  if (outflex)
  {
    //write out flexible residue data data
    for (unsigned j = 0, nr = results.size(); j < nr; j++) {
      results[j].writeFlex(outflex, outfext, j + 1);
    }
  }
  if (atomoutfile)
  {
    for (unsigned j = 0, m = results.size(); j < m; j++) {
      results[j].writeAtomValues(atomoutfile, &wt);
    }
  }
}

//function for the writing thread to write ligands in order to output file
void thread_a_writing(job_queue<writer_job>* writerq,
    global_state* gs,
    ozfile* outfile, std::string* outext, ozfile* outflex,
    std::string* outfext,
    int* nligs) {
  try {
    int nwritten = 0;
    boost::unordered_map<int, std::vector<result_info>*> proc_out;
    writer_job j;
    while (!writerq->wait_and_pop(j))
    {
      if (j.molid == nwritten) {
        write_out(*j.results, *outfile, *outext, *gs->settings, *gs->wt,
            *outflex, *outfext, *gs->atomoutfile);
        nwritten++;
        delete j.results;
        for (boost::unordered_map<int, std::vector<result_info>*>::iterator i;
            (i = proc_out.find(nwritten)) != proc_out.end();)
            {
          write_out(*i->second, *outfile, *outext, *gs->settings,
              *gs->wt, *outflex, *outfext, *gs->atomoutfile);
          nwritten++;
          delete i->second;
        }
      }
      else {
        proc_out[j.molid] = j.results;
      }
    }
  } catch (file_error& e)
  {
    std::cerr << "\n\nError: could not open \"" << e.name.string()
        << "\" for " << (e.in ? "reading" : "writing") << ".\n";
  } catch (boost::filesystem::filesystem_error& e)
  {
    std::cerr << "\n\nFile system error: " << e.what() << '\n';
  } catch (usage_error& e)
  {
    std::cerr << "\n\nUsage error: " << e.what() << "\n";
  }
}

int main(int argc, char* argv[]) {
  using namespace boost::program_options;
  const std::string version_string =
      std::string("gnina ") + GIT_TAG + " " + GIT_BRANCH + ":"+GIT_REV + "   Built " __DATE__ ".";
  const std::string error_message =
      "\n\n\
Please report this error at https://github.com/gnina/gnina/issues\n"
          "Please remember to include the following in your problem report:\n\
    * the EXACT error message,\n\
    * your version of the program,\n\
    * the type of computer system you are running it on,\n\
	* all command line options,\n\
	* configuration file (if used),\n\
    * ligand file as provided to gnina,\n\
    * receptor file as provided to gnina,\n\
	* output file (if any),\n\
	* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";

  const std::string cite_message =
      "              _             \n"
          "             (_)            \n"
          "   __ _ _ __  _ _ __   __ _ \n"
          "  / _` | '_ \\| | '_ \\ / _` |\n"
          " | (_| | | | | | | | | (_| |\n"
          "  \\__, |_| |_|_|_| |_|\\__,_|\n"
          "   __/ |                    \n"
          "  |___/                     \n"
          "\n" + version_string +
          "\ngnina is based on smina and AutoDock Vina.\nPlease cite appropriately.\n";

  try
  {
    std::string rigid_name, flex_name, config_name, log_name, atom_name;
    std::vector<std::string> ligand_names;
    std::string out_name;
    std::string outf_name;
    std::string ligand_names_file;
    std::string atomconstants_file;
    std::string custom_file_name;
    std::string usergrid_file_name;
    std::string flex_res;
    double flex_dist = -1.0;
    fl center_x = 0, center_y = 0, center_z = 0, size_x = 0, size_y = 0,
        size_z = 0;
    fl autobox_add = 4;
    bool autobox_extend = true;
    std::string autobox_ligand;
    std::string flexdist_ligand;
    std::string builtin_scoring;
    int flex_limit = -1;
    int flex_max = -1;
    int nflex = -1;
    bool nflex_hard_limit = true; // TODO@RMeli: Use for defining "soft" flexmax

    // -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 0.05846
    fl weight_gauss1 = -0.035579;
    fl weight_gauss2 = -0.005156;
    fl weight_repulsion = 0.840245;
    fl weight_hydrophobic = -0.035069;
    fl weight_hydrogen = -0.587439;
    fl weight_rot = 0.05846;
    fl user_grid_lambda;
    bool help = false, help_hidden = false, version = false;
    bool quiet = false;
    bool accurate_line = false;
    bool simple_ascent = false;
    bool flex_hydrogens = false;
    bool print_terms = false;
    bool print_atom_types = false;
    bool add_hydrogens = true;
    bool strip_hydrogens = false;
    bool full_flex_output = false;

    CovOptions copt;

    user_settings settings;
    cnn_options& cnnopts = settings.cnnopts;

    minimization_params minparms;
    ApproxType approx = LinearApprox;
    fl approx_factor = 32;

    positional_options_description positional; // remains empty

    options_description inputs("Input");
    inputs.add_options()
    ("receptor,r", value<std::string>(&rigid_name),
        "rigid part of the receptor")
    ("flex", value<std::string>(&flex_name),
        "flexible side chains, if any (PDBQT)")
    ("ligand,l", value<std::vector<std::string> >(&ligand_names),
        "ligand(s)")
    ("flexres", value<std::string>(&flex_res),
        "flexible side chains specified by comma separated list of chain:resid")
    ("flexdist_ligand", value<std::string>(&flexdist_ligand),
        "Ligand to use for flexdist")
    ("flexdist", value<double>(&flex_dist),
        "set all side chains within specified distance to flexdist_ligand to flexible")
    ("flex_limit", value<int>(&flex_limit),
        "Hard limit for the number of flexible residues")
    ("flex_max", value<int>(&flex_max),
        "Retain at at most the closest flex_max flexible residues");

    //options_description search_area("Search area (required, except with --score_only)");
    options_description search_area("Search space (required)");
    search_area.add_options()
    ("center_x", value<fl>(&center_x), "X coordinate of the center")
    ("center_y", value<fl>(&center_y), "Y coordinate of the center")
    ("center_z", value<fl>(&center_z), "Z coordinate of the center")
    ("size_x", value<fl>(&size_x), "size in the X dimension (Angstroms)")
    ("size_y", value<fl>(&size_y), "size in the Y dimension (Angstroms)")
    ("size_z", value<fl>(&size_z), "size in the Z dimension (Angstroms)")
    ("autobox_ligand", value<std::string>(&autobox_ligand),
        "Ligand to use for autobox")
    ("autobox_add", value<fl>(&autobox_add),
        "Amount of buffer space to add to auto-generated box (default +4 on all six sides)")
    ("autobox_extend", value<bool>(&autobox_extend)->default_value(true),
            "Expand the autobox if needed to ensure the input conformation of the ligand being docked can freely rotate within the box.")
    ("no_lig", bool_switch(&settings.no_lig)->default_value(false),
        "no ligand; for sampling/minimizing flexible residues");

    options_description covalent("Covalent docking");
    covalent.add_options()
    ("covalent_rec_atom", value<std::string>(&copt.covalent_rec_atom), "Receptor atom ligand is covalently bound to.  Can be specified as chain:resnum:atom_name or as x,y,z Cartesian coordinates.")
    ("covalent_lig_atom_pattern", value<std::string>(&copt.covalent_lig_atom_pattern), "SMARTS expression for ligand atom that will covalently bind protein.")
    ("covalent_lig_atom_position", value<std::string>(&copt.covalent_lig_atom_position), "Optional.  Initial placement of covalently bonding ligand atom in x,y,z Cartesian coordinates.  If not specified, OpenBabel's GetNewBondVector function will be used to position ligand.")
    ("covalent_fix_lig_atom_position", bool_switch(&copt.covalent_fix_lig_atom_position), "If covalent_lig_atom_position is specified, fix the ligand atom to this position as opposed to using this position to define the initial structure.")
    ("covalent_bond_order", value<int>(&copt.bond_order)->default_value(1), "Bond order of covalent bond. Default 1.")
    ("covalent_optimize_lig", bool_switch(&copt.covalent_optimize_lig), "Optimize the covalent complex of ligand and residue using UFF. This will change bond angles and lengths of the ligand.");


    options_description outputs("Output");
    outputs.add_options()
    ("out,o", value<std::string>(&out_name),
        "output file name, format taken from file extension")
    ("out_flex", value<std::string>(&outf_name),
        "output file for flexible receptor residues")
    ("log", value<std::string>(&log_name), "optionally, write log file")
    ("atom_terms", value<std::string>(&atom_name),
        "optionally write per-atom interaction term values")
    ("atom_term_data",
        bool_switch(&settings.include_atom_info)->default_value(false),
        "embedded per-atom interaction terms in output sd data")
    ("pose_sort_order",value<pose_sort_order>(&settings.sort_order)->default_value(CNNscore),
        "How to sort docking results: CNNscore (default), CNNaffinity, Energy")
    ("full_flex_output", bool_switch(&full_flex_output)->default_value(false), "Output entire structure for out_flex, not just flexible residues.");


    options_description scoremin("Scoring and minimization options");
    scoremin.add_options()
    ("scoring", value<std::string>(&builtin_scoring),
        ("specify alternative built-in scoring function: "+builtin_scoring_functions.names(" ")).c_str())
    ("custom_scoring", value<std::string>(&custom_file_name),
        "custom scoring function file")
    ("custom_atoms", value<std::string>(&atomconstants_file),
        "custom atom type parameters file")
    ("score_only", bool_switch(&settings.score_only)->default_value(false),
        "score provided ligand pose")
    ("local_only", bool_switch(&settings.local_only)->default_value(false),
        "local search only using autobox (you probably want to use --minimize)")
    ("minimize", bool_switch(&settings.dominimize)->default_value(false),
        "energy minimization")
    ("randomize_only", bool_switch(&settings.randomize_only),
        "generate random poses, attempting to avoid clashes")
    ("num_mc_steps", value<int>(&settings.num_mc_steps),
        "fixed number of monte carlo steps to take in each chain")
    ("max_mc_steps", value<int>(&settings.max_mc_steps),
        "cap on number of monte carlo steps to take in each chain")
    ("num_mc_saved", value<int>(&settings.num_mc_saved),
            "number of top poses saved in each monte carlo chain")
    ("temperature", value<fl>(&settings.temperature),
            "temperature for metropolis accept criterion")
    ("minimize_iters",
        value<unsigned>(&minparms.maxiters)->default_value(0),
        "number iterations of steepest descent; default scales with rotors and usually isn't sufficient for convergence")
    ("accurate_line", bool_switch(&accurate_line),
        "use accurate line search")
    ("simple_ascent", bool_switch(&simple_ascent), "use simple gradient ascent")
    ("minimize_early_term", bool_switch(&minparms.early_term),
        "Stop minimization before convergence conditions are fully met.")
    ("minimize_single_full", bool_switch(&minparms.single_min),
        "During docking perform a single full minimization instead of a truncated pre-evaluate followed by a full.")
    ("approximation", value<ApproxType>(&approx),
        "approximation (linear, spline, or exact) to use")
    ("factor", value<fl>(&approx_factor),
        "approximation factor: higher results in a finer-grained approximation")
    ("force_cap", value<fl>(&settings.forcecap),
        "max allowed force; lower values more gently minimize clashing structures")
    ("user_grid", value<std::string>(&usergrid_file_name),
        "Autodock map file for user grid data based calculations")
    ("user_grid_lambda", value<fl>(&user_grid_lambda)->default_value(-1.0),
        "Scales user_grid and functional scoring")
    ("print_terms", bool_switch(&print_terms),
        "Print all available terms with default parameterizations")
    ("print_atom_types", bool_switch(&print_atom_types),
        "Print all available atom types");

    options_description hidden("Hidden options for internal testing");
    hidden.add_options()
    ("verbosity", value<int>(&settings.verbosity)->default_value(1),
        "Adjust the verbosity of the output, default: 1")
    ("flex_hydrogens", bool_switch(&flex_hydrogens),
        "Enable torsions affecting only hydrogens (e.g. OH groups). This is stupid but provides compatibility with Vina.")
    ("outputmin", value<int>(&minparms.outputframes),
        "output minout.sdf of minimization with provided amount of interpolation")
    ("cnn_gradient_check",
        bool_switch(&cnnopts.gradient_check)->default_value(false),
        "Perform internal checks on gradient.")
    ("gpu_docking", bool_switch(&settings.gpu_docking), "Turn on GPU acceleration for non-CNN scoring operations.");


    options_description cnn("Convolutional neural net (CNN) scoring");
    cnn.add_options()
    ("cnn_scoring",value<cnn_scoring_level>(&cnnopts.cnn_scoring)->default_value(CNNrescore),
                    "Amount of CNN scoring: none, rescore (default), refinement, metrorescore (metropolis+rescore), metrorefine (metropolis+refine), all")
    ("cnn", value<std::vector<std::string> >(&cnnopts.cnn_model_names)->multitoken(),
        ("built-in model to use, specify PREFIX_ensemble to evaluate an ensemble of models starting with PREFIX: " + builtin_cnn_models()).c_str())
    ("cnn_model", value<std::vector<std::string>>(&cnnopts.cnn_models)->multitoken(),
        "caffe cnn model file; if not specified a default model will be used")
    ("cnn_weights", value<std::vector<std::string>>(&cnnopts.cnn_weights)->multitoken(),
        "caffe cnn weights file (*.caffemodel); if not specified default weights (trained on the default model) will be used")
    ("cnn_resolution", value<fl>(&cnnopts.resolution)->default_value(0.5),
        "resolution of grids, don't change unless you really know what you are doing")
    ("cnn_rotation", value<unsigned>(&cnnopts.cnn_rotations)->default_value(0),
        "evaluate multiple rotations of pose (max 24)")
    ("cnn_update_min_frame", value<bool>(&cnnopts.move_minimize_frame)->default_value(true),
        "During minimization, recenter coordinate frame as ligand moves")
    ("cnn_freeze_receptor", bool_switch(&cnnopts.fix_receptor),
        "Don't move the receptor with respect to a fixed coordinate system")
    ("cnn_mix_emp_force", bool_switch(&cnnopts.mix_emp_force)->default_value(false),
        "Merge CNN and empirical minus forces")
    ("cnn_mix_emp_energy", bool_switch(&cnnopts.mix_emp_energy)->default_value(false),
        "Merge CNN and empirical energy")
    ("cnn_empirical_weight", value<fl>(&cnnopts.empirical_weight)->default_value(1.0),
        "Weight for scaling and merging empirical force and energy ")
    ("cnn_outputdx", bool_switch(&cnnopts.outputdx),
        "Dump .dx files of atom grid gradient.")
    ("cnn_outputxyz", bool_switch(&cnnopts.outputxyz),
        "Dump .xyz files of atom gradient.")
    ("cnn_xyzprefix",
        value<std::string>(&cnnopts.xyzprefix)->default_value("gradient"),
        "Prefix for atom gradient .xyz files")
    ("cnn_center_x", value<fl>(&cnnopts.cnn_center[0]),
        "X coordinate of the CNN center")
    ("cnn_center_y", value<fl>(&cnnopts.cnn_center[1]),
        "Y coordinate of the CNN center")
    ("cnn_center_z", value<fl>(&cnnopts.cnn_center[2]),
        "Z coordinate of the CNN center")
    ("cnn_verbose", bool_switch(&cnnopts.verbose),
        "Enable verbose output for CNN debugging");

    options_description misc("Misc (optional)");
    misc.add_options()
    ("cpu", value<int>(&settings.cpu),
        "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
    ("seed", value<int>(&settings.seed), "explicit random seed")
    ("exhaustiveness",
        value<int>(&settings.exhaustiveness)->default_value(8),
        "exhaustiveness of the global search (roughly proportional to time)")
    ("num_modes", value<sz>(&settings.num_modes)->default_value(9),
        "maximum number of binding modes to generate")
    ("min_rmsd_filter", value<fl>(&settings.out_min_rmsd)->default_value(1.0),
        "rmsd value used to filter final poses to remove redundancy")
    ("quiet,q", bool_switch(&quiet), "Suppress output messages")
    ("addH", value<bool>(&add_hydrogens),
        "automatically add hydrogens in ligands (on by default)")
    ("stripH", value<bool>(&strip_hydrogens),
        "remove hydrogens from molecule _after_ performing atom typing for efficiency (off by default)")
    ("device", value<int>(&settings.device)->default_value(0),
        "GPU device to use")
    ("no_gpu", bool_switch(&settings.no_gpu), "Disable GPU acceleration, even if available.");


    options_description config("Configuration file (optional)");
    config.add_options()("config", value<std::string>(&config_name),
        "the above options can be put here");
    options_description info("Information (optional)");
    info.add_options()
    ("help", bool_switch(&help), "display usage summary")
    ("help_hidden", bool_switch(&help_hidden),
        "display usage summary with hidden options")
    ("version", bool_switch(&version), "display program version");

    options_description desc, desc_simple;
    desc.add(inputs).add(search_area).add(covalent).add(outputs).add(scoremin).add(cnn).
        add(hidden).add(misc).add(config).add(info);
    desc_simple.add(inputs).add(search_area).add(covalent).add(scoremin).add(cnn).
        add(outputs).add(misc).add(config).add(info);

    variables_map vm;
    try
    {
      store(
          command_line_parser(argc, argv).options(desc)
              .style(
              command_line_style::default_style
                  ^ command_line_style::allow_guessing)
              .positional(positional).run(), vm);
      notify(vm);
    } catch (boost::program_options::error& e)
    {
      std::cerr << "Command line parse error: " << e.what() << '\n'
          << "\nCorrect usage:\n" << desc_simple << '\n';
      return 1;
    }
    if (vm.count("config"))
        {
      try
      {
        ifile config_stream(config_name);
        store(parse_config_file(config_stream, desc), vm);
        notify(vm);
      } catch (boost::program_options::error& e)
      {
        std::cerr << "Configuration file parse error: " << e.what()
            << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
        return 1;
      }
    }
    if (help)
    {
      std::cout << desc_simple << '\n';
      return 0;
    }
    if (help_hidden)
    {
      std::cout << desc << '\n';
      return 0;
    }
    if (version)
    {
      std::cout << version_string << '\n';
      return 0;
    }

    tee log(quiet);
    if (vm.count("log") > 0)
      log.init(log_name);

    if (!atomconstants_file.empty())
      setup_atomconstants_from_file(atomconstants_file);

    if (print_terms)
    {
      custom_terms t;
      t.print_available_terms(std::cout);
      return 0;
    }

    if (print_atom_types)
    {
      print_atom_info(std::cout);
      return 0;
    }

    FLAGS_minloglevel = google::GLOG_ERROR; //don't spit out info messages
    // Google logging.
    ::google::InitGoogleLogging(argv[0]);
    // Provide a backtrace on segfault.
    ::google::InstallFailureSignalHandler();

#if (OB_VERSION > OB_VERSION_CHECK(2, 3, 2))
    OpenBabel::OBPlugin::LoadAllPlugins(); //for some reason loading on demand can be slow
#endif
    cnnopts.seed = settings.seed;

    OpenBabel::vector3 dummy; //openbabel uses system rand initialized with time seed
    dummy.randomUnitVector(); //this setups up the obrandom object with time
    srand(settings.seed); //so now it is safe(?) to set the system seed

    set_fixed_rotable_hydrogens(!flex_hydrogens);

    if (settings.dominimize) //set default settings for minimization
    {
      if (!vm.count("force_cap"))
        settings.forcecap = 10; //nice and soft

      if (minparms.maxiters == 0)
        minparms.maxiters = 10000; //will presumably converge
      settings.local_only = true;
      minparms.type = minimization_params::BFGSAccurateLineSearch;

      if (!vm.count("approximation"))
        approx = SplineApprox; //use high accuracy approximation for --minimize
      if (!vm.count("factor"))
        approx_factor = 10;
    }

    //check for GPU
    int cudaerror = cudaDeviceReset();
    if(cudaerror != cudaSuccess) {
      caffe::Caffe::set_cudnn(false); //if cudnn is on, won't fallback to cpu
    } else if(settings.no_gpu) {
      caffe::Caffe::set_cudnn(false);
    } else {
      //the following reserve a lot of memory, do we really need it??
      //cudaDeviceSetLimit(cudaLimitStackSize, 5120);
      caffe::Caffe::SetDevice(settings.device);
      caffe::Caffe::set_mode(caffe::Caffe::GPU);      
    }

    if (accurate_line)
    {
      minparms.type = minimization_params::BFGSAccurateLineSearch;
    }

    if (simple_ascent)
    {
      minparms.type = minimization_params::Simple;
    }

    bool search_box_needed = !(settings.score_only || settings.local_only); // randomize_only and local_only still need the search space; dkoes - for local get box from ligand
    bool output_produced = !settings.score_only;
    bool receptor_needed = !settings.randomize_only;


    if (cnnopts.cnn_scoring == CNNall) {
      cnnopts.move_minimize_frame = true;
    }

    if(cnnopts.cnn_scoring == CNNnone) {
      settings.sort_order = Energy;
    }

    if (receptor_needed)
    {
      if (vm.count("receptor") <= 0)
          {
        std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n"
            << desc_simple << '\n';
        return 1;
      }
    }

    if (ligand_names.size() == 0) {
      if (!settings.no_lig)
      {
        std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n"
            << desc_simple << '\n';
        return 1;
      }
      else //put in "fake" ligand
      {
        ligand_names.push_back("");
      }
    }
    else if (settings.no_lig) //ligand specified with no_lig
    {
      std::cerr << "Ligand specified with --no_lig.\n"
          << "\nCorrect usage:\n"
          << desc_simple << '\n';
      return 1;
    } else {
      for(const auto& lname : ligand_names) {
        if(!boost::filesystem::exists(lname)) {
          throw file_error(lname, true);
        }
      }
    }

    if (settings.exhaustiveness < 1)
      throw usage_error("exhaustiveness must be 1 or greater");
    if (settings.num_modes < 1)
      throw usage_error("num_modes must be 1 or greater");

    boost::optional<std::string> flex_name_opt;
    if (vm.count("flex"))
      flex_name_opt = flex_name;

    if (vm.count("flex") && !vm.count("receptor"))
      throw usage_error(
          "Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

    if(flex_limit > -1 && flex_max > -1){
      throw usage_error(
          "--flex_lim and --flex_max can't be used together.");
    }
    else if(flex_limit > -1){
      nflex = flex_limit;
      nflex_hard_limit = true;
    }
    else if(flex_max > -1){
      nflex = flex_max;
      nflex_hard_limit = false;
    }

    std::ofstream atomoutfile;
    if (vm.count("atom_terms") > 0)
      atomoutfile.open(atom_name.c_str());

    //output banner
    log << cite_message << '\n';

    if(cudaerror != cudaSuccess) {
      log << "WARNING: No GPU detected. CNN scoring will be slow.\n"
          "Recommend running with single model (--cnn crossdock_default2018)\n"
          "or without cnn scoring (--cnn_scoring=none).\n\n";
    }
    log << "Commandline:";
    for(unsigned i = 0; i < argc; i++) {
      log << " " << argv[i];
    }
    log << "\n";

    FlexInfo finfo(flex_res, flex_dist, flexdist_ligand, nflex, nflex_hard_limit, full_flex_output, log);
    CovInfo cinfo(copt, log);
    // dkoes - parse in receptor once
    MolGetter mols(rigid_name, flex_name, finfo, cinfo, add_hydrogens, strip_hydrogens, log);

    if (autobox_ligand.length() > 0) {
      setup_autobox(mols.getInitModel(),autobox_ligand, autobox_add,
          center_x, center_y, center_z, size_x, size_y, size_z);
    }

    if (search_box_needed && autobox_ligand.length() == 0)  {
      options_occurrence oo = get_occurrence(vm, search_area);
      if (!oo.all)
      {
        check_occurrence(vm, search_area);
        std::cerr << "\nCorrect usage:\n" << desc_simple << std::endl;
        return 1;
      }
      if (size_x <= 0 || size_y <= 0 || size_z <= 0)
        throw usage_error("Search space dimensions should be positive");
    }

    if (flex_dist > 0 && flexdist_ligand.size() == 0) {
      throw usage_error("Must specify flexdist_ligand with flex_dist");
    }

    if(nflex > 0 && flex_res.size() > 0){
      log << "WARNING: --flex_limit and --flexmax ignored with --flexres\n\n";
    }

    grid_dims gd; // n's = 0 via default c'tor
    grid_dims user_gd;
    grid user_grid;

    flv weights;

    //dkoes, set the scoring function
    custom_terms t;
    if (user_grid_lambda != -1.0)
    {
      t.set_scaling_factor(user_grid_lambda);
    }
    if (custom_file_name.size() > 0)
    {
      ifile custom_file(custom_file_name);
      t.add_terms_from_file(custom_file);
    }
    else if (builtin_scoring.size() > 0)
    {
      if (!builtin_scoring_functions.set(t, builtin_scoring))
      {
        throw usage_error(
            "Invalid built-in scoring function: " + builtin_scoring
                + ". Options are:\n" + builtin_scoring_functions.names("\n"));
      }
    }
    else
    {
      t.add("gauss(o=0,_w=0.5,_c=8)", -0.035579);
      t.add("gauss(o=3,_w=2,_c=8)", -0.005156);
      t.add("repulsion(o=0,_c=8)", 0.840245);
      t.add("hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069);
      t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439);
      t.add("num_tors_div", 5 * 0.05846 / 0.1 - 1);
    }

    if (settings.verbosity > 1)
      log << std::setw(12) << std::left << "Weights" << " Terms\n" << t << "\n";

    // Print out flexible residues 
    if(finfo.has_content()){
      finfo.print_flex();
    }
    
    // Print information about flexible residues use
    if (finfo.has_content() && cnnopts.cnn_scoring != CNNnone) {
      if(!cnnopts.fix_receptor && settings.verbosity > 1)
          log << "Receptor position and orientation are frozen.\n";
      cnnopts.fix_receptor = true; // Fix receptor position and orientation
    }

    if (usergrid_file_name.size() > 0) {
      ifile user_in(usergrid_file_name);
      fl ug_scaling_factor = 1.0;
      if (user_grid_lambda != -1.0)
          {
        ug_scaling_factor = 1 - user_grid_lambda;
      }
      setup_user_gd(user_gd, user_in);
      user_grid.init(user_gd, user_in, ug_scaling_factor); //initialize user grid
    }

    if (search_box_needed) {
      setup_grid_dims(center_x, center_y, center_z, size_x, size_y, size_z, gd);
    }

    if (settings.verbosity > 1) {
      log << "Using search box with center " << center_x << "," << center_y << "," << center_z << \
          " and size " << size_x << "," << size_y << "," << size_z << "\n";
    }

    if (vm.count("cpu") == 0) {
      settings.cpu = boost::thread::hardware_concurrency();
      if (settings.verbosity > 1) {
        if (settings.cpu > 0)
          log << "Detected " << settings.cpu << " CPU"
              << ((settings.cpu > 1) ? "s" : "") << '\n';
        else
          log << "Could not detect the number of CPUs, using 1\n";
      }
    }
    if (settings.cpu < 1)
      settings.cpu = 1;
    if (settings.verbosity > 1 && settings.exhaustiveness < settings.cpu)
      log  << "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";

    if (settings.verbosity <= 1) {
      OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obError);
    }
    //dkoes, hoist precalculation outside of loop
    weighted_terms wt(&t, t.weights());

    boost::shared_ptr<precalculate> prec;

    if (settings.gpu_docking || approx == GPU)  { //don't get a choice
      prec = boost::shared_ptr<precalculate>(
          new precalculate_gpu(wt, approx_factor));
    }
    else if (approx == SplineApprox)
      prec = boost::shared_ptr<precalculate>(
          new precalculate_splines(wt, approx_factor));
    else if (approx == LinearApprox)
      prec = boost::shared_ptr<precalculate>(
          new precalculate_linear(wt, approx_factor));
    else if (approx == Exact)
      prec = boost::shared_ptr<precalculate>(
          new precalculate_exact(wt));

    //setup single outfile
    using namespace OpenBabel;
    ozfile outfile;
    std::string outext;
    if (out_name.length() > 0) {
      outext = outfile.open(out_name);
    }

    ozfile outflex;
    std::string outfext;
    if (outf_name.length() > 0)
    {
      outfext = outflex.open(outf_name);
    }

    if (settings.score_only) //output header
    {
      std::vector<std::string> enabled_names = t.get_names(true);
      log << "## Name";
      VINA_FOR_IN(i, enabled_names)
      {
        log << " " << enabled_names[i];
      }
      for (unsigned i = 0, n = t.conf_independent_terms.size(); i < n;
          i++)
          {
        log << " " << t.conf_independent_terms[i].name;
      }
      log << "\n";
    }

    job_queue<worker_job> wrkq;
    job_queue<writer_job> writerq;
    int nligs = 0;
    size_t nthreads = settings.cpu;
    global_state gs(&settings, prec, &minparms, &wt, &user_grid,
        &log, &atomoutfile, cnnopts);
    boost::thread_group worker_threads;
    boost::timer::cpu_timer time;
    CNNScorer cnn_scorer(cnnopts); //shared network

    if (!settings.local_only)
      nthreads = 1; //docking is multithreaded already, don't add additional parallelism other than pipeline

    //launch worker threads to process ligands in the work queue
    for (int i = 0; i < nthreads; i++) {
      worker_threads.create_thread(boost::bind(threads_at_work, &wrkq,
          &writerq, &gs, &mols, &nligs, cnn_scorer));
    }

    //launch writer thread to write results wherever they go
    boost::thread writer_thread(thread_a_writing, &writerq, &gs, &outfile,
        &outext, &outflex, &outfext, &nligs);

    try {
      //loop over input ligands, adding them to the work queue
      for (unsigned l = 0, nl = ligand_names.size(); l < nl; l++) {
        doing(settings.verbosity, "Reading input", log);
        const std::string ligand_name = ligand_names[l];
        mols.setInputFile(ligand_name);

        unsigned i = 0;

        for (;;)  {
          model* m = new model;

          if (!mols.readMoleculeIntoModel(*m))  {
            delete m;
            break;
          }
          m->set_pose_num(i);
          m->gdata.device_on = settings.gpu_docking;
          m->gdata.device_id = settings.device;

          grid_dims gdbox(gd);
          if (settings.local_only)
          {
            gdbox = m->movable_atoms_box(autobox_add, box_granularity);
            bool skip = false;
            for(unsigned pos = 0; pos < 3; pos++) {
              if(gdbox.elems[pos].n*box_granularity > 100) {
                //we will run out of memory if the grid is too large
                log << "WARNING: Ligand " << i << " in " << ligand_name << " has an extent greater than 100A. Skipping.\n";
                skip = true;
                break;
              }
            }
            if(skip)
              continue;
          } else if(autobox_extend && !settings.no_lig) {
            //make sure every dimension is large enough for the ligand to fit
            fl maxdim = m->max_span(0);
            setup_grid_dims(center_x, center_y, center_z,
                size_x > maxdim ? size_x : maxdim,
                size_y > maxdim ? size_y : maxdim,
                size_z > maxdim ? size_z : maxdim, gdbox);
          }

          done(settings.verbosity, log);
          std::vector<result_info>* results = new std::vector<result_info>();
          worker_job j(i, m, results, gdbox);
          wrkq.push(j);

          i++;
          if (settings.no_lig)
            break;
        }
      }
    } catch (...)
    {
      //clean up threads before passing along exception
      wrkq.close(nthreads);
      worker_threads.join_all();
      writerq.close(1);
      writer_thread.join();
      cudaDeviceSynchronize();
      throw;
    }

    //join all the threads when their work is done
    wrkq.close(nthreads);
    worker_threads.join_all();
    writerq.close(1);
    writer_thread.join();

    sz free_byte = 0, total_byte = 0;
    if(settings.verbosity > 1 && cudaMemGetInfo( &free_byte, &total_byte ) == cudaSuccess) {
      double free_db = (double)free_byte ;
      double total_db = (double)total_byte ;
      double used_db = total_db - free_db ;
      log << "GPU memory usage: " << int(used_db/1024.0/1024.0) << " MB" << "\n";
    }

    cudaDeviceSynchronize();

    //std::cout << "Loop time " << time.elapsed().wall / 1000000000.0 << "\n";

  } catch (file_error& e)
  {
    std::cerr << "\n\nError: could not open \"" << e.name.string()
        << "\" for " << (e.in ? "reading" : "writing") << ".\n";
    return 1;
  } catch (boost::filesystem::filesystem_error& e)
  {
    std::cerr << "\n\nFile system error: " << e.what() << '\n';
    return 1;
  } catch (usage_error& e)
  {
    std::cerr << "\n\nUsage error: " << e.what() << "\n";
    return 1;
  } catch (parse_error& e)
  {
    std::cerr << "\n\nParse error on line " << e.line << " in file \""
        << e.file.string() << "\": " << e.reason << '\n';
    return 1;
  } catch (std::bad_alloc&)
  {
    std::cerr << "\n\nError: insufficient memory!\n";
    return 1;
  } catch (scoring_function_error e)
  {
    std::cerr << "\n\nError with scoring function specification.\n";
    std::cerr << e.msg << "[" << e.name << "]\n";
    return 1;
  } catch(std::runtime_error e)
  {
    std::cerr << "\nRuntime Error\n";
    std::cerr << e.what() << "\n";
    sz free_byte = 0, total_byte = 0;
    cudaMemGetInfo( &free_byte, &total_byte ) ;

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;
    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
             used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

    return 1;
  }

// Errors that shouldn't happen:

  catch (std::exception& e)
  {
    std::cerr << "\n\nAn error occurred: " << e.what() << ". "
        << error_message;
    return 1;
  } catch (internal_error& e)
  {
    std::cerr << "\n\nAn internal error occurred in " << e.file << "("
        << e.line << "). " << error_message;
    return 1;
  } catch (...)
  {
    std::cerr << "\n\nAn unknown error occurred. " << error_message;
    return 1;
  }
}
