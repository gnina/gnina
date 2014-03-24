#include <iostream>
#include <string>
#include <exception>
#include <vector> // ligand paths#include <cmath> // for ceila#include <boost/program_options.hpp>#include <boost/filesystem/fstream.hpp>#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/convenience.hpp> // filesystem::basename#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?#include <boost/lexical_cast.hpp>#include "parse_pdbqt.h"#include "parallel_mc.h"
#include "file.h"
#include "cache.h"
#include "non_cache.h"
#include "naive_non_cache.h"
#include "non_cache_gpu.h"
#include "parse_error.h"
#include "everything.h"
#include "weighted_terms.h"
#include "quasi_newton.h"
#include "tee.h"
#include "custom_terms.h"
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
#include <boost/shared_ptr.hpp>
#include "coords.h"
#include "obmolopener.h"
#include "gpucode.h"
#include "precalculate_gpu.h"
#include <boost/timer.hpp>
#include <boost/algorithm/string.hpp>
#include "array3d.h"
#include "grid.h"

using namespace boost::iostreams;
using boost::filesystem::path;

struct resultInfo
{
	fl energy;
	fl rmsd;
	std::string mol;
	std::string atominfo;

	resultInfo() :
			energy(0), rmsd(-1)
	{
	}
	resultInfo(fl e, fl r, const std::string& m) :
			energy(e), rmsd(r), mol(m)
	{
	}

	//write a table (w/header) of per atom values to out
	void writeAtomValues(std::ostream& out, const weighted_terms *wt) const
			{
		const terms *t = wt->unweighted_terms();
		out << "atomid el pos";
		std::vector<std::string> names = t->get_names(true);
		for (unsigned j = 0, m = names.size(); j < m; j++)
		{
			out << " " << names[j];
		}
		out << "\n";
		out << atominfo;
	}

	//computes per-atom term values and formats them into the atominfo string
	void setAtomValues(const model& m, const weighted_terms *wt)
	{
		std::vector<flv> values;
		const terms *t = wt->unweighted_terms();
		t->evale_robust(m, values);
		std::stringstream str;
		vecv coords = m.get_ligand_coords();
		assert(values.size() == coords.size());
		for (unsigned i = 0, n = values.size(); i < n; i++)
		{
			assert(values[i].size() == t->size());
			str << m.ligand_atom_str(i) << " ";
			coords[i].print(str);
			for (unsigned j = 0, m = values[i].size(); j < m; j++)
			{
				str << " " << values[i][j] * wt->weight(j);
			}
			str << "\n";
		}
		str << "END\n";
		atominfo = str.str();
	}
};

path make_path(const std::string& str)
{
	return path(str);
}

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

fl do_randomization(model& m, std::ostream& out, const vec& corner1,
		const vec& corner2, int seed, int verbosity, tee& log)
{
	conf init_conf = m.get_initial_conf();
	rng generator(static_cast<rng::result_type>(seed));
	if (verbosity > 1)
	{
		log << "Using random seed: " << seed;
		log.endl();
	}
	const sz attempts = 10000;
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
		}
	}
	m.set(best_conf);
	if (verbosity > 1)
	{
		log << "Clash penalty: " << best_clash_penalty; // FIXME rm?
		log.endl();
	}
	m.write_structure(out);
	return best_clash_penalty;
}

void refine_structure(model& m, const precalculate& prec, non_cache& nc,
		output_type& out, const vec& cap, const minimization_params& minparm,
		grid& user_grid)
{
	change g(m.get_size());
	quasi_newton quasi_newton_par(minparm);
	const fl slope_orig = nc.getSlope();
	//try 5 times to get ligand into box
	VINA_FOR(p, 5)
	{
		nc.setSlope(100 * std::pow(10.0, 2.0 * p));
		quasi_newton_par(m, prec, nc, out, g, cap, user_grid); //quasi_newton operator
		m.set(out.c); // just to be sure
		if (nc.within(m))
			break;
	}
	out.coords = m.get_heavy_atom_movable_coords();
	if (!nc.within(m))
		out.e = max_fl;
	nc.setSlope(slope_orig);
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
	VINA_FOR_IN(i, in)
		add_to_output_container(tmp, in[i], min_rmsd, in.size());
	return tmp;
}

//dkoes - return all energies and rmsds to original conf with result
void do_search(model& m, const boost::optional<model>& ref,
		const weighted_terms& sf, const precalculate& prec, const igrid& ig,
		non_cache& nc, // nc.slope is changed
		const vec& corner1, const vec& corner2,
		const parallel_mc& par, fl energy_range, sz num_modes, int seed,
		int verbosity, bool score_only, bool local_only, bool compute_atominfo,
		fl out_min_rmsd, tee& log,
		const terms *t, grid& user_grid, std::vector<resultInfo>& results)
{
	boost::timer time;

	precalculate_exact exact_prec(sf); //use exact computations for final score
	conf_size s = m.get_size();
	conf c = m.get_initial_conf();
	fl e = max_fl;
	fl rmsd = 0;
	const vec authentic_v(1000, 1000, 1000);
	if (score_only)
	{
		fl intramolecular_energy = m.eval_intramolecular(exact_prec,
				authentic_v, c);
		naive_non_cache nnc(&exact_prec); // for out of grid issues
		e = m.eval_adjusted(sf, exact_prec, nnc, authentic_v, c,
				intramolecular_energy, user_grid);

		log << "Affinity: " << std::fixed << std::setprecision(5) << e
				<< " (kcal/mol)";
		log.endl();
		std::vector<flv> atominfo;
		flv term_values = t->evale_robust(m);
		log << "Intramolecular energy: " << std::fixed << std::setprecision(5)
				<< intramolecular_energy << "\n";

		log
		<< "Term values, before weighting:\n";
		log << std::setprecision(5);
		log << "## " << boost::replace_all_copy(m.get_name()," ","_");

		VINA_FOR_IN(i, term_values)
		{
			log << ' ' << term_values[i];
		}

		conf_independent_inputs in(m);
		const flv nonweight(1, 1.0);
		for (unsigned i = 0, n = t->conf_independent_terms.size(); i < n; i++)
		{
			flv::const_iterator pos = nonweight.begin();
			log << " " <<
					t->conf_independent_terms[i].eval(in, (fl) 0.0, pos);
		}
		log << '\n';

		results.push_back(resultInfo(e, -1, ""));
		if (compute_atominfo)
			results.back().setAtomValues(m, &sf);
	}
	else if (local_only)
	{
		vecv origcoords = m.get_heavy_atom_movable_coords();
		output_type out(c, e);
		doing(verbosity, "Performing local search", log);
		refine_structure(m, prec, nc, out, authentic_v, par.mc.ssd_par.minparm,
				user_grid);
		done(verbosity, log);

		//be as exact as possible for final score
		naive_non_cache nnc(&exact_prec); // for out of grid issues

		fl intramolecular_energy = m.eval_intramolecular(exact_prec,
				authentic_v,
				out.c);
		e = m.eval_adjusted(sf, exact_prec, nnc, authentic_v, out.c,
				intramolecular_energy, user_grid);

		vecv newcoords = m.get_heavy_atom_movable_coords();
		assert(newcoords.size() == origcoords.size());
		for (unsigned i = 0, n = newcoords.size(); i < n; i++)
		{
			rmsd += (newcoords[i] - origcoords[i]).norm_sqr();
		}
		rmsd /= newcoords.size();
		rmsd = sqrt(rmsd);
		log << "Affinity: " << std::fixed << std::setprecision(5) << e << "  "
				<< intramolecular_energy
				<< " (kcal/mol)\nRMSD: " << rmsd;
		;
		log.endl();
		if (!nc.within(m))
			log
			<< "WARNING: not all movable atoms are within the search space\n";

		output_container out_cont;
		out_cont.push_back(new output_type(out));
		std::stringstream str;
		write_all_output(m, out_cont, 1, str);
		done(verbosity, log);
		results.push_back(resultInfo(e, rmsd, str.str()));
		if (compute_atominfo)
			results.back().setAtomValues(m, &sf);
	}
	else
	{
		rng generator(static_cast<rng::result_type>(seed));
		log << "Using random seed: " << seed;
		log.endl();
		output_container out_cont;
		doing(verbosity, "Performing search", log);
		par(m, out_cont, prec, ig, corner1, corner2, generator, user_grid);
		done(verbosity, log);
		doing(verbosity, "Refining results", log);
		VINA_FOR_IN(i, out_cont)
			refine_structure(m, prec, nc, out_cont[i], authentic_v,
					par.mc.ssd_par.minparm, user_grid);

		if (!out_cont.empty())
		{
			out_cont.sort();
			const fl best_mode_intramolecular_energy = m.eval_intramolecular(
					prec, authentic_v, out_cont[0].c);

			VINA_FOR_IN(i, out_cont)
				if (not_max(out_cont[i].e))
					out_cont[i].e = m.eval_adjusted(sf, prec, nc, authentic_v,
							out_cont[i].c, best_mode_intramolecular_energy,
							user_grid);
			// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
			out_cont.sort();
		}
		out_cont = remove_redundant(out_cont, out_min_rmsd);

		done(verbosity, log);

		log.setf(std::ios::fixed, std::ios::floatfield);
		log.setf(std::ios::showpoint);
		log << '\n';
		log << "mode |   affinity | dist from best mode\n";
		log << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
		log << "-----+------------+----------+----------\n";

		model best_mode_model = m;
		if (!out_cont.empty())
			best_mode_model.set(out_cont.front().c);

		sz how_many = 0;
		VINA_FOR_IN(i, out_cont)
		{
			if (how_many >= num_modes || !not_max(out_cont[i].e)
					|| out_cont[i].e > out_cont[0].e + energy_range)
				break; // check energy_range sanity FIXME
			++how_many;
			log << std::setw(4) << i + 1 << "    " << std::setw(9)
					<< std::setprecision(1) << out_cont[i].e; // intermolecular_energies[i];
			m.set(out_cont[i].c);
			const model& r = ref ? ref.get() : best_mode_model;
			const fl lb = m.rmsd_lower_bound(r);
			const fl ub = m.rmsd_upper_bound(r);
			log << "  " << std::setw(9) << std::setprecision(3) << lb << "  "
					<< std::setw(9) << std::setprecision(3) << ub; // FIXME need user-readable error messages in case of failures

			log.endl();

			//dkoes - setup resultInfo
			std::stringstream str;
			m.write_model(str, i + 1, "");
			results.push_back(resultInfo(out_cont[i].e, -1, str.str()));
			if (compute_atominfo)
				results.back().setAtomValues(m, &sf);

		}
		done(verbosity, log);

		if (how_many < 1)
		{
			log
					<< "WARNING: Could not find any conformations completely within the search space.\n"
					<< "WARNING: Check that it is large enough for all movable atoms, including those in the flexible side chains.";
			log.endl();
		}
	}
	//std::cout << "Refine time " << time.elapsed() << "\n";
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

void main_procedure(model& m, precalculate& prec,
		const boost::optional<model>& ref, // m is non-const (FIXME?)
		bool score_only, bool local_only,
		bool randomize_only, bool no_cache, bool compute_atominfo, bool gpu_on,
		const grid_dims& gd,
		int exhaustiveness, minimization_params minparm,
		const weighted_terms& wt, int cpu, int seed,
		int verbosity, sz num_modes, fl energy_range, fl out_min_rmsd, tee& log,
		std::vector<resultInfo>& results, grid& user_grid)
{
	doing(verbosity, "Setting up the scoring function", log);

	done(verbosity, log);

	vec corner1(gd[0].begin, gd[1].begin, gd[2].begin);
	vec corner2(gd[0].end, gd[1].end, gd[2].end);

	parallel_mc par;
	sz heuristic = m.num_movable_atoms()
			+ 10 * m.get_size().num_degrees_of_freedom();
	par.mc.num_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
	par.mc.ssd_par.evals = unsigned((25 + m.num_movable_atoms()) / 3);
	if (minparm.maxiters == 0)
		minparm.maxiters = par.mc.ssd_par.evals;
	par.mc.ssd_par.minparm = minparm;
	par.mc.min_rmsd = 1.0;
	par.mc.num_saved_mins = num_modes > 20 ? num_modes : 20; //dkoes, support more than 20
	par.mc.hunt_cap = vec(10, 10, 10);
	par.num_tasks = exhaustiveness;
	par.num_threads = cpu;
	par.display_progress = true;

	/*
	 std::cout << score_only << "\n";
	 std::cout << local_only << "\n";
	 std::cout << no_cache << "\n";
	 */

	szv_grid_cache gridcache(m, prec.cutoff_sqr());
	const fl slope = 1e6; // FIXME: too large? used to be 100
	if (randomize_only)
	{
		std::stringstream str;
		fl e = do_randomization(m, str, corner1, corner2, seed, verbosity, log);
		results.push_back(resultInfo(e, 0, str.str()));
		return;
	}
	else
	{

		non_cache *nc = NULL;
#ifdef SMINA_GPU
		if(gpu_on)
		{
			precalculate_gpu *gprec = dynamic_cast<precalculate_gpu*>(&prec);
			if(!gprec) abort();
			nc = new non_cache_gpu(gridcache, gd, gprec, slope);
		}
		else
#endif
		{
			nc = new non_cache(gridcache, gd, &prec, slope);
		}
		if (no_cache)
		{
			do_search(m, ref, wt, prec, *nc, *nc, corner1, corner2, par,
					energy_range, num_modes, seed, verbosity, score_only,
					local_only, compute_atominfo, out_min_rmsd, log,
					wt.unweighted_terms(), user_grid,
					results);
		}
		else
		{
			bool cache_needed = !(score_only || randomize_only || local_only);
			if (cache_needed)
				doing(verbosity, "Analyzing the binding site", log);
			cache c("scoring_function_version001", gd, slope);
			if (cache_needed)
			{
				std::vector<smt> atom_types_needed;
				m.get_movable_atom_types(atom_types_needed);
				c.populate(m, prec, atom_types_needed, user_grid);
			}
			if (cache_needed)
				done(verbosity, log);
			do_search(m, ref, wt, prec, c, *nc, corner1, corner2, par,
					energy_range, num_modes, seed, verbosity, score_only,
					local_only, compute_atominfo, out_min_rmsd, log,
					wt.unweighted_terms(), user_grid, results);
		}
		delete nc;
	}
}

struct usage_error: public std::runtime_error
{
	usage_error(const std::string& message) :
			std::runtime_error(message)
	{
	}
};

struct options_occurrence
{
	bool some;
	bool all;
	options_occurrence() :
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

//generate a box around the provided ligand padded by autobox_add
//if centers are always overwritten, but if sizes are non zero they are preserved
void setup_autobox(const std::string& autobox_ligand, fl autobox_add,
		fl& center_x,
		fl& center_y, fl& center_z, fl& size_x, fl& size_y, fl& size_z)
{
	using namespace OpenBabel;
	obmol_opener opener;

	//a ligand file can be provided from which to autobox
	OBConversion conv;
	opener.openForInput(conv, autobox_ligand);

	OBMol mol;
	center_x = center_y = center_z = 0;
	fl min_x = HUGE_VAL, min_y = HUGE_VAL, min_z = HUGE_VAL;
	fl max_x = -HUGE_VAL, max_y = -HUGE_VAL, max_z = -HUGE_VAL;
	fl num = 0;
	while (conv.Read(&mol)) //openbabel divides separate residues into multiple molecules
	{
		unsigned n;

		FOR_ATOMS_OF_MOL(a, mol)
				{
			center_x += a->x();
			center_y += a->y();
			center_z += a->z();
			min_x = std::min(min_x, (fl) a->x());
			min_y = std::min(min_y, (fl) a->y());
			min_z = std::min(min_z, (fl) a->z());
			max_x = std::max(max_x, (fl) a->x());
			max_y = std::max(max_y, (fl) a->y());
			max_z = std::max(max_z, (fl) a->z());
			num++;
		}
	}
	center_x /= num;
	center_y /= num;
	center_z /= num;
	if (size_x == 0)
		size_x = (max_x - min_x) + autobox_add;
	if (size_y == 0)
		size_y = (max_y - min_y) + autobox_add;
	if (size_z == 0)
		size_z = (max_z - min_z) + autobox_add;

	else
	{
		std::cerr << "Unable to read  " << autobox_ligand << "\n";
		exit(-1);
	}
}

//several built-in functions I've designed
void setup_dkoes_terms(custom_terms& t, bool dkoes_score, bool dkoes_score_old,
		bool dkoes_fast)
{
	if (dkoes_score + dkoes_score_old + dkoes_fast != 1)
		throw usage_error("dkoes scoring options are mutually exclusive");

	if (dkoes_score)
	{
		t.add("vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.009900);
		t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.153055);
		t.add("ad4_solvation(d-sigma=3.6,_s/q=0.01097,_c=8)", 0.048934);
		t.add("num_tors_sqr", 0.317267);
		t.add("constant_term", -2.469020);
		/* trained with openbabel partial charges
		 weights.push_back(0.010764); //vdw
		 weights.push_back(-0.156861); //hbond
		 weights.push_back(0.062407); //desolvation
		 weights.push_back(0.337036); //tors sqr
		 weights.push_back(-2.231827); //constant
		 */
	}
	else if (dkoes_score_old)
	{
		t.add("vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.010607);
		t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", 0.197201);
		t.add("num_tors_sqr", .285035);
		t.add("constant_term", -2.585651);

		//desolvation wasn't actually getting counted, but this is the constant
//			weights.push_back(0.044580); //desolvation

	}
	else if (dkoes_fast)
	{
		t.add("vdw(i=4,_j=8,_s=0,_^=100,_c=8)", 0.008962);
		t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", 0.387739);
		t.add("num_tors_sqr", .285035);
		t.add("constant_term", -2.467357);
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

//enum options and their parsers
enum ApproxType
{
	LinearApprox, SplineApprox, Exact, GPU
};

std::istream& operator>>(std::istream& in, ApproxType& type)
{
	using namespace boost::program_options;

	std::string token;
	in >> token;
	if (token == "spline")
		type = SplineApprox;
	else if (token == "linear")
		type = LinearApprox;
	else if (token == "exact")
		type = Exact;
#ifdef SMINA_GPU
	else if (token == "gpu")
	type = GPU;
#endif
	else
		throw validation_error(validation_error::invalid_option_value);
	return in;
}

//helper function for setting molecular data that will wrap data in remark
//if output format is pdb; copy by value because we might change things
static void setMolData(OpenBabel::OBFormat *format, OpenBabel::OBMol& mol,
		std::string attr, std::string value)
{
	using namespace OpenBabel;
	OBPairData* sddata = new OBPairData();
	if (strcmp(format->GetID(), "pdb") == 0 ||
			strcmp(format->GetID(), "ent") == 0 ||
			strcmp(format->GetID(), "pdbqt") == 0)
	{
		//note that pdbqt currently ignores these.. hopefully this will be fixed
		//put attribute name into value so attr can be REMARK
		value = " " + attr + " " + value;
		attr = "REMARK";
	}

	sddata->SetAttribute(attr);
	sddata->SetValue(value);
	mol.SetData(sddata);
}

#ifdef SMINA_GPU

//set the default device to device and exit with error if there are any problems
static void initializeCUDA(int device)
{
	cudaSetDevice(device);
	cudaError_t error;
	cudaDeviceProp deviceProp;
	error = cudaGetDevice(&device);

	if (error != cudaSuccess)
	{
		std::cerr << "cudaGetDevice returned error code " << error << "\n";
		exit(-1);
	}

	error = cudaGetDeviceProperties(&deviceProp, device);

	if (deviceProp.computeMode == cudaComputeModeProhibited)
	{
		std::cerr << "Error: device is running in <Compute Mode Prohibited>, no threads can use ::cudaSetDevice().\n";
		exit(-1);
	}

	if (error != cudaSuccess)
	{
		std::cerr << "cudaGetDeviceProperties returned error code " << error << "\n";
		exit(-1);
	}
}
#endif

//create the initial model from the specified receptor files
//mostly because Matt kept complaining about it, this will automatically create
//pdbqts if necessary using open babel
static void create_init_model(const std::string& rigid_name,
		const std::string& flex_name, model& initm, tee& log)
{
	if (rigid_name.size() > 0)
	{
		std::istream *rigidin_ptr = NULL; //either file or string from openbabel conversion
		ifile rigidin(rigid_name);
		std::stringstream str_rigid;
		rigidin_ptr = &rigidin;
		if (boost::filesystem::extension(rigid_name) != ".pdbqt")
		{
			using namespace OpenBabel;
			obmol_opener fileopener;
			OBConversion conv;
			conv.SetOutFormat("PDBQT");
			conv.AddOption("r", OBConversion::OUTOPTIONS); //rigid molecule, otherwise really slow and useless analysis is triggered
			fileopener.openForInput(conv, rigid_name);
			OBMol rec;
			if (!conv.Read(&rec))
				throw file_error(rigid_name, true);

			rec.AddHydrogens();
			std::string recstr = conv.WriteString(&rec);
			str_rigid.str(recstr);
			rigidin_ptr = &str_rigid;
		}

		if (flex_name.size() > 1
				&& boost::filesystem::extension(flex_name)
						!= ".pdbqt")
			log
					<< "WARNING: flexible receptor does not appear to be in PDBQT format\n";

		if (flex_name.size() > 0) //have flexible residues
		{
			ifile flexin(make_path(flex_name));
			initm = parse_receptor_pdbqt(rigid_name, *rigidin_ptr,
					flex_name, flexin);
		}
		else //just rigid
		{
			initm = parse_receptor_pdbqt(rigid_name, *rigidin_ptr);
		}
	}
}

int main(int argc, char* argv[])
{
	using namespace boost::program_options;
	const std::string version_string =
			"Smina "__DATE__".  Based on AutoDock Vina 1.1.2.";
	const std::string error_message =
			"\n\n\
Please report this error at http://smina.sf.net\n"
					"Please remember to include the following in your problem report:\n\
    * the EXACT error message,\n\
    * your version of the program,\n\
    * the type of computer system you are running it on,\n\
	* all command line options,\n\
	* configuration file (if used),\n\
    * ligand file as PDBQT,\n\
    * receptor file as PDBQT,\n\
	* flexible side chains file as PDBQT (if used),\n\
	* output file as PDBQT (if any),\n\
    * input (if possible),\n\
	* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";

	const std::string cite_message =
			"   _______  _______ _________ _        _______ \n"
					"  (  ____ \\(       )\\__   __/( (    /|(  ___  )\n"
					"  | (    \\/| () () |   ) (   |  \\  ( || (   ) |\n"
					"  | (_____ | || || |   | |   |   \\ | || (___) |\n"
					"  (_____  )| |(_)| |   | |   | (\\ \\) ||  ___  |\n"
					"        ) || |   | |   | |   | | \\   || (   ) |\n"
					"  /\\____) || )   ( |___) (___| )  \\  || )   ( |\n"
					"  \\_______)|/     \\|\\_______/|/    )_)|/     \\|\n"
					"\n\nsmina is based off AutoDock Vina. Please cite appropriately.\n";

	try
	{
		std::string rigid_name, flex_name, config_name, log_name, atom_name;
		std::vector<std::string> ligand_names;
		std::string out_name;
		std::string ligand_names_file;
		std::string custom_file_name;
		std::string usergrid_file_name;
		fl center_x = 0, center_y = 0, center_z = 0, size_x = 0, size_y = 0,
				size_z = 0;
		fl autobox_add = 8;
		fl out_min_rmsd = 1;
		std::string autobox_ligand;
		int cpu = 0, seed, exhaustiveness, verbosity, num_modes = 9,
				device = 0;
		fl energy_range = 2.0;

		// -0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 0.05846
		fl weight_gauss1 = -0.035579;
		fl weight_gauss2 = -0.005156;
		fl weight_repulsion = 0.840245;
		fl weight_hydrophobic = -0.035069;
		fl weight_hydrogen = -0.587439;
		fl weight_rot = 0.05846;
		fl user_grid_lambda;
		bool score_only = false, local_only = false, randomize_only = false,
				help = false, help_hidden = false, version = false;
		bool dominimize = false;
		bool dkoes_score = false;
		bool ad4_score = false;
		bool dkoes_score_old = false;
		bool dkoes_fast = false;
		bool quiet = false;
		bool accurate_line = false;
		bool flex_hydrogens = false;
		bool include_atom_terms = false;
		bool gpu_on = false;
		bool print_terms = false;
		bool add_hydrogens = true;
		bool no_lig = false;
		minimization_params minparms;
		ApproxType approx = LinearApprox;
		fl approx_factor = 32;

		positional_options_description positional; // remains empty

		options_description inputs("Input");
		inputs.add_options()
		("receptor,r", value<std::string>(&rigid_name),
				"rigid part of the receptor (PDBQT)")
		("flex", value<std::string>(&flex_name),
				"flexible side chains, if any (PDBQT)")
		("ligand,l", value<std::vector<std::string> >(&ligand_names),
				"ligand(s)");
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
				"Amount of buffer space to add to auto-generated box (default 8)")
		("no_lig", bool_switch(&no_lig),
				"no ligand; for sampling/minimizing flexible residues");

		//options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
		options_description outputs("Output (optional)");
		outputs.add_options()
		("out,o", value<std::string>(&out_name),
				"output file name, format taken from file extension")
		("log", value<std::string>(&log_name), "optionally, write log file")
		("atom_terms", value<std::string>(&atom_name),
				"optionally write per-atom interaction term values")
		("atom_term_data", bool_switch(&include_atom_terms),
				"embedded per-atom interaction terms in output sd data");

		options_description scoremin("Scoring and minimization options");
		scoremin.add_options()
		("custom_scoring", value<std::string>(&custom_file_name),
				"custom scoring function file")
		("score_only", bool_switch(&score_only), "score provided ligand pose")
		("local_only", bool_switch(&local_only),
				"local search only using autobox (you probably want to use --minimize)")
		("minimize", bool_switch(&dominimize),
				"energy minimization")
		("randomize_only", bool_switch(&randomize_only),
				"generate random poses, attempting to avoid clashes")
		("minimize_iters",
				value<unsigned>(&minparms.maxiters)->default_value(0),
				"number iterations of steepest descent; default scales with rotors and usually isn't sufficient for convergence")
		("accurate_line", bool_switch(&accurate_line),
				"use accurate line search")
		("minimize_early_term", bool_switch(&minparms.early_term),
				"Stop minimization before convergence conditions are fully met.")
		("approximation", value<ApproxType>(&approx),
				"approximation (linear, spline, or exact) to use")
		("factor", value<fl>(&approx_factor),
				"approximation factor: higher results in a finer-grained approximation")
		("print_terms", bool_switch(&print_terms),
				"Print all available terms with default parameterizations");

		options_description hidden("Hidden options for internal testing");
		hidden.add_options()
		("dkoes_scoring", bool_switch(&dkoes_score),
				"Use my custom scoring function")
		("dkoes_scoring_old", bool_switch(&dkoes_score_old),
				"Use old (vdw+hbond) scoring function")
		("dkoes_fast", bool_switch(&dkoes_fast), "VDW+nrot only")
		("ad4_scoring", bool_switch(&ad4_score),
				"Approximation of Autodock 4 scoring")
		("user_grid", value<std::string>(&usergrid_file_name),
				"Autodock map file for user grid data based calculations, not implemented yet")
		("verbosity", value<int>(&verbosity)->default_value(1),
				"Adjust the verbosity of the output, default: 1")
		("user_grid_lambda", value<fl>(&user_grid_lambda)->default_value(-1.0),
						"Scales user_grid and functional scoring");

		options_description misc("Misc (optional)");
		misc.add_options()
		("cpu", value<int>(&cpu),
				"the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
		("seed", value<int>(&seed), "explicit random seed")
		("exhaustiveness", value<int>(&exhaustiveness)->default_value(8),
				"exhaustiveness of the global search (roughly proportional to time)")
		("num_modes", value<int>(&num_modes)->default_value(9),
				"maximum number of binding modes to generate")
		("energy_range", value<fl>(&energy_range)->default_value(3.0),
				"maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
		("min_rmsd_filter", value<fl>(&out_min_rmsd)->default_value(1.0),
				"rmsd value used to filter final poses to remove redundancy")
		("quiet,q", bool_switch(&quiet), "Suppress output messages")
		("addH", value<bool>(&add_hydrogens),
				"automatically add hydrogens in ligands (on by default)")
		("flex_hydrogens", bool_switch(&flex_hydrogens),
				"Enable torsions effecting only hydrogens (e.g. OH groups). This is stupid but provides compatibility with Vina.")
				#ifdef SMINA_GPU
				("device", value<int>(&device)->default_value(0), "GPU device to use")
				("gpu", bool_switch(&gpu_on), "Turn on GPU acceleration")
#endif
				;

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
		desc.add(inputs).add(search_area).add(outputs).add(scoremin).
				add(hidden).add(misc).add(config).add(info);
		desc_simple.add(inputs).add(search_area).add(scoremin).
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
				path name = make_path(config_name);
				ifile config_stream(name);
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
		if (print_terms)
		{
			custom_terms t;
			t.print_available_terms(std::cout);
			return 0;
		}
#ifdef SMINA_GPU
		initializeCUDA(device);
#endif

		set_fixed_rotable_hydrogens(!flex_hydrogens);

		if (dominimize) //set default settings for minimization
		{
			if (minparms.maxiters == 0)
				minparms.maxiters = 10000; //will presumably converge
			local_only = true;
			minparms.type = minimization_params::BFGSAccurateLineSearch;

			if (!vm.count("approximation"))
				approx = SplineApprox;
			if (!vm.count("factor"))
				approx_factor = 10;
		}

		if (accurate_line)
		{
			minparms.type = minimization_params::BFGSAccurateLineSearch;
		}

		bool search_box_needed = !(score_only || local_only); // randomize_only and local_only still need the search space; dkoes - for local get box from ligand
		bool output_produced = !score_only;
		bool receptor_needed = !randomize_only;

		if (receptor_needed)
		{
			if (vm.count("receptor") <= 0)
			{
				std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n"
						<< desc_simple << '\n';
				return 1;
			}
		}

		if (ligand_names.size() == 0)
		{
			if (!no_lig)
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
		else if (no_lig) //ligand specified with no_lig
		{
			std::cerr << "Ligand specified with --no_lig.\n"
					<< "\nCorrect usage:\n"
					<< desc_simple << '\n';
			return 1;
		}
		if (cpu < 1)
			cpu = 1;
		if (vm.count("seed") == 0)
			seed = auto_seed();
		if (exhaustiveness < 1)
			throw usage_error("exhaustiveness must be 1 or greater");
		if (num_modes < 1)
			throw usage_error("num_modes must be 1 or greater");
		sz max_modes_sz = static_cast<sz>(num_modes);

		boost::optional<std::string> flex_name_opt;
		if (vm.count("flex"))
			flex_name_opt = flex_name;

		if (vm.count("flex") && !vm.count("receptor"))
			throw usage_error(
					"Flexible side chains are not allowed without the rest of the receptor"); // that's the only way parsing works, actually

		tee log(quiet);
		if (vm.count("log") > 0)
			log.init(log_name);

		std::ofstream atomoutfile;
		if (vm.count("atom_terms") > 0)
			atomoutfile.open(atom_name.c_str());

		if (autobox_ligand.length() > 0)
		{
			setup_autobox(autobox_ligand, autobox_add,
					center_x, center_y, center_z,
					size_x, size_y, size_z);
		}
		if (search_box_needed && autobox_ligand.length() == 0)
		{
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

		log << cite_message << '\n';

		grid_dims gd; // n's = 0 via default c'tor
		grid_dims user_gd;
		grid user_grid;

		flv weights;

		//dkoes, set the scoring function
		custom_terms t;
		if(user_grid_lambda != -1.0){
			t.set_scaling_factor(user_grid_lambda);
		}
		if (custom_file_name.size() > 0)
		{
			ifile custom_file(make_path(custom_file_name));
			t.add_terms_from_file(custom_file);
		}
		else if (dkoes_score || dkoes_score_old || dkoes_fast)
		{
			//my own built-in scoring functions
			setup_dkoes_terms(t, dkoes_score, dkoes_score_old, dkoes_fast);
		}
		else if (ad4_score)
		{
			t.add("vdw(i=6,_j=12,_s=0,_^=100,_c=8)", 0.1560);
			t.add("non_dir_h_bond_lj(o=-0.7,_^=100,_c=8)", -0.0974);
			t.add("ad4_solvation(d-sigma=3.5,_s/q=0.01097,_c=8)", 0.1159);
			t.add("electrostatic(i=1,_^=100,_c=8)", 0.1465);
			t.add("num_tors_add", 0.2744);
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
				
		log << std::setw(12) << std::left << "Weights" << " Terms\n" << t
				<< "\n";

		if (usergrid_file_name.size() > 0)
		{
			ifile user_in(make_path(usergrid_file_name));
			fl ug_scaling_factor = 1.0;
			if(user_grid_lambda != -1.0){
				ug_scaling_factor = 1 - user_grid_lambda;
			}
			setup_user_gd(user_gd, user_in);
			user_grid.init(user_gd, user_in, ug_scaling_factor); //initialize user grid
		}

		const fl granularity = 0.375;
		if (search_box_needed)
		{
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

		if (vm.count("cpu") == 0)
		{
			unsigned num_cpus = boost::thread::hardware_concurrency();
			if (verbosity > 1)
			{
				if (num_cpus > 0)
					log << "Detected " << num_cpus << " CPU"
							<< ((num_cpus > 1) ? "s" : "") << '\n';
				else
					log << "Could not detect the number of CPUs, using 1\n";
			}
			if (num_cpus > 0)
				cpu = num_cpus;
			else
				cpu = 1;
		}
		if (cpu < 1)
			cpu = 1;
		if (verbosity > 1 && exhaustiveness < cpu)
			log
					<< "WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs\n";

		//dkoes - parse in receptor once
		model initm;

		create_init_model(rigid_name, flex_name, initm, log);

		//dkoes, hoist precalculation outside of loop
		weighted_terms wt(&t, t.weights());

		boost::shared_ptr<precalculate> prec;

		if (gpu_on || approx == GPU)
		{ //don't get a choice
#ifdef SMINA_GPU
		prec = boost::shared_ptr<precalculate>(new precalculate_gpu(wt, approx_factor));
#endif
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
		obmol_opener outfileopener;
		OBConversion outconv;
		if (out_name.length() > 0)
		{
			outfileopener.openForOutput(outconv, out_name);
			VINA_CHECK(outconv.SetInFormat("PDBQT"));
		}

		if(score_only) //output header
		{
			std::vector<std::string> enabled_names = t.get_names(true);
			log << "## Name";
			VINA_FOR_IN(i, enabled_names)
			{
				log << " " << enabled_names[i];
			}
			for (unsigned i = 0, n = t.conf_independent_terms.size(); i < n; i++)
			{
				log << " " << t.conf_independent_terms[i].name;
			}
			log << "\n";
		}
		//loop over input ligands
		for (unsigned l = 0, nl = ligand_names.size(); l < nl; l++)
		{
			doing(verbosity, "Reading input", log);
			const std::string& ligand_name = ligand_names[l];
			boost::filesystem::path lpath(ligand_name);

			//parse with open babel
			OBConversion conv;
			obmol_opener infileopener;

			if (ligand_name.size() > 0) //is zero if no_lig
			{
				infileopener.openForInput(conv, ligand_name);
				VINA_CHECK(conv.SetOutFormat("PDBQT"));
			}

			//process input molecules one at a time
			OBMol mol;
			unsigned i = 0;
			while (no_lig || conv.Read(&mol))
			{
				model m = initm;
				std::string name = mol.GetTitle();
				m.set_name(name);
				try
				{
					//this is suboptimal: do not read/write pdbqt with openbabel
					//because it will lose information about rigid bonds
					if (no_lig)
					{
						no_lig = false; //only enter loop once
					}
					else if (lpath.extension() == ".pdbqt")
					{
						m.append(parse_ligand_pdbqt(lpath));
					}
					else
					{
						if (add_hydrogens)
							mol.AddHydrogens(); //needed for atom typing
						std::string pdbqt = conv.WriteString(&mol);
						std::stringstream pdbqtStream(pdbqt);

						m.append(parse_ligand_stream_pdbqt(ligand_name,
								pdbqtStream));
					}
				} catch (parse_error& e)
				{
					std::cerr << "\n\nParse error with molecule "
							<< mol.GetTitle() << " in file \""
							<< e.file.string() << "\": " << e.reason
							<< '\n';
					continue;
				}

				if (local_only)
				{
					//dkoes - for convenience get box from model
					gd = m.movable_atoms_box(autobox_add, granularity);
				}

				boost::optional<model> ref;
				done(verbosity, log);

				std::stringstream output;
				std::vector<resultInfo> results;

				main_procedure(m, *prec, ref, score_only,
						local_only, randomize_only,
						false, // no_cache == false
						atomoutfile.is_open() || include_atom_terms, gpu_on,
						gd, exhaustiveness, minparms, wt, cpu, seed, verbosity,
						max_modes_sz, energy_range, out_min_rmsd, log, results,
						user_grid);

				if (outconv.GetOutStream() != NULL)
				{
					//write out molecular data
					for (unsigned j = 0, m = results.size(); j < m; j++)
					{
						if (results[j].mol.length() > 0)
							outconv.ReadString(&mol, results[j].mol); //otherwise keep orig mol
						mol.DeleteData(OBGenericDataType::PairData); //remove remarks

						setMolData(outconv.GetOutFormat(), mol,
								"minimizedAffinity",
								boost::lexical_cast<std::string>(
										results[j].energy));

						if (results[j].rmsd >= 0)
						{
							setMolData(outconv.GetOutFormat(), mol,
									"minimizedRMSD",
									boost::lexical_cast<std::string>(
											results[j].rmsd));
						}
						if (include_atom_terms)
						{
							std::stringstream astr;
							results[j].writeAtomValues(astr, &wt);
							setMolData(outconv.GetOutFormat(), mol,
									"atomic_interaction_terms", astr.str());
						}
						mol.SetTitle(name); //otherwise lose space separated names
						outconv.SetOutputIndex(j + 2); //workaround openbabel bug #859
						outconv.Write(&mol);
					}
				}
				if (atomoutfile)
				{
					for (unsigned j = 0, m = results.size(); j < m; j++)
					{
						results[j].writeAtomValues(atomoutfile, &wt);
					}
				}
				mol.Clear();
				i++;
			}
		}
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
		std::cerr << "\n\nUsage error: " << e.what() << ".\n";
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
