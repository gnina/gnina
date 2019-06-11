/*
 * MinimizationQuery.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#include <sstream>

#include "MinimizationQuery.h"
#include "conf.h"
#include "non_cache.h"
#include "quasi_newton.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/unordered_set.hpp>
#include <boost/timer/timer.hpp>

using namespace boost;

MinimizationParameters::MinimizationParameters()
    : wt(NULL), nthreads(1) {
  //default settings
  minparms.maxiters = 10000;
  minparms.type = minimization_params::BFGSAccurateLineSearch;

  //vina scoring function
  t.add("gauss(o=0,_w=0.5,_c=8)", -0.035579);
  t.add("gauss(o=3,_w=2,_c=8)", -0.005156);
  t.add("repulsion(o=0,_c=8)", 0.840245);
  t.add("hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069);
  t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439);
  t.add("num_tors_div", 5 * 0.05846 / 0.1 - 1);

  wt = new weighted_terms(&t, t.weights());
  prec = new precalculate_splines(*wt, 10.0);
  exact_prec = new precalculate_exact(*wt);
  nnc = new naive_non_cache(exact_prec);
}

MinimizationParameters::~MinimizationParameters() {
  if (wt) delete wt;
  if (prec) delete prec;
  if (exact_prec) delete exact_prec;
  if (nnc) delete nnc;
}

MinimizationQuery::~MinimizationQuery() {
  checkThread();
  assert(minimizationSpawner == NULL); //should not be deleted while minimization is running
  for (unsigned i = 0, n = allResults.size(); i < n; i++) {
    delete allResults[i];
  }
  allResults.clear();
}

void MinimizationQuery::checkThread() {
  //clean up minimizer threads
  if (minimizationSpawner != NULL) {
    if (minimizationSpawner->timed_join(posix_time::time_duration(0, 0, 0))) {
      delete minimizationSpawner;
      minimizationSpawner = NULL; // thread is done
    }
  }
}

//return true if down minimizing
bool MinimizationQuery::finished() {
  checkThread();
  return isFinished; //do not want to return true before minimization even starts
}

//execute the query - set up the minimization threads and communication queue
void MinimizationQuery::execute(bool block/*=false*/) {
  assert(minimizationSpawner == NULL);
  minimizationSpawner = new boost::thread(thread_startMinimization, this);

  if (block) minimizationSpawner->join();
}

//launch and wait for the minimization threads
void MinimizationQuery::thread_startMinimization(MinimizationQuery *query) {
  //generate a thread for each set of databases
  thread_group minthreads;

  timer::cpu_timer mintime;

  for (unsigned t = 0; t < query->minparm.nthreads; t++) {
    minthreads.add_thread(new boost::thread(thread_minimize, query));
  }
  minthreads.join_all(); //block till all done

  query->minTime = mintime.elapsed().wall / 1e9;
  query->io->close();
  query->isFinished = true;
}

//thread safe minimization of m
//allocates and returns a result structure, caller takes responsibility for memory
MinimizationQuery::Result* MinimizationQuery::minimize(model& m) {
  static const grid empty_grid;
  static const fl autobox_add = 4;
  static const fl granularity = 0.375;
  vec authentic_v(10, 10, 10); //"soft" minimization

  if (isFrag) authentic_v = vec(100, 100, 100); //harder since we can't resolve clashes with fixed parts

  vecv origcoords = m.get_heavy_atom_movable_coords();
  fl e = max_fl;

  //do minimization
  grid_dims gd = m.movable_atoms_box(autobox_add, granularity);
  szv_grid_cache gridcache(m, minparm.prec->cutoff_sqr());
  non_cache nc(gridcache, gd, minparm.prec);
  conf c = m.get_initial_conf(nc.move_receptor());
  output_type out(c, e);
  change g(m.get_size(), nc.move_receptor());

  //regular minimization
  quasi_newton quasi_newton_par(minparm.minparms);

  //it rarely takes more than one try
  fl slope = 10;
  unsigned i;
  for (i = 0; i < 3; i++) {
    nc.setSlope(slope);
    quasi_newton_par(m, *minparm.prec, nc, out, g, authentic_v, empty_grid);
    m.set(out.c); // just to be sure
    if (nc.within(m)) break;
    slope *= 10;
  }

  if (i == 3) //couldn't stay in box
    out.e = max_fl;

  if (isFrag) {
    e = m.eval_flex(*minparm.exact_prec, authentic_v, out.c, numProteinAtoms);
  } else { //standard ligand stuff
    fl intramolecular_energy = m.eval_intramolecular(*minparm.exact_prec,
        authentic_v, out.c);
    e = m.eval_adjusted(*minparm.wt, *minparm.exact_prec, *minparm.nnc,
        authentic_v, out.c, intramolecular_energy, empty_grid);
  }

  vecv newcoords = m.get_heavy_atom_movable_coords();
  assert(newcoords.size() == origcoords.size());

  fl rmsd = 0;
  for (unsigned i = 0, n = newcoords.size(); i < n; i++) {
    rmsd += (newcoords[i] - origcoords[i]).norm_sqr();
  }

  if (newcoords.size() > 0) { //not totally rigid
    rmsd /= newcoords.size();
    rmsd = sqrt(rmsd);
  }

  //construct result
  Result *result = new Result();
  result->score = e;
  result->rmsd = rmsd;
  result->name = m.get_name();
  stringstream str;

  if (m.num_ligands() > 0) {
    m.write_sdf(str);
  } else { //if no ligs, assume optimizing "residue"
    m.write_flex_sdf(str);
  }
  result->sdf = str.str();

  return result;
}

//read a chunk of ligands at a time
//return false if there's definitely nothing else to read
bool MinimizationQuery::thread_safe_read(vector<LigandData>& ligands) {
  ligands.clear();
  boost::unique_lock<boost::mutex> lock(io_mutex);

  if (!io_strm) return false;

  boost::archive::binary_iarchive serialin(io_strm,
      boost::archive::no_header | boost::archive::no_tracking);

  for (unsigned i = 0; i < chunk_size; i++) {
    try {
      LigandData data;

      if (hasReorient) data.reorient.read(io_strm);
      serialin >> data.numtors;
      serialin >> data.p;
      serialin >> data.c;
      data.origpos = io_position;
      io_position++;

      ligands.push_back(data);
    } catch (boost::archive::archive_exception& e) {
      return ligands.size() > 0; //need to minimize last set of ligands
    } catch (...) {
      return false; //give up
    }
  }
  return !!io_strm;
}

//read chunks of ligands, minimize them, and store the result
void MinimizationQuery::thread_minimize(MinimizationQuery* q) {
  unsigned cnt = 0;
  try {
    vector<LigandData> ligands;
    while (q->thread_safe_read(ligands)) {
      if (q->stopQuery) //cancelled
        return;
      else
        if (ligands.size() == 0) //nothing to do yet
          usleep(10);
        else //let's minimize!
        {
          vector<Result*> results;
          for (unsigned i = 0, n = ligands.size(); i < n; i++) {
            //construct model
            cnt++;
            LigandData& l = ligands[i];
            model m = q->initm;

            if (q->hasReorient) l.reorient.reorient(l.p);

            non_rigid_parsed nr;
            pdbqt_initializer tmp;

            if (q->isFrag) {
              //treat as residue
              postprocess_residue(nr, l.p, l.c);
            } else {
              postprocess_ligand(nr, l.p, l.c, l.numtors);
            }

            tmp.initialize_from_nrp(nr, l.c, !q->isFrag);
            tmp.initialize(nr.mobility_matrix());
            m.set_name(l.c.sdftext.name);

            m.append(tmp.m);

            Result *result = q->minimize(m);
            result->orig_position = l.origpos;
            if (result != NULL) results.push_back(result);
          }

          //add computed results
          boost::lock_guard<shared_mutex> lock(q->results_mutex);
          for (unsigned i = 0, n = results.size(); i < n; i++) {
            results[i]->position = q->allResults.size();
            q->allResults.push_back(results[i]);
          }
        }

    }
  } catch (...) //don't die
  {
    q->cancel();
  }
}

//output the mol at position pos
void MinimizationQuery::outputMol(unsigned pos, ostream& out) {
  Result *res = NULL;
  //thread safe grab from allResults
  results_mutex.lock_shared();
  if (pos < allResults.size()) res = allResults[pos];
  results_mutex.unlock_shared();

  if (res != NULL) //empty result on error
  {
    out << res->sdf;
    out << "$$$$\n";
  }
}

//comparison object for results
class MinimizationQuery::ResultsSorter {
    const MinimizationFilters& filter;

  public:
    ResultsSorter(const MinimizationFilters& f)
        : filter(f) {
    }

    bool operator()(const Result* lhs, const Result* rhs) const {
      bool res = false;
      switch (filter.sort) {
      case MinimizationFilters::Score:
        res = lhs->score < rhs->score;
        break;
      case MinimizationFilters::RMSD:
        res = lhs->rmsd < rhs->rmsd;
        break;
      case MinimizationFilters::OrigPos:
        res = lhs->orig_position < rhs->orig_position;
        break;
      }

      if (filter.reverseSort) res = !res;
      return res;
    }

};
//copies allResults (safely) into results which is then sorted and filter according
//to filter; returns the total number of results _before_filtering
//does not truncate results by start/num of filter
unsigned MinimizationQuery::loadResults(const MinimizationFilters& filter,
    vector<Result*>& results) {
  results_mutex.lock_shared();
  results = allResults;
  results_mutex.unlock_shared();

  //filter
  unsigned total = results.size();
  unsigned i = 0;
  while (i < results.size()) {
    Result *res = results[i];
    if (res->rmsd > filter.maxRMSD || res->score > filter.maxScore) {
      //needs to be removed
      swap(results[i], results.back());
      results.pop_back();
      //do NOT increment i, it's a new one
    } else {
      i++; //go to next
    }
  }

  //now sort
  ResultsSorter sorter(filter);
  sort(results.begin(), results.end(), sorter);

  //uniquification, take best after sort
  if (filter.unique) {
    unordered_set<string> seen;
    vector<Result*> tmpres;
    tmpres.reserve(results.size());

    for (unsigned i = 0, n = results.size(); i < n; i++) {
      Result *res = results[i];
      if (seen.count(res->name) == 0) //first time we've seen it
          {
        tmpres.push_back(res);
        seen.insert(res->name);
      }
    }
    swap(results, tmpres);

  }
  return total;
}

//output text formated data
void MinimizationQuery::outputData(const MinimizationFilters& f, ostream& out) {
  checkThread();
  vector<Result*> results;
  unsigned total = loadResults(f, results);

  //first line is status header with doneness and number done and filtered number
  out << finished() << " " << total << " " << results.size() << " " << minTime
      << "\n";

  unsigned end = f.start + f.num;
  if (end > results.size() || f.num == 0) end = results.size();
  for (unsigned i = f.start; i < end; i++) {
    Result *res = results[i];
    out << res->position << "," << res->orig_position << "," << res->name << ","
        << res->score << "," << res->rmsd << "\n";
  }
}

//output json formated data, based off of datatables, does not include opening/closing brackets
void MinimizationQuery::outputJSONData(const MinimizationFilters& f, int draw,
    ostream& out) {
  checkThread();
  vector<Result*> results;
  unsigned total = loadResults(f, results);

  //first line is status header with doneness and number done and filtered number
  out << "{\n";
  out << "\"finished\": " << finished() << ",\n";
  out << "\"recordsTotal\": " << total << ",\n";
  out << "\"recordsFiltered\": " << results.size() << ",\n";
  out << "\"time\": " << minTime << ",\n";
  out << "\"draw\": " << draw << ",\n";
  out << "\"data\": [\n";

  unsigned end = f.start + f.num;
  if (end > results.size() || f.num == 0) end = results.size();
  for (unsigned i = f.start; i < end; i++) {
    Result *res = results[i];
    out << "[" << res->position << "," << res->orig_position << ",\""
        << res->name << "\"," << res->score << "," << res->rmsd << "]";
    if (i != end - 1) out << ",";
    out << "\n";
  }
  out << "]}\n";
}

//write out all results in sdf.gz format
void MinimizationQuery::outputMols(const MinimizationFilters& f, ostream& out) {
  vector<Result*> results;
  loadResults(f, results);

  //gzip output
  boost::iostreams::filtering_stream<boost::iostreams::output> strm;
  strm.push(boost::iostreams::gzip_compressor());
  strm.push(out);
  for (unsigned i = 0, n = results.size(); i < n; i++) {
    Result *r = results[i];
    strm << r->sdf;
    //include sddata for affinity/rmsd
    strm << ">  <minimizedAffinity>\n";
    strm << std::fixed << std::setprecision(5) << r->score << "\n\n";

    if (r->rmsd >= 0) {
      strm << ">  <minimizedRMSD>\n";
      strm << std::fixed << std::setprecision(5) << r->rmsd << "\n\n";
    }
    strm << "$$$$\n";
  }
}
