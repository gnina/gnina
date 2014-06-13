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

using namespace boost;

MinimizationParameters::MinimizationParameters() :
		wt(NULL), nthreads(1)
{
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

MinimizationParameters::~MinimizationParameters()
{
	if (wt)
		delete wt;
	if (prec)
		delete prec;
	if (exact_prec)
		delete exact_prec;
	if (nnc)
		delete nnc;
}

MinimizationQuery::~MinimizationQuery()
{
	checkThread();
	assert(minimizationSpawner == NULL); //should not be deleted while minimization is running
	for (unsigned i = 0, n = allResults.size(); i < n; i++)
	{
		delete allResults[i];
	}
	allResults.clear();
}

void MinimizationQuery::checkThread()
{
	//clean up minimizer threads
	if (minimizationSpawner != NULL)
	{
		if (minimizationSpawner->timed_join(posix_time::time_duration(0, 0, 0)))
		{
			delete minimizationSpawner;
			minimizationSpawner = NULL; // thread is done
		}
	}
}

//return true if down minimizing
bool MinimizationQuery::finished()
{
	checkThread();
	return minimizationSpawner == NULL;
}

//execute the query - set up the minimization threads and communication queue
void MinimizationQuery::execute(bool block/*=false*/)
{
	assert(minimizationSpawner == NULL);
	minimizationSpawner = new thread(thread_startMinimization, this);

	if (block)
		minimizationSpawner->join();
}

//launch and wait for the minimization threads
void MinimizationQuery::thread_startMinimization(MinimizationQuery *query)
{
	//generate a thread for each set of databases
	thread_group minthreads;

	for (unsigned t = 0; t < query->minparm.nthreads; t++)
	{
		minthreads.add_thread(new thread(thread_minimize, query));
	}
	minthreads.join_all(); //block till all done
	query->io->close();
}

//thread safe minimization of m
//allocates and returns a result structure, caller takes responsibility for memory
MinimizationQuery::Result* MinimizationQuery::minimize(model& m)
{
	static const vec authentic_v(10, 10, 10);
	static const grid empty_grid;
	static const fl autobox_add = 8;
	static const fl granularity = 0.375;

	vecv origcoords = m.get_heavy_atom_movable_coords();
	fl e = max_fl;
	conf c = m.get_initial_conf();
	output_type out(c, e);

	//do minimization
	grid_dims gd = m.movable_atoms_box(autobox_add, granularity);
	change g(m.get_size());
	quasi_newton quasi_newton_par(minparm.minparms);
	szv_grid_cache gridcache(m, minparm.prec->cutoff_sqr());

	non_cache nc(gridcache, gd, minparm.prec);
	//it rarely takes more than one try
	fl slope = 10;
	unsigned i;
	for (i = 0; i < 3; i++)
	{
		nc.setSlope(slope);
		quasi_newton_par(m, *minparm.prec, nc, out, g, authentic_v, empty_grid);
		m.set(out.c); // just to be sure
		if (nc.within(m))
			break;
		slope *= 10;
	}

	if (i == 3) //couldn't stay in box
		out.e = max_fl;

	fl intramolecular_energy = m.eval_intramolecular(*minparm.exact_prec,
			authentic_v, out.c);
	e = m.eval_adjusted(*minparm.wt, *minparm.exact_prec, *minparm.nnc,
			authentic_v, out.c,
			intramolecular_energy, empty_grid);

	vecv newcoords = m.get_heavy_atom_movable_coords();
	assert(newcoords.size() == origcoords.size());

	fl rmsd = 0;
	for (unsigned i = 0, n = newcoords.size(); i < n; i++)
	{
		rmsd += (newcoords[i] - origcoords[i]).norm_sqr();
	}
	rmsd /= newcoords.size();
	rmsd = sqrt(rmsd);

	//construct result
	Result *result = new Result();
	result->score = e;
	result->rmsd = rmsd;
	result->name = m.get_name();
	stringstream str;
	m.write_sdf(str);
	result->sdf = str.str();

	return result;
}

//read a chunk of ligands at a time
//return false if there's definitely nothing else to read
bool MinimizationQuery::thread_safe_read(vector<LigandData>& ligands)
{
	ligands.clear();
	unique_lock<mutex> lock(io_mutex);

	if (!*io)
		return false;

	for (unsigned i = 0; i < chunk_size; i++)
	{
		try
		{
			LigandData data;

			boost::archive::binary_iarchive serialin(*io,
					boost::archive::no_header | boost::archive::no_tracking);
			serialin >> data.numtors;
			serialin >> data.p;
			serialin >> data.c;
		}
		catch (boost::archive::archive_exception& e)
		{
			return ligands.size() > 0; //need to minimize last set of ligands
		}
	}
}

//read chunks of ligands, minimize them, and store the result
void MinimizationQuery::thread_minimize(MinimizationQuery* q)
{
	vector<LigandData> ligands;
	while (q->thread_safe_read(ligands))
	{
		if (q->stopQuery) //cancelled
			return;
		else if (ligands.size() == 0) //nothing to do yet
			usleep(10);
		else //let's minimize!
		{
			vector<Result*> results;
			for (unsigned i = 0, n = ligands.size(); i < n; i++)
			{
				//construct model
				LigandData& l = ligands[i];
				model m = q->initm;
				non_rigid_parsed nr;
				postprocess_ligand(nr, l.p, l.c, l.numtors);

				pdbqt_initializer tmp;
				tmp.initialize_from_nrp(nr, l.c, true);
				tmp.initialize(nr.mobility_matrix());
				m.set_name(l.c.sdftext.name);

				if (q->hasReorient)
					l.reorient.reorient(tmp.m.coordinates());

				m.append(tmp.m);

				Result *result = q->minimize(m);
				if (result != NULL)
					results.push_back(result);
			}

			//add computed results
			lock_guard<shared_mutex> lock(q->results_mutex);
			for (unsigned i = 0, n = results.size(); i < n; i++)
			{
				q->allResults.push_back(results[i]);
			}
		}

	}
}

