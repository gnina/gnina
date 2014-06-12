/*
 * MinimizationQuery.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#include "MinimizationQuery.h"
#include <sstream>

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
}

MinimizationParameters::~MinimizationParameters()
{
	if(wt) delete wt;
	if(prec) delete prec;
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
//allocates and returns a result structure
MinimizationQuery::Result* MinimizationQuery::minimize(model& m)
{

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
				m.append(tmp.m);

				Result *result = q->minimize(m);
				if(result != NULL)
					results.push_back(result);
			}

			//add computed results
			lock_guard<shared_mutex> lock(q->results_mutex);
			for(unsigned i = 0, n = results.size(); i < n; i++)
			{
				q->allResults.push_back(results[i]);
			}
		}

	}
}

