/*
 * QueryManager.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#include "QueryManager.h"
#include "Reorienter.h"
#include "MinimizationQuery.h"
#include <boost/algorithm/string.hpp>

using namespace boost;

//add a query, return zero if unsuccessful
unsigned QueryManager::add(unsigned oldqid, stream_ptr io)
{
	mu.lock();
	//explicitly remove old query
	if (oldqid > 0 && queries.count(oldqid) > 0)
	{
		QueryPtr oldq = queries[oldqid];
		oldq->cancel();
		if (oldq->finished())
			queries.erase(oldqid); //removates shared ptr reference from map if not still minimizing
		//otherwise need to wait to purge
	}
	mu.unlock();

	//read receptor info and rotation/translation info, but leave ligand for minimizer to stream
	string str;
	*io >> str;
	if (str != "receptor")
	{
		cerr << "No receptor\n";
		*io << "ERROR\nNo receptor\n";
		return 0;
	}
	unsigned rsize = 0; //size of receptor string, must be in pdbqt
	*io >> rsize;
	if (rsize == 0)
	{
		cerr << "invalid receptor size\n";
		*io << "ERROR\nInvalid receptor size\n";
		return 0;
	}
	getline(*io, str); //absorb newline
	string recstr(rsize, '\0'); //note that c++ strings are built with null at the end
	io->read(&recstr[0], rsize);
	//does the ligand data have to be reoriented?
	bool hasR = false;
	*io >> hasR;
	//absorb newline before ligand data
	getline(*io, str);

	//attempt to create query
	QueryPtr q;
	try
	{
		q = QueryPtr(new MinimizationQuery(minparm, recstr, io, hasR));
	} catch (parse_error& pe) //couldn't read receptor
	{
		cerr << "couldn't read receptor\n";
		*io << "ERROR\n" << pe.reason << "\n";
		return 0;
	}  catch (...) { //output below
	}

	if (q)
	{
		//protect access to nextID and queries
		mu.lock();
		unsigned id = nextID++;
		queries[id] = q;
		mu.unlock();

		q->execute(); //don't wait for result

		return id;
	}
	else //error,
	{
		*io << "ERROR\nCould not construct minimization query\n";
		return 0;
	}
}

//count types of queries
void QueryManager::getCounts(unsigned& active, unsigned& inactive,
		unsigned& defunct)
{
	active = inactive = defunct = 0;
	boost::lock_guard<boost::mutex> L(mu);

	for (QueryMap::iterator itr = queries.begin(), end = queries.end(); itr
			!= end; itr++)
	{
		QueryPtr q = itr->second;
		if (!q->finished())
			active++;
		else if (q->idle() > timeout)
			defunct++;
		else
			inactive++;
	}
}

QueryPtr QueryManager::get(unsigned qid)
{
	boost::unique_lock<boost::mutex>(mu);

	if (queries.count(qid) == 0)
		return QueryPtr();
	QueryPtr q = queries[qid];
	if (q == NULL)
		return QueryPtr();

	q->access();
	return q;
}

//remove stale queries
unsigned QueryManager::purgeOldQueries()
{
	boost::unique_lock<boost::mutex>(mu);
	vector<unsigned> toErase;
	for (QueryMap::iterator itr = queries.begin(), end = queries.end(); itr
			!= end; itr++)
	{
		QueryPtr q = itr->second;
		if (q->idle() > timeout || q->cancelled())
		{
			if (q->finished())
			{
				queries[itr->first] = QueryPtr(); //should remove shared ptr reference
				toErase.push_back(itr->first); //this way can bypass iterator invalidation issues
			}
			else
			{
				q->cancel();
			}
		}
	}

	for (unsigned i = 0, n = toErase.size(); i < n; i++)
	{
		queries.erase(toErase[i]);
	}

	return toErase.size();
}
