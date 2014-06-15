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
		*io << "ERROR\nNo receptor\n";
		return 0;
	}
	unsigned rsize = 0; //size of receptor string, must be in pdbqt
	*io >> rsize;
	if (rsize == 0)
	{
		*io << "ERROR\nInvalid receptor size\n";
		return 0;
	}

	string recstr(rsize + 1, '\0');
	io->read(&recstr[0], rsize);

	//does the ligand data have to be reoriented?
	bool hasR = false;
	*io >> hasR;

	//absorb newline before ligand data
	getline(*io, str);

	QueryPtr q;
	//protect access to nextID and queries
	mu.lock();
	unsigned id = nextID++;
	q = QueryPtr(new MinimizationQuery(minparm, recstr, io, hasR));
	queries[id] = q;
	mu.unlock();

	q->execute(); //don't wait for result

	return id;
}

//count types of queries
void QueryManager::getCounts(unsigned& active, unsigned& inactive,
		unsigned& defunct)
{
	active = inactive = defunct = 0;
	lock_guard<mutex> L(mu);

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
	unique_lock<mutex>(mu);

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
	unique_lock<mutex>(mu);
	vector<unsigned> toErase;
	for (QueryMap::iterator itr = queries.begin(), end = queries.end(); itr
			!= end; itr++)
	{
		QueryPtr q = itr->second;
		if (q->idle() > timeout  || q->cancelled())
		{
			if(q->finished())
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

	for(unsigned i = 0, n = toErase.size(); i < n; i++)
	{
		queries.erase(toErase[i]);
	}

	return toErase.size();
}
