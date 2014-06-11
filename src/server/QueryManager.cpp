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


//add a query, return zero if unsuccessful
unsigned QueryManager::add(unsigned oldqid, stream_ptr io)
{

	lock_guard<boost::mutex> lock(mu); //protect the query map

	//explicitly remove old query
	if (oldqid > 0 && queries.count(oldqid) > 0)
	{
		QueryPtr oldq = queries[oldqid];
		oldq->cancel();
		queries.erase(oldqid); //removates shared ptr reference from map
	}

	//read receptor info and rotation/translation info, but leave ligand for minimizer to stream
	string str;
	*io >> str;
	if(str != "receptor")
	{
		*io << "ERROR\nNo receptor\n";
		return 0;
	}
	unsigned rsize = 0; //size of receptor string, must be in pdbqt
	*io >> rsize;
	if(rsize == 0)
	{
		*io << "ERROR\nInvalid receptor size\n";
		return 0;
	}

	string recstr(rsize+1, '\0');
	io->read(&recstr[0], rsize);

	//now rotation and translation matrix
	Reorienter reorient;
	reorient.read(*io);

	//absorb newline before ligand data
	getline(*io,str);

	unsigned id = nextID++;
	queries[id] = QueryPtr(new MinimizationQuery(recstr, reorient,io));
	queries[id]->execute(); //don't wait for result
	return id;
}
