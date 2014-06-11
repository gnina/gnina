/*
 * MinimizationQuery.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#include "MinimizationQuery.h"
#include "parse_pdbqt.h"
#include "parsing.h"
#include <sstream>

MinimizationQuery::MinimizationQuery(const string& recstr, const Reorienter& re,
		stream_ptr data) :
		valid(true), stopQuery(false), lastAccessed(time(NULL)), reorient(re), io(
				data)
{
	//create the initial model
	stringstream rec(recstr);
	initm = parse_receptor_pdbqt("rigid.pdbqt", rec);

}

MinimizationQuery::~MinimizationQuery()
{
	for (unsigned i = 0, n = allResults.size(); i < n; i++)
	{
		delete allResults[i];
	}
	allResults.clear();
	processedResults.clear();

}

