/*
 * MinimizationQuery.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 */

#ifndef MINIMIZATIONQUERY_H_
#define MINIMIZATIONQUERY_H_

#include <vector>

#include "Reorienter.h"
#include "server_common.h"
#include "model.h"

class MinimizationQuery
{
	bool valid;
	bool stopQuery;
	time_t lastAccessed;

	Reorienter reorient;
	stream_ptr io;

	model initm;

	//holds the result of minimization
	struct Result {
		double score;
		double rmsd;
		string name;
		string sdf;
	};

	vector<Result*> processedResults; //after sorting and filtering, what the user sees
	vector<Result*> allResults; //order doesn't change, minimizers add to this
public:

	//criteria for filtering and (maybe eventually) sorting the data
	struct Filters {
		double maxScore;
		double maxRMSD;

		Filters(): maxScore(HUGE_VAL), maxRMSD(HUGE_VAL) {}
	};



	MinimizationQuery(const string& recstr, const Reorienter& re, stream_ptr data);

	~MinimizationQuery();

	bool isValid() const
	{
		return valid;
	}

	void execute(bool block = true);

	//all of the result/output functions can be called while an asynchronous
	//query is running

	//return all current results

	void outputData(const Filters& dp, ostream& out);
	//write out all results in sdf.gz format
	void outputMols(const Filters& dp, ostream& out);
	//output single mol in sdf format; referenced using current position in processed results array
	void outputMol(unsigned index, ostream& out);

	//attempt to cancel,
	void cancel();
	bool finished(); //okay to deallocate, user may still care though
	bool cancelled() { return stopQuery; } //user no longer cares
	void access()
	{
		lastAccessed = time(NULL);
		stopQuery = false;
	}
	const time_t idle()
	{
		return time(NULL) - lastAccessed;
	} //time since last access


};

#endif /* MINIMIZATIONQUERY_H_ */
