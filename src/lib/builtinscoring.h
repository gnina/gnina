/*
 * builtinscoring.h
 *
 *  Created on: Oct 15, 2015
 *      Author: dkoes
 *
 *  Singleton class for setting up builtin scoring terms
 */

#include <vector>
#include <boost/unordered_map.hpp>
#include "custom_terms.h"

class builtin_scoring_function
{
	struct singleterm
	{
		const char* term;
		double weight;

		singleterm(): term(NULL), weight(0) {}
		singleterm(const char* t, double w): term(t), weight(w) {}
	};

	//map from names to sets of terms
	boost::unordered_map<std::string, std::vector<singleterm> > functions;

	void add(const std::string& name, const char *term, double w)
	{
		functions[name].push_back(singleterm(term,w));
	}

public:
	builtin_scoring_function();

	//set t to the terms specified by builtin name
	//return true if successful
	bool set(custom_terms& t, const std::string& name)
	{
		if(functions.count(name) == 0)
			return false;

		std::vector<singleterm>& terms = functions[name];
		for(unsigned i = 0, n = terms.size(); i < n; i++)
		{
			t.add(terms[i].term, terms[i].weight);
		}
		return true;
	}
};

extern builtin_scoring_function builtin_scoring_functions;
