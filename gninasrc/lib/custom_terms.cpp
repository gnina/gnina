/*
 * Copyright 2012 David Koes and University of Pittsburgh
 *
 * This implements a customizable scoring function class that can be
 * initialized using a user provided file.
 */

#include "custom_terms.h"
#include "everything.h"
#include <boost/lexical_cast.hpp>
#include <algorithm>

using namespace boost;


//parse the name of a term and add it with the proper parameterization
void custom_terms::add(const std::string& name, fl weight)
{
	try
	{
		for(unsigned i = 0, n = creators.size(); i < n; i++)
		{
			term *t = creators[i]->createFrom(name);
			if(t != NULL)
			{
				TermKind kind = t->kind();
				term_weights[(unsigned)kind].push_back(weight * custom_terms::scaling_factor);
				terms::add(1, t);
				return;
			}
		}
	}
	catch(bad_lexical_cast& be)
	{
		throw scoring_function_error(name,"Could not convert parameters. ");
	}
	throw scoring_function_error(name, "Unknown term ");
}

//return weights in correct order
flv custom_terms::weights() const
{
	flv ret;
	for(int index = LastTermKind-1; index > BaseTerm; index--)
	{
		//revers order is what is expected
		ret.insert(ret.end(), term_weights[index].begin(), term_weights[index].end());
	}

	return ret;
}

//add weights and terms from file,
//each line has the weight then the term name, whitespace separated, anything
//after the term name is ignored
void custom_terms::add_terms_from_file(std::istream& in)
{
	std::string line;
	while(getline(in, line))
	{
		if(line.size() < 1 || line[0] == '#')
			continue; //support comments and blank lines
		std::stringstream str(line);
		fl w =0 ;
		std::string name;
		if((str >> w >> name) && name.length() > 0)
		{
			add(name, w);
		}
	}
}

//print out terms in format that can be read back in
void custom_terms::print(std::ostream& out) const
{
	flv ret;
	unsigned pad = 12;

	for(int index = LastTermKind-1; index > ConfIndependent; index--)
	{
		//revers order is what is expected
		ret.insert(ret.end(), term_weights[index].begin(), term_weights[index].end());
	}

	std::vector<std::string> names = get_names(true);
	assert(ret.size() == names.size());
	for(unsigned i = 0, n = ret.size(); i < n; i++)
	{
		out << std::setw(pad) << ret[i] << " " << names[i] << "\n";
	}

	//now conf_indep, which are separate for some reason
	std::vector<std::string> conf_indep_names;
	conf_independent_terms.get_names(true, conf_indep_names);
	assert(conf_indep_names.size() == term_weights[ConfIndependent].size());
	for(unsigned i = 0, n = conf_indep_names.size(); i < n; i++)
	{
		out << std::setw(pad) << term_weights[ConfIndependent][i] << " " << conf_indep_names[i] << "\n";
	}
}

//print all the term creators, not those created
void custom_terms::print_available_terms(std::ostream& out) const
{
	for(unsigned i = 0, n = creators.size(); i < n; i++)
	{
		out << creators[i]->name << "\n";
	}
}


std::ostream& operator<<(std::ostream& out, const custom_terms& t)
{
	t.print(out);
	return out;
}

void custom_terms::set_scaling_factor(fl sf){
	scaling_factor = sf;
}
