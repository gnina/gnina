/*
 * Copyright 2012 David Koes and University of Pittsburgh
 *
 * This implements a customizable scoring function class that can be
 * initialized using a user provided file.
 */

#include "custom_terms.h"
#include "everything.h"
#include <boost/lexical_cast.hpp>

using namespace boost;

custom_terms::custom_terms()
{
	ad4_solvation_re.assign("ad4_solvation\\(d-sigma=(\\S+),_s/q=(\\S+),_q=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	electrostatic_re.assign("electrostatic\\(i=(\\S+),_^=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	gauss_re.assign("gauss\\(o=(\\S+),_w=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	hydrophobic_re.assign("hydrophobic\\(g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	ligand_length_re.assign("ligand_length",boost::regex::perl);
	non_dir_h_bond_re.assign("non_dir_h_bond\\(g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	non_hydrophobic_re.assign("non_hydrophobic\\(g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	num_re.assign("num_(\\S+)",boost::regex::perl);
	repulsion_re.assign("repulsion\\(o=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
	vdw_re.assign("vdw\\(i=(\\S+),_j=(\\S+),_s=(\\S+),_^=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
}

//parse the name of a term and add it with the proper parameterization
void custom_terms::add(const std::string& name, fl weight)
{
	smatch match;
	try
	{

	if(regex_match(name, match, vdw_re))
	{
		fl i = lexical_cast<fl>(match[1]);
		fl j = lexical_cast<fl>(match[2]);
		fl s = lexical_cast<fl>(match[3]);
		fl cap = lexical_cast<fl>(match[4]);
		fl c = lexical_cast<fl>(match[5]);
		if(i == 4.0 && j == 8)
			terms::add(1, new vdw<4,8>(s, cap, c));
		else if(i == 6 && j == 12)
			terms::add(1, new vdw<6,12>(s, cap, c));
		else
			throw scoring_function_error(name,"Unsupported LJ exponents: try <4,8> or <6,12>.");
		usable_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, repulsion_re))
	{
		fl o = lexical_cast<fl>(match[1]);
		fl c = lexical_cast<fl>(match[2]);
		terms::add(1, new repulsion(o, c));
		usable_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, num_re))
	{
		std::string counter = match[1];
		if(counter == "tors_div")
			terms::add(1, new num_tors_div());
		else if(counter == "heavy_atoms_div")
			terms::add(1,new num_heavy_atoms_div());
		else if(counter == "heavy_atoms")
			terms::add(1,new num_heavy_atoms());
		else if(counter == "tors_add")
			terms::add(1,new num_tors_add());
		else if(counter == "tors_sqr")
			terms::add(1,new num_tors_sqr());
		else if(counter == "tors_sqrt")
			terms::add(1,new num_tors_sqrt());
		else if(counter == "hydrophobic_atoms")
			terms::add(1,new num_hydrophobic_atoms());
		else
			throw scoring_function_error(name, "Unknown counter");

		conf_independent_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, non_hydrophobic_re))
	{
		fl g = lexical_cast<fl>(match[1]);
		fl b = lexical_cast<fl>(match[2]);
		fl c = lexical_cast<fl>(match[3]);
		terms::add(1, new non_hydrophobic(g, b, c));
		usable_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, non_dir_h_bond_re))
	{
		fl g = lexical_cast<fl>(match[1]);
		fl b = lexical_cast<fl>(match[2]);
		fl c = lexical_cast<fl>(match[3]);
		terms::add(1, new non_dir_h_bond(g, b, c));
		usable_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, ligand_length_re))
	{
		terms::add(1, new ligand_length());
		conf_independent_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, hydrophobic_re))
	{
		fl g = lexical_cast<fl>(match[1]);
		fl b = lexical_cast<fl>(match[2]);
		fl c = lexical_cast<fl>(match[3]);
		terms::add(1, new hydrophobic(g, b, c));
		usable_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, gauss_re))
	{
		fl o = lexical_cast<fl>(match[1]);
		fl w = lexical_cast<fl>(match[2]);
		fl c = lexical_cast<fl>(match[3]);
		terms::add(1, new gauss(o, w, c));
		usable_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, electrostatic_re))
	{
		fl i = lexical_cast<fl>(match[1]);
		fl cap = lexical_cast<fl>(match[2]);
		fl c = lexical_cast<fl>(match[3]);
		if(i == 1)
			terms::add(1, new electrostatic<1>(cap,c));
		else if(i == 2)
			terms::add(1, new electrostatic<2>(cap,c));
		else
			throw scoring_function_error(name,"Invalid exponent: 1 or 2 only");

		distance_additive_weights.push_back(weight);
		return;
	}
	if(regex_match(name, match, ad4_solvation_re))
	{
		fl sigma = lexical_cast<fl>(match[1]);
		fl w = lexical_cast<fl>(match[2]);
		bool q = lexical_cast<bool>(match[3]);
		fl c = lexical_cast<fl>(match[4]);
		terms::add(1, new ad4_solvation(sigma, w,  q, c));
		distance_additive_weights.push_back(weight);
		return;
	}

	}
	catch(bad_lexical_cast be)
	{
		throw scoring_function_error(name,"Could not convert parameters. ");
	}
}

//return weights in correct order
flv custom_terms::weights() const
{
	flv ret;
	ret.insert(ret.end(), usable_weights.begin(), usable_weights.end());
	ret.insert(ret.end(), distance_additive_weights.begin(), distance_additive_weights.end());
	ret.insert(ret.end(), additive_weights.begin(), additive_weights.end());
	ret.insert(ret.end(), intermolecular_weights.begin(), intermolecular_weights.end());
	ret.insert(ret.end(), conf_independent_weights.begin(), conf_independent_weights.end());

	return ret;
}
