/* GPL 2.0
 * Copyright 2012 David Koes and University of Pittsburgh
 *
 * This implements a customizable scoring function class that can be
 * initialized using a user provided file.
 */
#ifndef SMINA_CUSTOM_TERMS_H
#define SMINA_CUSTOM_TERMS_H

#include "terms.h"
#include "int_pow.h"
#include "everything.h"
#include <string>
#include <boost/regex.hpp>
#include "file.h"

//thrown when can't parse name of term
struct scoring_function_error {
	std::string name;
	std::string msg;
	scoring_function_error(const std::string& n, const std::string& m="") : name(n), msg(m) {}
};

//terms that can be dynamically built
//also, keeps track of weights (in right order)
struct custom_terms: public terms
{
private:
	//weights must match type-specific ordering of terms
	flv distance_additive_weights;
	flv usable_weights;
	flv additive_weights;
	flv intermolecular_weights;
	flv conf_independent_weights;

	//construct regular expressions for matchiner term names once
	boost::regex ad4_solvation_re;
	boost::regex constant_re;
	boost::regex electrostatic_re;
	boost::regex gauss_re;
	boost::regex hydrophobic_re;
	boost::regex ligand_length_re;
	boost::regex non_dir_h_bond_re;
	boost::regex non_dir_h_bond_quadratic_re;
	boost::regex non_dir_h_bond_lj_re;
	boost::regex non_hydrophobic_re;
	boost::regex num_re; //counts
	boost::regex repulsion_re;
	boost::regex vdw_re;

public:
	custom_terms();//creates empty set, inits regexes
	void add(const std::string& name, fl weight);
	flv weights() const;
	void add_terms_from_file(std::istream& in);
	void print(std::ostream& out) const;

	friend std::ostream& operator<<(std::ostream& out, const custom_terms& t);
};

#endif
