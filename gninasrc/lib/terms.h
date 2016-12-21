/*

 Copyright (c) 2006-2010, The Scripps Research Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 Author: Dr. Oleg Trott <ot14@columbia.edu>,
 The Olson Lab,
 The Scripps Research Institute

 */

#ifndef VINA_TERMS_H
#define VINA_TERMS_H

#include <boost/ptr_container/ptr_vector.hpp> 
#include <boost/regex.hpp>
#include "model.h"
#include "result_components.h"

//thrown when can't parse name of term
struct scoring_function_error {
	std::string name;
	std::string msg;
	scoring_function_error(const std::string& n, const std::string& m="") : name(n), msg(m) {}
};


enum TermKind {BaseTerm, ConfIndependent, InterMolecular, Additive,  DistanceAdditive, ChargeDependent, ChargeIndependent, LastTermKind};

struct term
{
	boost::regex rexpr;
	std::string name;
	virtual ~term()
	{
	}

	//every term must be able to create a new term from string description
	//of the paramerterized term
	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return BaseTerm;
	}
};


struct distance_additive: public term
{
	fl cutoff;

	distance_additive(): cutoff(0) {}
	distance_additive(fl cutoff_) :
			cutoff(cutoff_)
	{
	}
	virtual fl eval(const atom_base& a, const atom_base& b, fl r) const = 0;
	virtual ~distance_additive()
	{
	}

	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return DistanceAdditive;
	}
};


struct charge_dependent: public distance_additive
{
	charge_dependent() {}
	charge_dependent(fl cut): distance_additive(cut) {}
	virtual ~charge_dependent() {}

	//unique to charge_dependent, return comonents for given types and distance
	virtual result_components eval_components(smt t1, smt t2, fl r) const = 0;

	fl eval(const atom_base& a, const atom_base& b, fl r) const
	{
		result_components c = eval_components(a.sm, b.sm, r);
		return c.eval(a,b);
	}

	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return ChargeDependent;
	}

};

//these terms can be fully precalculated just from atom types
struct charge_independent: public distance_additive
{
	charge_independent() {}
	charge_independent(fl cutoff_) :
			distance_additive(cutoff_)
	{
	}
	fl eval(const atom_base& a, const atom_base& b, fl r) const
	{
		return eval(a.get(), b.get(), r);
	}
	virtual fl eval(smt t1, smt t2, fl r) const
	{
		VINA_CHECK(false);
		return 0;
	}
	virtual ~charge_independent()
	{
	}

	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return ChargeIndependent;
	}
};

struct additive: public term
{
	fl cutoff;
	additive() :
			cutoff(max_fl)
	{
	}
	virtual fl eval(const model& m, const atom_index& i,
			const atom_index& j) const = 0;
	virtual ~additive()
	{
	}

	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return Additive;
	}

};

struct intermolecular: public term
{
	virtual fl eval(const model& m) const = 0;

	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return InterMolecular;
	}

};

struct conf_independent_inputs
{
	fl num_tors;
	fl num_rotors;
	fl num_heavy_atoms;
	fl num_hydrophobic_atoms;
	fl ligand_max_num_h_bonds;
	fl num_ligands;
	fl ligand_lengths_sum;
	operator flv() const;
	conf_independent_inputs(const model& m);
	std::vector<std::string> get_names() const;
	conf_independent_inputs();
	private:
	unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i) const; // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned atom_rotors(const model& m, const atom_index& i) const; // the number of rotatable bonds to heavy ligand atoms

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version)
	{
		ar & num_tors;
		ar & num_rotors;
		ar & num_heavy_atoms;
		ar & num_hydrophobic_atoms;
		ar & ligand_max_num_h_bonds;
		ar & num_ligands;
		ar & ligand_lengths_sum;
	}
};

struct conf_independent: public term
{
	virtual fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& it) const = 0;
	virtual sz size() const = 0; // how many parameters does it take

	//assume these are unparamerterized and so just have to match regex
	virtual term* createFrom(const std::string& name) const = 0;

	virtual TermKind kind() const {
		return ConfIndependent;
	}

};

template<typename T>
struct term_set
{
	std::vector<bool> enabled;
	boost::ptr_vector<T> fun; // FIXME? const T?
	void add(unsigned e, T* f)
	{ // FIXME? const T* ?
		enabled.push_back(e > 0);
		fun.push_back(f);
	}
	sz num_enabled() const
	{
		sz tmp = 0;
		VINA_FOR_IN(i, enabled)
		if(enabled[i])
		++tmp;
		return tmp;
	}
	void get_names(bool enabled_only, std::vector<std::string>& out) const
	{ // appends to "out"
		VINA_CHECK(enabled.size() == fun.size());
		VINA_FOR_IN(i, fun)
		if(!enabled_only || enabled[i])
		out.push_back(fun[i].name);
	}
	void filter(flv::const_iterator& in, flv& out) const
	{
		VINA_CHECK(enabled.size() == fun.size());
		VINA_FOR_IN(i, enabled)
		{
			if(enabled[i])
			out.push_back(*in);
			++in;
		}
	}
	fl max_cutoff() const
	{
		fl tmp = 0;
		VINA_FOR_IN(i, fun)
		tmp = (std::max)(tmp, fun[i].cutoff);
		return tmp;
	}
	sz size() const
	{	return fun.size();}
	const T& operator[](sz i) const
	{	return fun[i];}
};

struct factors
{
	flv e; // external
	flv i; // internal
	sz size() const
	{
		return e.size() + i.size();
	}
	//sz num_weights() const { return (std::max)(e.size(), i.size()); } // FIXME? compiler bug? getting warnings here
	sz num_weights() const
	{
		return (e.size() > i.size()) ? e.size() : i.size();
	}
	fl eval(const flv& weights, bool include_internal) const;
	private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version)
	{
		ar & e;
		ar & i;
	}
};

struct terms
{
	term_set<charge_independent> charge_independent_terms;
	term_set<charge_dependent> charge_dependent_terms;
	term_set<distance_additive> distance_additive_terms;
	term_set<additive> additive_terms;
	term_set<intermolecular> intermolecular_terms;
	term_set<conf_independent> conf_independent_terms;

	// the class takes ownership of the pointer with 'add'
	void add(unsigned e, charge_independent* p)
	{
		charge_independent_terms.add(e, p);
	}
	void add(unsigned e, charge_dependent* p)
	{
		charge_dependent_terms.add(e, p);
	}
	void add(unsigned e, distance_additive* p)
	{
		distance_additive_terms.add(e, p);
	}
	void add(unsigned e, additive* p)
	{
		additive_terms.add(e, p);
	}
	void add(unsigned e, intermolecular* p)
	{
		intermolecular_terms.add(e, p);
	}
	void add(unsigned e, conf_independent* p)
	{
		conf_independent_terms.add(e, p);
	}

	void add(unsigned e, term *p)
	{
		switch(p->kind())
		{
		case ChargeIndependent:
			charge_independent_terms.add(e, dynamic_cast<charge_independent*>(p));
			break;
		case ChargeDependent:
			charge_dependent_terms.add(e, dynamic_cast<charge_dependent*>(p));
			break;
		case DistanceAdditive:
			distance_additive_terms.add(e, dynamic_cast<distance_additive*>(p));
			break;
		case Additive:
			additive_terms.add(e, dynamic_cast<additive*>(p));
			break;
		case InterMolecular:
			intermolecular_terms.add(e,dynamic_cast<intermolecular*>(p));
			break;
		case ConfIndependent:
			conf_independent_terms.add(e, dynamic_cast<conf_independent*>(p));
			break;
		default:
			VINA_CHECK(false);
			break;
		}
	}

	std::vector<std::string> get_names(bool enabled_only) const; // does not include conf-independent
	sz size_internal() const;
	sz size() const
	{
		return size_internal() + intermolecular_terms.size();
	}
	sz size_conf_independent(bool enabled_only) const; // number of parameters does not necessarily equal the number of operators
	fl max_r_cutoff() const;
	flv evale_robust(const model& m, std::vector<flv>& per_atom) const;
	flv evale_robust(const model& m) const { std::vector<flv> pa; return evale_robust(m, pa); }
	fl eval_conf_independent(const conf_independent_inputs& in, fl x,
			flv::const_iterator& it) const;
	flv filter_external(const flv& v) const;
	flv filter_internal(const flv& v) const;
	factors filter(const factors& f) const;
	void display_info() const;
	private:
	void eval_additive_aux(const model& m, const atom_index& i,
			const atom_index& j, fl r, flv& out) const; // out is added to

};

#endif
