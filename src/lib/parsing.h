/*
 * parsing.h
 *
 *  Created on: Jun 4, 2014
 *      Author: dkoes
 *
 *  Export various pdbqt parsing routines.
 */

#ifndef PARSING_H_
#define PARSING_H_

#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <boost/utility.hpp> // for noncopyable
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include "parse_pdbqt.h"
#include "atom_constants.h"
#include "convert_substring.h"
#include "parse_error.h"

struct movable_atom : public atom {
	vec relative_coords;
	movable_atom() {}
	movable_atom(const atom& a, const vec& relative_coords_) : atom(a) {
		relative_coords = relative_coords_;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & relative_coords;
        ar & boost::serialization::base_object<atom>(*this);
	}
};

struct rigid {
	atomv atoms;
};

typedef std::vector<movable_atom> mav;

struct non_rigid_parsed {
	vector_mutable<ligand> ligands;
	vector_mutable<residue> flex;

	mav atoms;
	atomv inflex;

	distance_type_matrix atoms_atoms_bonds;
	matrix<distance_type> atoms_inflex_bonds;
	distance_type_matrix inflex_inflex_bonds;

	distance_type_matrix mobility_matrix() const {
		distance_type_matrix tmp(atoms_atoms_bonds);
		tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
		return tmp;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & ligands;
		ar & flex;
		ar & atoms;
		ar & inflex;
		ar & atoms_atoms_bonds;
		ar & atoms_inflex_bonds;
		ar & inflex_inflex_bonds;
	}

};

struct parsed_atom : public atom {
	unsigned number;
	parsed_atom(smt sm_, fl charge_, const vec& coords_, unsigned number_) : number(number_) {
		sm = sm_;
		charge = charge_;
		coords = coords_;
	}
};

struct atom_reference {
	sz index;
	bool inflex;
	atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

struct parsing_struct {
	// start reading after this class
	template<typename T> // T == parsing_struct
		struct node_t {
			sz context_index;
			parsed_atom a;
			std::vector<T> ps;
			node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

			// inflex atom insertion
			void insert_inflex(non_rigid_parsed& nr) {
				VINA_FOR_IN(i, ps)
					ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
				nr.inflex.push_back(a);
			}
			void insert_immobiles_inflex(non_rigid_parsed& nr) {
				VINA_FOR_IN(i, ps)
					ps[i].insert_immobile_inflex(nr);
			}

			// insertion into non_rigid_parsed
			void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
				VINA_FOR_IN(i, ps)
					ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
				vec relative_coords; relative_coords = a.coords - frame_origin;
				c[context_index].second = nr.atoms.size();
				nr.atoms.push_back(movable_atom(a, relative_coords));
			}
			void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
				VINA_FOR_IN(i, ps)
					ps[i].insert_immobile(nr, c, frame_origin);
			}
		};

	typedef node_t<parsing_struct> node;
	boost::optional<sz> immobile_atom; // which of `atoms' is immobile, if any
	boost::optional<atom_reference> axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
	boost::optional<atom_reference> axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
	std::vector<node> atoms;

	void add(const parsed_atom& a, const context& c) {
		VINA_CHECK(c.size() > 0);
		atoms.push_back(node(a, c.size()-1));
	}
	const vec& immobile_atom_coords() const {
		VINA_CHECK(immobile_atom);
		VINA_CHECK(immobile_atom.get() < atoms.size());
		return atoms[immobile_atom.get()].a.coords;
	}
	// inflex insertion
	void insert_immobile_inflex(non_rigid_parsed& nr) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.inflex.size(), true);
			atoms[immobile_atom.get()].insert_inflex(nr);
		}
	}

	// insertion into non_rigid_parsed
	void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.atoms.size(), false);
			atoms[immobile_atom.get()].insert(nr, c, frame_origin);
		}
	}

	bool essentially_empty() const { // no sub-branches besides immobile atom, including sub-sub-branches, etc
		VINA_FOR_IN(i, atoms) {
			if(immobile_atom && immobile_atom.get() != i)
				return false;
			const node& nd = atoms[i];
			if(!nd.ps.empty())
				return false; // FIXME : iffy
		}
		return true;
	}

	//return true if all the mobile atoms in this branch are hydrogens
	bool mobile_hydrogens_only() const {
		if(!get_fixed_rotable_hydrogens())
			return false;
		VINA_FOR_IN(i, atoms) {
			if(!atoms[i].ps.empty())
				return false; //must be terminal
			if(immobile_atom && immobile_atom.get() != i)
			{
				//mobile atom
				if(!atoms[i].a.is_hydrogen())
					return false;
			}
		}
		return true;
	}

	//move the atoms of from into this
	void mergeInto(const parsing_struct& from)
	{
		VINA_FOR_IN(i, from.atoms)
		{
			atoms.push_back(node(from.atoms[i].a, from.atoms[i].context_index));
		}
	}
};

extern void add_context(context& c, const std::string& str);

#endif /* PARSING_H_ */
