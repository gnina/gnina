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
	atmidx  number;
	parsed_atom(): number(0) {}
	parsed_atom(smt sm_, fl charge_, const vec& coords_, unsigned number_) : number(number_) {
		sm = sm_;
		charge = charge_;
		coords = coords_;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		//ar & number; //not needed after parsing
        ar & boost::serialization::base_object<atom>(*this);
	}
};

struct atom_reference {
	atmidx index;
	bool inflex;

	atom_reference(): index(USHRT_MAX), inflex(0) {}
	atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}

	bool valid() { return index != USHRT_MAX; } //so we can avoid optionals
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & index;
		ar & inflex;
	}
};

struct parsing_struct {
	// start reading after this class
	template<typename T> // T == parsing_struct
		struct node_t {
			atmidx pdbqt_context_index; //index into pdbqt context (lines of pdbqt file for reinsertion of coordinates)
			atmidx sdf_context_index; //index into sdf file (what atom index is represented)
			parsed_atom a;
			std::vector<T> ps;
			node_t(): pdbqt_context_index(0), sdf_context_index(0) {} //for serialization
			node_t(const parsed_atom& a_, sz context_index_, sz sdf_index) : pdbqt_context_index(context_index_), sdf_context_index(sdf_index), a(a_) {}


			// inflex atom insertion
			void insert_inflex(non_rigid_parsed& nr, context& c) {
				VINA_FOR_IN(i, ps)
					ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
				c.set(pdbqt_context_index, sdf_context_index, nr.inflex.size(), true);
				nr.inflex.push_back(a);
			}
			void insert_immobiles_inflex(non_rigid_parsed& nr, context& c) {
				VINA_FOR_IN(i, ps)
					ps[i].insert_immobile_inflex(nr, c);
			}

			// insertion into non_rigid_parsed
			void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
				VINA_FOR_IN(i, ps)
					ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
				vec relative_coords; relative_coords = a.coords - frame_origin;
				c.set(pdbqt_context_index, sdf_context_index, nr.atoms.size());
				nr.atoms.push_back(movable_atom(a, relative_coords));
			}
			void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
				VINA_FOR_IN(i, ps)
					ps[i].insert_immobile(nr, c, frame_origin);
			}
		};

	typedef node_t<parsing_struct> node;
	boost::optional<atmidx> immobile_atom; // which of `atoms' is immobile, if any
	atom_reference axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
	atom_reference axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
	std::vector<node> atoms;

	void add(const parsed_atom& a, const context& c, int sdfcontextpos) {
		atoms.push_back(node(a, c.pdbqtsize()-1, sdfcontextpos));
	}
	void add(const parsed_atom& a, const context& c) {
		add(a, c, 0);
	}
	const vec& immobile_atom_coords() const {
		VINA_CHECK(immobile_atom);
		VINA_CHECK(immobile_atom.get() < atoms.size());
		return atoms[immobile_atom.get()].a.coords;
	}
	// inflex insertion
	void insert_immobile_inflex(non_rigid_parsed& nr, context& c) {
		if(!atoms.empty()) {
			VINA_CHECK(immobile_atom);
			VINA_CHECK(immobile_atom.get() < atoms.size());
			axis_end = atom_reference(nr.inflex.size(), true);
			atoms[immobile_atom.get()].insert_inflex(nr, c);
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
			atoms.push_back(node(from.atoms[i].a, from.atoms[i].pdbqt_context_index, from.atoms[i].sdf_context_index));
		}
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & immobile_atom;
		ar & axis_begin;
		ar & axis_end;
		ar & atoms;
	}

};

//STL containers ignore above (seems like a bug)
template<class Archive>
void serialize(Archive & ar, parsing_struct::node_t<parsing_struct>& node, const unsigned int version)
{
//	ar & node.pdbqt_context_index; //not part of smina format
	ar & node.sdf_context_index;
	ar & node.a;
	ar & node.ps;
}

extern void add_pdbqt_context(context& c, const std::string& str);
extern void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof);
extern void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c);

struct pdbqt_initializer {
	model m;
	void initialize_from_rigid(const rigid& r) { // static really
		VINA_CHECK(m.grid_atoms.empty());
		m.grid_atoms = r.atoms;
	}
	void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
		VINA_CHECK(m.ligands.empty());
		VINA_CHECK(m.flex   .empty());

		m.ligands = nrp.ligands;
		m.flex    = nrp.flex;

		VINA_CHECK(m.atoms.empty());

		sz n = nrp.atoms.size() + nrp.inflex.size();
		m.atoms.reserve(n);
		m.coords.reserve(n);

		VINA_FOR_IN(i, nrp.atoms) {
			const movable_atom& a = nrp.atoms[i];
			atom b = static_cast<atom>(a);
			b.coords = a.relative_coords;
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_FOR_IN(i, nrp.inflex) {
			const atom& a = nrp.inflex[i];
			atom b = a;
			b.coords = zero_vec; // to avoid any confusion; presumably these will never be looked at
			m.atoms.push_back(b);
			m.coords.push_back(a.coords);
		}
		VINA_CHECK(m.coords.size() == n);

		m.internal_coords.resize(m.coords.size(), zero_vec); // FIXME

		m.minus_forces = m.coords;
		m.m_num_movable_atoms = nrp.atoms.size();

		if(is_ligand) {
			VINA_CHECK(m.ligands.size() == 1);
			m.ligands.front().cont = c;
		}
		else
			m.flex_context = c;

	}
	void initialize(const distance_type_matrix& mobility) {
		m.initialize(mobility);
	}
};

namespace boost {
namespace serialization {
//default all our classes to not have version info
template <class T>
struct implementation_level_impl< const T >
{
    template<class U>
    struct traits_class_level {
        typedef BOOST_DEDUCED_TYPENAME U::level type;
    };

    typedef mpl::integral_c_tag tag;

    typedef
        BOOST_DEDUCED_TYPENAME mpl::eval_if<
            is_base_and_derived<boost::serialization::basic_traits, T>,
            traits_class_level< T >,
        //else
        BOOST_DEDUCED_TYPENAME mpl::eval_if<
            is_fundamental< T >,
            mpl::int_<primitive_type>,
        //else
        BOOST_DEDUCED_TYPENAME mpl::eval_if<
            mpl::or_<is_class< T >, is_array< T> >,
            mpl::int_<object_serializable>,
        //else
        BOOST_DEDUCED_TYPENAME mpl::eval_if<
            is_enum< T >,
                mpl::int_<primitive_type>,
        //else
            mpl::int_<not_serializable>
        >
        >
        >
        >::type type;
    BOOST_STATIC_CONSTANT(int, value = type::value);
};

//STL containers ignore above (seems like a bug)
template<class Archive, class T>
void serialize(Archive & ar, std::vector<T>  & v, const unsigned int version)
{
	atmidx sz = v.size(); //try to save some space
	assert(v.size() < USHRT_MAX);
    ar & sz;

    //depending on whether we are storing or loading, sz may change
    v.resize(sz);
    for(unsigned i = 0, n = v.size(); i < n; i++)
    	ar & v[i];
}

}
}

#endif /* PARSING_H_ */
