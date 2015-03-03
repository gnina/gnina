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

#ifndef VINA_MODEL_H
#define VINA_MODEL_H

#include <boost/optional.hpp> // for context
#include <boost/serialization/optional.hpp>
#include <boost/serialization/utility.hpp>
#include <string>
#include "file.h"
#include "tree.h"
#include "matrix.h"
#include "precalculate.h"
#include "igrid.h"
#include "grid_dim.h"
#include "grid.h"

struct interacting_pair {
	smt t1;
	smt t2;
	sz a;
	sz b;
	interacting_pair(): t1(smina_atom_type::Hydrogen), t2(smina_atom_type::Hydrogen), a(0), b(0) {}
	interacting_pair(smt t1_, smt t2_, sz a_, sz b_) : t1(t1_), t2(t2_), a(a_), b(b_) {}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & t1;
		ar & t2;
		ar & a;
		ar & b;
	}
};

typedef std::vector<interacting_pair> interacting_pairs;

typedef std::pair<std::string, boost::optional<sz> > parsed_line;
typedef std::vector<parsed_line> pdbqtcontext;

// dkoes - as an alternative to pdbqt, this stores information
//in an sdf friendly format
struct sdfcontext {

	struct sdfatom { //atom info
		atmidx index; //this is set after parsing and corresponds to the model's atom index
		//the sdf index is just the index into the atoms array plus one
		char elem[2]; //element symbol, note not necessarily null terminated
		bool inflex; //true if atom is nonmoving atom in flex - which means we need an offset to get to the coordinate

		sdfatom():index(0), inflex(false) { elem[0] = elem[1] = 0;}
		sdfatom(const char* nm): index(0), inflex(false)
		{
			elem[0] = elem[1] = 0;
			if(nm) {
				elem[0] = nm[0];
				if(nm[0]) elem[1] = nm[1];
			}
		}

		template<class Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & elem;
			//do NOT export index since this is not set until model creation
			//same for inflex
		}
	};
	struct sdfbond { //bond connectivity and type
		atmidx a;
		atmidx b;
		unsigned char type;

		sdfbond(): a(0), b(0), type(0) {}
		sdfbond(unsigned a_, unsigned b_, unsigned t): a(a_), b(b_), type(t) {}
		template<class Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & a;
			ar & b;
			ar & type;
		}
	};

	struct sdfprop { //property (CHG or ISO) info
		atmidx atom;
		char type; // 'c' or 'i'
		char value;

		sdfprop(): atom(0), type(0), value(0) {}
		sdfprop(unsigned short atm, char t, char v): atom(atm), type(t), value(v) {}

		template<class Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & atom;
			ar & type;
			ar & value;
		}
	};
	std::string name; //molecule name
	std::vector<sdfatom> atoms; //index should match index into coords
	std::vector<sdfbond> bonds;
	std::vector<sdfprop> properties; //CHG and ISO go here


	void dump(std::ostream& out) const;
	void write(const vecv& coords, sz nummove, std::ostream& out) const; //output sdf with provided coords
	bool valid() const {return atoms.size() > 0; }
	sz size() const {return atoms.size(); }
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & name;
		ar & atoms;
		ar & bonds;
		ar & properties;
	}
};

class appender; //in model.cpp
//dkoes - the context consists of the original molecular data with references
//(eventually) to atom indices so we can re-insert atom coordinates
//typically, molecules are converted to pdbqt, however this isn't the friendliest
//or most efficient data format, so we also support an sdf context
struct context {
	pdbqtcontext pdbqttext;
	sdfcontext sdftext;

	void writePDBQT(const vecv& coords, std::ostream& out) const;
	void writeSDF(const vecv& coords, sz nummove, std::ostream& out) const { sdftext.write(coords, nummove, out); }
	void update(const appender& transform);
	void set(sz pdbqtindex, sz sdfindex, sz atomindex, bool inf = false);

	sz pdbqtsize() const { return pdbqttext.size(); }
	sz sdfsize() const { return sdftext.size(); }

	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & pdbqttext;
		ar & sdftext;
	}
};

struct ligand : public flexible_body, atom_range {
	unsigned degrees_of_freedom; // can be different from the apparent number of rotatable bonds, because of the disabled torsions
	interacting_pairs pairs;
	context cont;
	ligand(): degrees_of_freedom(0) {}
	ligand(const flexible_body& f, unsigned degrees_of_freedom_) : flexible_body(f), atom_range(0, 0), degrees_of_freedom(degrees_of_freedom_) {}
	void set_range();

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & degrees_of_freedom;
		ar & pairs;
		ar & cont;
        ar & boost::serialization::base_object<flexible_body>(*this);
        ar & boost::serialization::base_object<atom_range>(*this);
	}
};

struct residue : public main_branch {
	residue() {} //serialization
	residue(const main_branch& m) : main_branch(m) {}
};

enum distance_type {DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE};
typedef strictly_triangular_matrix<distance_type> distance_type_matrix;

struct non_cache; // forward declaration
struct naive_non_cache; // forward declaration
struct cache; // forward declaration
struct szv_grid; // forward declaration
struct terms; // forward declaration
struct conf_independent_inputs; // forward declaration
struct pdbqt_initializer; // forward declaration - only declared in parse_pdbqt.cpp
struct model_test;

struct model {
	void append(const model& m);

	sz num_movable_atoms() const { return m_num_movable_atoms; }
	sz num_internal_pairs() const;
	sz num_other_pairs() const { return other_pairs.size(); }
	sz num_ligands() const { return ligands.size(); }
	sz num_flex() const { return flex.size(); }
	sz ligand_degrees_of_freedom(sz ligand_number) const { return ligands[ligand_number].degrees_of_freedom; }
	sz ligand_longest_branch(sz ligand_number) const;
	sz ligand_length(sz ligand_number) const;
	void get_movable_atom_types(std::vector<smt>& movingtypes) const;

	void set_name(const std::string& n) { name = n; }
	const std::string& get_name() const { return name; }

	conf_size get_size() const;
	conf get_initial_conf() const; // torsions = 0, orientations = identity, ligand positions = current

	grid_dims movable_atoms_box(fl add_to_each_dimension, fl granularity = 0.375) const;

	void write_flex  (const path& name, const std::string& remark) const { write_context(flex_context, name, remark); }

	void write_flex  (std::ostream& out) const {
		write_context(flex_context, out);
	}

	void write_flex_sdf( std::ostream& out) const {
		flex_context.writeSDF(coords, m_num_movable_atoms, out);
	}
	void dump_flex_sdf( std::ostream& out) const {
		flex_context.sdftext.dump(out);
	}
	void write_ligand(std::ostream& out) const {
		VINA_FOR_IN(i, ligands)
			write_context(ligands[i].cont, out);
	}
	void write_structure(std::ostream& out) const {
		VINA_FOR_IN(i, ligands)
			write_context(ligands[i].cont, out);
		if(num_flex() > 0) // otherwise remark is written in vain
			write_context(flex_context, out);
	}

	//write ligand data as sdf (no flex); return true if successful
	bool write_sdf(std::ostream& out) const {
		if(ligands.size() > 0 && ligands[0].cont.sdftext.valid()) {
			ligands[0].cont.writeSDF(coords,m_num_movable_atoms,out);
			return true;
		}
		return false;
	}
	void write_structure(std::ostream& out, const std::string& remark) const {
		out << remark;
		write_structure(out);
	}
	void write_structure(const path& name) const { ofile out(name); write_structure(out); }
	void write_model(std::ostream& out, sz model_number, const std::string& remark = "") const {
		out << "MODEL " << model_number << '\n';
		write_structure(out, remark);
		out << "ENDMDL\n";
	}
	void seti(const conf& c);
	void sete(const conf& c);
	void set (const conf& c);

	std::string ligand_atom_str(sz i, sz lig=0) const;
	fl gyration_radius(sz ligand_number) const; // uses coords

	const atom_base& movable_atom  (sz i) const { assert(i < m_num_movable_atoms); return  atoms[i]; }
	const vec&       movable_coords(sz i) const { assert(i < m_num_movable_atoms); return coords[i]; }

	const vec& atom_coords(const atom_index& i) const;
	fl distance_sqr_between(const atom_index& a, const atom_index& b) const;
	bool atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const; // there is an atom closer to both a and b then they are to each other and immobile relative to them

	distance_type distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const;

	// clean up
	fl evali     (const precalculate& p,                  const vec& v                          		) const;
	fl evale     (const precalculate& p, const igrid& ig, const vec& v                          		) const;
	fl eval      (const precalculate& p, const igrid& ig, const vec& v, const conf& c, const grid& user_grid	);
	fl eval_deriv(const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g, const grid& user_grid);

	fl eval_flex(const precalculate& p, const vec& v, const conf& c, unsigned maxGridAtom=0);
	fl eval_intramolecular(const precalculate& p, const vec& v, const conf& c);
	fl eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy, const grid& user_grid);


	fl rmsd_lower_bound(const model& m) const; // uses coords
	fl rmsd_upper_bound(const model& m) const; // uses coords
	fl rmsd_ligands_upper_bound(const model& m) const; // uses coords

	void verify_bond_lengths() const;
	void about() const;

	vecv get_ligand_internal_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			tmp.push_back(internal_coords[i]);
		return tmp;
	}

	vecv get_ligand_coords() const { // FIXME rm
		VINA_CHECK(ligands.size() == 1);
		vecv tmp;
		const ligand& lig = ligands.front();
		VINA_RANGE(i, lig.begin, lig.end)
			tmp.push_back(coords[i]);
		return tmp;
	}

	vecv& coordinates() { //return reference to all coords
		return coords;
	}

	void dump_coords(std::ostream& out) const {
		VINA_FOR(i, coords.size()) {
			out << i << " " << coords[i][0] << "," << coords[i][1] << "," << coords[i][2] << "\n";
		}
	}
	vecv get_heavy_atom_movable_coords() const { // FIXME mv
		vecv tmp;
		VINA_FOR(i, num_movable_atoms())
			if(!atoms[i].is_hydrogen())
				tmp.push_back(coords[i]);
		return tmp;
	}
	void check_internal_pairs() const;
	void print_stuff() const; // FIXME rm

	fl clash_penalty() const;

	model() : m_num_movable_atoms(0) {};

private:
	//my, aren't we friendly!
	friend struct non_cache;
	friend struct non_cache_gpu;
	friend struct naive_non_cache;
	friend struct cache;
	friend struct szv_grid;
	friend class szv_grid_cache;
	friend struct terms;
	friend struct conf_independent_inputs;
	friend struct appender_info;
	friend struct pdbqt_initializer;
	friend struct model_test;

	const atom& get_atom(const atom_index& i) const { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }
	      atom& get_atom(const atom_index& i)       { return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]); }

	void write_context(const context& c, std::ostream& out) const;
	void write_context(const context& c, std::ostream& out, const std::string& remark) const {
		out << remark;
	}
	void write_context(const context& c, const path& name) const {
		ofile out(name);
		write_context(c, out);
	}
	void write_context(const context& c, const path& name, const std::string& remark) const {
		ofile out(name);
		write_context(c, out, remark);
	}
	fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const; // actually static
	
	atom_index sz_to_atom_index(sz i) const; // grid_atoms, atoms
	bool bonded_to_HD(const atom& a) const;
	bool bonded_to_heteroatom(const atom& a) const;
	sz find_ligand(sz a) const;
	void bonded_to(sz a, sz n, szv& out) const;
	szv bonded_to(sz a, sz n) const;

	void assign_bonds(const distance_type_matrix& mobility); // assign bonds based on relative mobility, distance and covalent length
	void assign_types();
	void initialize_pairs(const distance_type_matrix& mobility);
	void initialize(const distance_type_matrix& mobility);
	fl clash_penalty_aux(const interacting_pairs& pairs) const;

	fl eval_interacting_pairs(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords) const;
	fl eval_interacting_pairs_deriv(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces) const;

	vecv internal_coords;
	vecv coords;
	vecv minus_forces; //I believe this contains the accumulated directional deltas for each atom

	atomv grid_atoms;
	atomv atoms; // movable, inflex

	vector_mutable<ligand> ligands;
	vector_mutable<residue> flex;
	context flex_context;
	interacting_pairs other_pairs;  // all except internal to one ligand: ligand-other ligands; ligand-flex/inflex; flex-flex/inflex

	sz m_num_movable_atoms;

	std::string name;
};

#endif
