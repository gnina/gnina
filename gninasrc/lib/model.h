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
#include <boost/serialization/utility.hpp>
#include <openbabel/mol.h>
#include <string>
#include "optional_serialization.h"
#include "file.h"
#include "tree.h"
#include "tree_gpu.h"
#include "matrix.h"
#include "precalculate.h"
#include "igrid.h"
#include "grid_dim.h"
#include "grid.h"
#include "gpucode.h"
#include "interacting_pairs.h"
#include "user_opts.h"

typedef std::vector<interacting_pair> interacting_pairs;

typedef std::pair<std::string, boost::optional<sz> > parsed_line;
typedef std::vector<parsed_line> pdbqtcontext;
struct parallel_mc_task;
void test_eval_intra();
struct parallel_mc_task;

struct gpu_data {
    //put all pointers to gpu data into a struct that provides a lightweight 
    //model-like object for the GPU
    atom_params *coords;
    vec *atom_coords;
    force_energy_tup *minus_forces;
    tree_gpu *treegpu;
    interacting_pair *interacting_pairs;
    // all except internal to one ligand: ligand-other ligands;
    // ligand-flex/inflex; flex-flex/inflex
    interacting_pair *other_pairs;
    size_t* dfs_order_bfs_indices;
    size_t* bfs_order_dfs_indices;
    float *scratch; //single value for returning total energy

    unsigned coords_size;
    unsigned atom_coords_size;
    unsigned forces_size;
    unsigned pairs_size;
    unsigned other_pairs_size;
    bool device_on; //use gpu for docking, not just CNN
    int device_id;

    //TODO delete
    size_t nlig_roots;

    gpu_data()
        : coords(NULL), atom_coords(NULL), minus_forces(NULL), treegpu(NULL),
            interacting_pairs(NULL), other_pairs(NULL),
            dfs_order_bfs_indices(NULL), bfs_order_dfs_indices(NULL),
            scratch(NULL), coords_size(0), atom_coords_size(0), forces_size(0),
            pairs_size(0), other_pairs_size(0), device_on(false), device_id(0) {
    }

    template<typename infoT>
    __host__  __device__ fl eval_interacting_pairs_deriv_gpu(const infoT& info,
        fl v, interacting_pair* pairs, unsigned pairs_sz) const;

    template<typename infoT> __device__ fl eval_deriv_gpu(const infoT& info,
        const vec& v, const conf_gpu& c, change_gpu& g);

    size_t node_idx_dfs2bfs(const size_t node_idx);
    fl eval(const GPUNonCacheInfo& info, const float v);
    fl eval_intramolecular(const GPUNonCacheInfo& info, const float v);
    //copy relevant data to gpu buffers
    void copy_to_gpu(model& m);
    //copy back relevant data from gpu buffers
    void copy_from_gpu(model& m);

    size_t node_idx_cpu2gpu(size_t cpu_idx) const;

  private:
    gpu_data(const gpu_data&) = default;
    template<typename T> friend struct quasi_newton_aux_gpu;
    friend parallel_mc_task;
};

// dkoes - as an alternative to pdbqt, this stores information
//in an sdf friendly format
struct sdfcontext {

    struct sdfatom { //atom info
        atmidx index = 0; //this is set after parsing and corresponds to the model's atom index
                      //the sdf index is just the index into the atoms array plus one
        char elem[2] = {0,0}; //element symbol, note not necessarily null terminated
        bool inflex = false; //true if atom is nonmoving atom in flex - which means we
                     //need an offset to get to the coordinate
        bool iscovlig = false; //true if atom is part of a covalent ligand

        sdfatom()
            : index(0), inflex(false) {
          elem[0] = elem[1] = 0;
        }
        sdfatom(const char* nm, bool iscov = false)
            : index(0), inflex(false), iscovlig(iscov) {
          elem[0] = elem[1] = 0;
          if (nm) {
            elem[0] = nm[0];
            if (nm[0]) elem[1] = nm[1];
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

        sdfbond()
            : a(0), b(0), type(0) {
        }
        sdfbond(unsigned a_, unsigned b_, unsigned t)
            : a(a_), b(b_), type(t) {
        }
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

        sdfprop()
            : atom(0), type(0), value(0) {
        }
        sdfprop(unsigned short atm, char t, char v)
            : atom(atm), type(t), value(v) {
        }

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
    std::string datastr; //retained properties

    void dump(std::ostream& out) const;
    //output sdf with provided coords
    void write(const vecv& coords, std::ostream& out, bool covonly=false) const;
    bool valid() const {
      return atoms.size() > 0;
    }
    sz size() const {
      return atoms.size();
    }

    void set_inflex_indices(sz nummove);

    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & name;
      ar & atoms;
      ar & bonds;
      ar & properties;
      ar & datastr;
    }
};

class appender;
//in model.cpp
//dkoes - the context consists of the original molecular data with references
//(eventually) to atom indices so we can re-insert atom coordinates
//typically, molecules are converted to pdbqt, however this isn't the friendliest
//or most efficient data format, so we also support an sdf context
struct context {
    pdbqtcontext pdbqttext;
    sdfcontext sdftext;
    bool has_cov_lig = false;

    void writePDBQT(const vecv& coords, std::ostream& out) const;
    void writeSDF(const vecv& coords, std::ostream& out, bool covonly=false) const {
      sdftext.write(coords, out, covonly);
    }
    void update(const appender& transform);
    void set(sz pdbqtindex, sz sdfindex, sz atomindex, bool inf = false);

    sz pdbqtsize() const {
      return pdbqttext.size();
    }
    sz sdfsize() const {
      return sdftext.size();
    }

    void set_inflex_indices(sz nummove) {
      sdftext.set_inflex_indices(nummove);
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & pdbqttext;
      ar & sdftext;
    }
};

struct ligand : public flexible_body, atom_range {
    // can be different from the apparent number of rotatable bonds,
    // because of the disabled torsions
    unsigned degrees_of_freedom;
    interacting_pairs pairs;
    context cont;
    ligand()
        : degrees_of_freedom(0) {
    }
    ligand(const flexible_body& f, unsigned degrees_of_freedom_)
        : flexible_body(f), atom_range(0, 0),
            degrees_of_freedom(degrees_of_freedom_) {
    }
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
    residue() {
    } //serialization
    residue(const main_branch& m)
        : main_branch(m) {
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & boost::serialization::base_object<main_branch>(*this);
    }
};

enum distance_type {
  DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE
};
typedef strictly_triangular_matrix<distance_type> distance_type_matrix;

struct non_cache;
// forward declaration
struct naive_non_cache;
// forward declaration
struct cache;
// forward declaration
struct szv_grid;
// forward declaration
struct terms;
// forward declaration
struct conf_independent_inputs;
// forward declaration
struct pdbqt_initializer;
// forward declaration - only declared in parse_pdbqt.cpp
struct model_test;

struct model {

    model()
        : m_num_movable_atoms(0), hydrogens_stripped(false) {
    }
    ;
    ~model() {
    }
    ;

    model(const model& m)
        : tree_width(m.tree_width), coords(m.coords), 
            extra_box_coords(m.extra_box_coords), ligands(m.ligands),
            minus_forces(m.minus_forces),
            m_num_movable_atoms(m.m_num_movable_atoms), atoms(m.atoms),
            grid_atoms(m.grid_atoms), other_pairs(m.other_pairs),
            hydrogens_stripped(m.hydrogens_stripped),
            internal_coords(m.internal_coords), flex(m.flex),
            flex_context(m.flex_context), name(m.name), pose_num(m.pose_num) {
    }

    void append(const model& m);
    void strip_hydrogens();

    sz num_movable_atoms() const {
      return m_num_movable_atoms;
    }
    sz num_internal_pairs() const;
    sz num_other_pairs() const {
      return other_pairs.size();
    }
    sz num_ligands() const {
      return ligands.size();
    }
    sz num_flex() const {
      return flex.size();
    }
    sz ligand_degrees_of_freedom(sz ligand_number) const {
      return ligands[ligand_number].degrees_of_freedom;
    }
    sz ligand_longest_branch(sz ligand_number) const;
    sz ligand_length(sz ligand_number) const;
    unsigned tree_width;
    void get_movable_atom_types(std::vector<smt>& movingtypes) const;

    void set_name(const std::string& n) {
      name = n;
    }
    const std::string& get_name() const {
      return name;
    }

    void set_pose_num(int n) {
      pose_num = n;
    }
    int get_pose_num() const {
      return pose_num;
    }

    conf_size get_size() const;
    // torsions = 0, orientations = identity, ligand positions = current
    conf get_initial_conf(bool enable_receptor) const;

    grid_dims movable_atoms_box(fl add_to_each_dimension,
        fl granularity = 0.375) const;

    void write_flex(std::ostream& out, bool& issdf, bool covonly = false) const;

    void dump_flex_sdf(std::ostream& out) const {
      flex_context.sdftext.dump(out);
    }

    bool flex_has_covalent() const {
      return flex_context.has_cov_lig;
    }

    void write_ligand(std::ostream& out, bool& issdf) const {
      issdf = false;
      VINA_FOR_IN(i, ligands) {
        if (ligands[i].cont.sdftext.valid()) {
          ligands[0].cont.writeSDF(coords, out);
          issdf = true;
        } else {
          write_context(ligands[i].cont, out);
        }
      }
    }

    void write_rigid_xyz(std::ostream& out, const vec& center) const;
    void write_structure(std::ostream& out) const {
      VINA_FOR_IN(i, ligands)
        write_context(ligands[i].cont, out);
      if (num_flex() > 0) // otherwise remark is written in vain
        write_context(flex_context, out);
    }

    void write_structure(std::ostream& out, const std::string& remark) const {
      out << remark;
      write_structure(out);
    }
    void write_structure(const path& name) const {
      ofile out(name);
      write_structure(out);
    }
    void write_model(std::ostream& out, sz model_number,
        const std::string& remark = "") const {
      out << "MODEL " << model_number << '\n';
      write_structure(out, remark);
      out << "ENDMDL\n";
    }
    void seti(const conf& c);
    void sete(const conf& c);
    void set(const conf& c);

    std::string ligand_atom_str(sz i, sz lig = 0) const;
    fl gyration_radius(sz ligand_number) const; // uses coords
    fl max_span(sz ligand_number) const; // uses coords

    const atom_base& movable_atom(sz i) const {
      assert(i < m_num_movable_atoms);
      return atoms[i];
    }
    const vec& movable_coords(sz i) const {
      assert(i < m_num_movable_atoms);
      return coords[i];
    }
    vec& movable_minus_forces(sz i) {
      assert(i < m_num_movable_atoms);
      return minus_forces[i];
    }

    const vec& atom_coords(const atom_index& i) const;
    fl distance_sqr_between(const atom_index& a, const atom_index& b) const;
    // there is an atom closer to both a and b then they are to each other
    // and immobile relative to them
    bool atom_exists_between(const distance_type_matrix& mobility,
        const atom_index& a, const atom_index& b,
        const szv& relevant_atoms) const;

    distance_type distance_type_between(const distance_type_matrix& mobility,
        const atom_index& i, const atom_index& j) const;

    // clean up
    fl evali(const precalculate& p, const vec& v);
    fl evale(const precalculate& p, const igrid& ig, const vec& v);
    fl eval(const precalculate& p, const igrid& ig, const vec& v, const conf& c,
        const grid& user_grid);
    fl eval_deriv(const precalculate& p, const igrid& ig, const vec& v,
        const conf& c, change& g, const grid& user_grid);
    fl eval_intra(const precalculate& p, const vec& v);
    fl eval_flex(const precalculate& p, const vec& v, const conf& c,
        unsigned maxGridAtom = 0);
    fl eval_intramolecular(const precalculate& p, const vec& v, const conf& c);
    fl eval_adjusted(const scoring_function& sf, const precalculate& p,
        const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy,
        const grid& user_grid);

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

    const vecv& coordinates() const { //return reference to all coords
      return coords;
    }

    void dump_coords(std::ostream& out) const {
      VINA_FOR(i, coords.size()) {
        out << i << " " << coords[i][0] << "," << coords[i][1] << ","
            << coords[i][2] << "\n";
      }
    }
    vecv get_heavy_atom_movable_coords() const { // FIXME mv
      vecv tmp;
      VINA_FOR(i, num_movable_atoms())
        if (!atoms[i].is_hydrogen()) tmp.push_back(coords[i]);
      return tmp;
    }
    void check_internal_pairs() const;
    void print_stuff() const; // FIXME rm
    void print_counts(unsigned nrec_atoms) const;

    fl clash_penalty() const;

    const atomv& get_fixed_atoms() const {
      return grid_atoms;
    }
    const atomv& get_movable_atoms() const {
      return atoms;
    }

    vecv get_extra_box_coords() const {
      //make sure the covalent residue is in the box by providing its coordinates
      return extra_box_coords;
    }

    void clear_minus_forces();
    void add_minus_forces(const std::vector<gfloat3>& forces);
    void sub_minus_forces(const std::vector<gfloat3>& forces);
    void scale_minus_forces(fl scale);
    void round_minus_forces();

    fl get_minus_forces_sum_magnitude() const;

    void set_rigid(const OpenBabel::OBMol& r) {
      rigid = r;
    }

    //allocate gpu memory, model must be setup
    //also copies over data that does not change during minimization
    //if model changes, must re-initialize
    void initialize_gpu();

    bool gpu_initialized() const {
      return gdata.coords != NULL;
    } //true if good to go
    //deallocate gpu memory
    void deallocate_gpu();

    vecv coords;
    vecv extra_box_coords;
    vecv minus_forces;
    gpu_data gdata;
    vector_mutable<ligand> ligands;
    sz m_num_movable_atoms;
    atomv atoms; // movable, inflex
    atomv grid_atoms;
    interacting_pairs other_pairs;

    //for cnn, allow rigid body movement of receptor
    rigid_change rec_change; //set by non_cache/cnn scoring
    rigid_conf rec_conf;

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
    friend class appender;
    friend struct pdbqt_initializer;
    friend struct model_test;
    friend void test_eval_intra();

    const atom& get_atom(const atom_index& i) const {
      return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]);
    }

    atom& get_atom(const atom_index& i) {
      return (i.in_grid ? grid_atoms[i.i] : atoms[i.i]);
    }

    void write_context(const context& c, std::ostream& out) const;
    void write_context(const context& c, std::ostream& out,
        const std::string& remark) const {
      out << remark;
    }
    void write_context(const context& c, const path& name) const {
      ofile out(name);
      write_context(c, out);
    }
    void write_context(const context& c, const path& name,
        const std::string& remark) const {
      ofile out(name);
      write_context(c, out, remark);
    }
    // actually static
    fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const;

    atom_index sz_to_atom_index(sz i) const; // grid_atoms, atoms
    bool bonded_to_HD(const atom& a) const;
    bool bonded_to_heteroatom(const atom& a) const;
    sz find_ligand(sz a) const;
    void bonded_to(sz a, sz n, szv& out) const;
    szv bonded_to(sz a, sz n) const;

    // assign bonds based on relative mobility, distance and covalent length
    void assign_bonds(const distance_type_matrix& mobility);
    void assign_types();
    void initialize_pairs(const distance_type_matrix& mobility);
    void initialize(const distance_type_matrix& mobility);
    fl clash_penalty_aux(const interacting_pairs& pairs) const;

    fl eval_interacting_pairs(const precalculate& p, fl v,
        const interacting_pairs& pairs, const vecv& coords) const;
    fl eval_interacting_pairs_deriv(const precalculate& p, fl v,
        const interacting_pairs& pairs, const vecv& coords, vecv& forces) const;

    bool hydrogens_stripped;
    vecv internal_coords;
    /* TODO:reprivate */
    /* vecv coords; */
    //This contains the accumulated directional deltas for each atom
    /* vecv minus_forces; */
    /*sz m_num_movable_atoms; */
    /*atomv atoms; // movable, inflex*/
    /*atomv grid_atoms;*/

    vector_mutable<residue> flex;
    context flex_context;
    // all except internal to one ligand: ligand-other ligands;
    // ligand-flex/inflex; flex-flex/inflex
    // interacting_pairs other_pairs; 

    OpenBabel::OBMol rigid; //for full_flex_output

    std::string name;
    int pose_num;
};

#endif
