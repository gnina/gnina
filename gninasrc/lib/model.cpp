#include "model.h"
#include "common.h"
#include "file.h"
#include "curl.h"
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include "non_cache_gpu.h"
#include "loop_timer.h"

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

struct appender_info {
    sz grid_atoms_size;
    sz m_num_movable_atoms;
    sz atoms_size;

    appender_info(const model& m)
        : grid_atoms_size(m.grid_atoms.size()),
            m_num_movable_atoms(m.m_num_movable_atoms),
            atoms_size(m.atoms.size()) {
    }
};

class appender {
    appender_info a_info;
    appender_info b_info;
    sz new_grid_index(sz x) const {
      return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
    }
  public:
    bool is_a;

    appender(const model& a, const model& b)
        : a_info(a), b_info(b), is_a(true) {
    }

    sz operator()(sz x) const { // transform coord index
      if (is_a) {
        if (x < a_info.m_num_movable_atoms)
          return x; // a-movable unchanged
        else
          return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
      } else {
        if (x < b_info.m_num_movable_atoms)
          return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
        else
          return x + a_info.atoms_size; // all a's spliced before b-inflex
      }
    }
    atom_index operator()(const atom_index& x) const { // transform atom_index
      atom_index tmp(x);
      if (tmp.in_grid)
        tmp.i = new_grid_index(tmp.i);
      else
        tmp.i = operator()(tmp.i);
      return tmp;
    }

    // type-directed old -> new transformations
    void update(interacting_pair& ip) const {
      ip.a = operator()(ip.a);
      ip.b = operator()(ip.b);
    }
    void update(vec& v) const { // coordinates & forces - do nothing
    }
    void update(ligand& lig) const {
      lig.transform(*this); // ligand as an atom_range subclass
      transform_ranges(lig, *this);
      VINA_FOR_IN(i, lig.pairs)
        this->update(lig.pairs[i]);

      lig.cont.update(*this);
    }
    void update(residue& r) const {
      transform_ranges(r, *this);
    }
    void update(parsed_line& p) const {
      if (p.second) p.second = (*this)(p.second.get());
    }
    void update(atom& a) const {
      VINA_FOR_IN(i, a.bonds) {
        bond& b = a.bonds[i];
        b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
      }
    }

    // ligands, flex, flex_context, atoms; also used for other_pairs
    template<typename T, typename A>
    void append(std::vector<T, A>& a, const std::vector<T, A>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
      sz a_sz = a.size();
      vector_append(a, b);

      is_a = true;
      VINA_FOR(i, a_sz)
        update(a[i]);

      is_a = false;
      VINA_RANGE(i, a_sz, a.size())
        update(a[i]);
    }

    //add b to a
    void append(context& a, const context& b) {
      append(a.pdbqttext, b.pdbqttext);

      if (a.sdftext.valid() && b.sdftext.valid()) {
        //dkoes - I don't think I will need to do this, so not bothering to implement it for now
        abort();
      } else
        if (!a.sdftext.valid()) {
          //but I do need to be able to add flex res to an existing model
          a.sdftext = b.sdftext;
        }
    }

    // internal_coords, coords, minus_forces, atoms
    template<typename T, typename A>
    void coords_append(std::vector<T, A>& a, const std::vector<T, A>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
      std::vector<T, A> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

      is_a = true;
      VINA_FOR_IN(i, a)
        update(a[i]);

      is_a = false;
      VINA_FOR_IN(i, b_copy)
        update(b_copy[i]);

      // interleave
      typedef typename std::vector<T, A>::const_iterator const cci;
      cci b1 = b_copy.begin();
      cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;
      cci b3 = b_copy.end();

      a.insert(a.begin() + a_info.m_num_movable_atoms, b1, b2);
      a.insert(a.end(), b2, b3);
    }
};

//dkoes - modify atom indices according to transform in this context
//I'm not sure this is actually needed/relevant - merging ligands will generally
//just break things
void context::update(const appender& transform) {
  VINA_FOR_IN(i, pdbqttext)
    transform.update(pdbqttext[i]); // parsed_line update

  for (unsigned i = 0, n = sdftext.atoms.size(); i < n; i++) {
    sdftext.atoms[i].index = transform(sdftext.atoms[i].index);
  }
}

//associate pdfbt/sdf data with appropriate atom index from model
void context::set(sz pdbqtindex, sz sdfindex, sz atomindex, bool inf) {
  if (pdbqtindex < pdbqttext.size()) {
    if (!inf) //inf is true if a nonmoving atom, which will just use the old coordinates
      pdbqttext[pdbqtindex].second = atomindex;

  }
  if (sdfindex < sdftext.atoms.size()) {
    sdftext.atoms[sdfindex].index = atomindex;
    sdftext.atoms[sdfindex].inflex = inf;
  }
}

void model::append(const model& m) {
  deallocate_gpu();
  appender t(*this, m);

  hydrogens_stripped |= m.hydrogens_stripped;
  t.append(other_pairs, m.other_pairs);

  VINA_FOR_IN(i, atoms)
    VINA_FOR_IN(j, m.atoms) {
      if (i >= m_num_movable_atoms && j >= m.m_num_movable_atoms) continue; // no need for inflex-inflex interactions

      const atom& a = atoms[i];
      const atom& b = m.atoms[j];

      smt t1 = a.get();
      smt t2 = b.get();
      sz n = num_atom_types();

      if (t1 < n && t2 < n) {
        t.is_a = true;
        sz new_i = t(i);
        t.is_a = false;
        sz new_j = t(j);
        other_pairs.push_back(interacting_pair(t1, t2, new_i, new_j));
      }
    }

  VINA_CHECK(minus_forces.size() == coords.size());
  VINA_CHECK(m.minus_forces.size() == m.coords.size());

  t.coords_append(internal_coords, m.internal_coords);
  t.coords_append(coords, m.coords);
  t.coords_append(minus_forces, m.minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

  t.append(ligands, m.ligands);
  t.append(flex, m.flex);
  t.append(flex_context, m.flex_context);

  t.append(grid_atoms, m.grid_atoms);
  t.coords_append(atoms, m.atoms);

  m_num_movable_atoms += m.m_num_movable_atoms;

//initialize_gpu();

}

///////////////////  end  MODEL::APPEND /////////////////////////

//helper functions for strip_hydrogesn

//remove any pairs with atoms that aren't in atommap
static void striph_interacting_pairs(const std::vector<int>& atommap,
    interacting_pairs& ip) {
  interacting_pairs newpairs;
  newpairs.reserve(ip.size());
  VINA_FOR_IN(i, ip) {
    interacting_pair p = ip[i];
    assert(p.a < atommap.size());
    assert(p.b < atommap.size());
    if (atommap[p.a] >= 0 && atommap[p.b] >= 0) {
      p.a = atommap[p.a];
      p.b = atommap[p.b];
      newpairs.push_back(p);
    }
  }
  ip.swap(newpairs);
}

template<typename T>
static void striph_coord_vec(const std::vector<int>& atommap, unsigned natm,
    std::vector<T>& vec) {
  //remove elements of vec at indices that aren't remmaped in atommap
  std::vector<T> newvec;
  newvec.resize(natm);
  assert(vec.size() == atommap.size());
  VINA_FOR_IN(i, vec) {
    int newi = atommap[i];
    if (newi >= 0) {
      //these should be one after another..
      newvec[newi] = vec[i];
    }
  }
  vec.swap(newvec);
}

//recursively reduce atom_ranges
template<typename T>
static void striph_ranges(const std::vector<int>& atommap, T& t) {
  t.node.reduce(atommap);
  VINA_FOR_IN(i, t.children) {
    striph_ranges(atommap, t.children[i]);
  }
}

static void striph_context(const std::vector<int>& atommap, context& c) {

  //first pdbqt, change atom indices
  pdbqtcontext newpc;
  newpc.reserve(c.pdbqttext.size());
  VINA_FOR_IN(i, c.pdbqttext) {
    parsed_line p = c.pdbqttext[i];
    if (p.second) {
      int newi = atommap[p.second.get()];
      if (newi >= 0) {
        p.second = newi;
        newpc.push_back(p);
      }
    } else {
      newpc.push_back(p);
    }
  }
  c.pdbqttext.swap(newpc);

  //sdf if a little trickier..
  if (c.sdftext.atoms.size() > 0) {

    std::vector<sdfcontext::sdfatom> newatoms;
    newatoms.reserve(c.sdftext.atoms.size());
    std::vector<sdfcontext::sdfbond> newbonds;
    newbonds.reserve(c.sdftext.bonds.size());
    std::vector<sdfcontext::sdfprop> newprops;
    newprops.reserve(c.sdftext.properties.size());

    //sdf atoms have their own mapping..
    std::vector<int> sdfmap(c.sdftext.atoms.size(), -1); //sdfbond indices are into the sdf atoms array, keep updated map

    VINA_FOR_IN(i, c.sdftext.atoms) { //note these atoms are out of order
      sdfcontext::sdfatom a = c.sdftext.atoms[i];
      int newi = atommap[a.index];
      if (newi >= 0) {
        a.index = newi;
        sdfmap[i] = newatoms.size();
        newatoms.push_back(a);
      }
    }

    VINA_FOR_IN(i, c.sdftext.bonds) {
      sdfcontext::sdfbond bond = c.sdftext.bonds[i];
      int newa = sdfmap[bond.a];
      int newb = sdfmap[bond.b];
      if (newa >= 0 && newb >= 0) {
        bond.a = newa;
        bond.b = newb;
        newbonds.push_back(bond);
      }
    }

    VINA_FOR_IN(i, c.sdftext.properties) {
      sdfcontext::sdfprop p = c.sdftext.properties[i];
      int newi = sdfmap[p.atom];
      if (newi >= 0) {
        p.atom = newi;
        newprops.push_back(p);
      }
    }

    c.sdftext.atoms.swap(newatoms);
    c.sdftext.bonds.swap(newbonds);
    c.sdftext.properties.swap(newprops);
  }
}

static void striph_ligand(const std::vector<int>& atommap, ligand& lig) {
  lig.reduce(atommap); // ligand as an atom_range subclass
  striph_ranges(atommap, lig);
  striph_interacting_pairs(atommap, lig.pairs);
  striph_context(atommap, lig.cont);
}

static void striph_update_bonds(const std::vector<int>& gridatommap,
    const std::vector<int>& atommap, atom& a) {

  std::vector<bond> newbonds;
  newbonds.reserve(a.bonds.size());
  VINA_FOR_IN(i, a.bonds) {
    bond b = a.bonds[i];
    atom_index& ai = b.connected_atom_index;
    int newi = -1;
    if (ai.in_grid)
      newi = gridatommap[ai.i];
    else
      newi = atommap[ai.i];
    if (newi >= 0) {
      ai.i = newi;
      newbonds.push_back(b);
    }
  }
  a.bonds.swap(newbonds);
}

//Remove hydrogens from model in-place.  Must be called after final assignment of atom types.
void model::strip_hydrogens() {
  deallocate_gpu();
  hydrogens_stripped = true;
  sz N = num_atom_types();

  //update grid atoms
  std::vector<int> gridatommap(grid_atoms.size(), -1);
  atomv newgridatoms;
  newgridatoms.reserve(grid_atoms.size());
  VINA_FOR_IN(i, grid_atoms) {
    smt typ = grid_atoms[i].get();
    if (typ < N && !is_hydrogen(typ)) {
      gridatommap[i] = newgridatoms.size();
      newgridatoms.push_back(grid_atoms[i]);
    }
  }

  //atoms vector - several other arrays use the same indexing so maintain a
  std::vector<int> atommap(atoms.size(), -1); //mapping from atoms array index to new indexing
  atomv newatoms;
  newatoms.reserve(atoms.size());
  sz n_good_atoms = 0;
  sz n_good_moveable = 0;

  VINA_FOR_IN(i, atoms) {
    smt typ = atoms[i].get();
    if (typ < N && !is_hydrogen(typ)) {
      atommap[i] = newatoms.size();
      newatoms.push_back(atoms[i]);
      n_good_atoms++;
      if (i < m_num_movable_atoms) {
        //flex has some non-moveable atoms still in atoms
        n_good_moveable++;
      }
    }
  }

  //now update bonds
  VINA_FOR_IN(i, newgridatoms) {
    striph_update_bonds(gridatommap, atommap, newgridatoms[i]);
  }

  VINA_FOR_IN(i, newatoms) {
    striph_update_bonds(gridatommap, atommap, newatoms[i]);
  }

  //set atoms arrays
  grid_atoms.swap(newgridatoms);
  atoms.swap(newatoms);

  m_num_movable_atoms = n_good_moveable;

  //both atoms in other pairs need to okay
  striph_interacting_pairs(atommap, other_pairs);

  striph_coord_vec(atommap, n_good_atoms, internal_coords);
  striph_coord_vec(atommap, n_good_atoms, coords);
  striph_coord_vec(atommap, n_good_atoms, minus_forces); //really could just resize..

  VINA_FOR_IN(i, ligands) {
    ligand& lig = ligands[i];
    striph_ligand(atommap, lig);
  }

  VINA_FOR_IN(r, flex) {
    residue& res = flex[r];
    striph_ranges(atommap, res);
  }

  striph_context(atommap, flex_context);

}

struct branch_metrics {
    sz length;
    sz corner2corner;
    branch_metrics()
        : length(0), corner2corner(0) {
    }
};

template<typename T>
branch_metrics get_branch_metrics(const T& t) {
  branch_metrics tmp;
  if (!t.children.empty()) {
    sz corner2corner_max = 0;
    szv lengths;
    VINA_FOR_IN(i, t.children) {
      branch_metrics res = get_branch_metrics(t.children[i]);
      if (corner2corner_max < res.corner2corner)
        corner2corner_max = res.corner2corner;
      lengths.push_back(res.length + 1); // FIXME? weird compiler warning (sz -> unsigned)
    }
    std::sort(lengths.begin(), lengths.end());

    tmp.length = lengths.back();

    tmp.corner2corner = tmp.length;
    if (lengths.size() >= 2) tmp.corner2corner += lengths[lengths.size() - 1];

    if (tmp.corner2corner < corner2corner_max) tmp.corner2corner =
        corner2corner_max;
  }
  return tmp;
}

sz model::ligand_longest_branch(sz ligand_number) const {
  return get_branch_metrics(ligands[ligand_number]).length;
}

sz model::ligand_length(sz ligand_number) const {
  return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

template<typename T>
atom_range get_atom_range(const T& t) {
  atom_range tmp = t.node;
  VINA_FOR_IN(i, t.children) {
    atom_range r = get_atom_range(t.children[i]);
    if (tmp.begin > r.begin) tmp.begin = r.begin;
    if (tmp.end < r.end) tmp.end = r.end;
  }
  return tmp;
}

void ligand::set_range() {
  atom_range tmp = get_atom_range(*this);
  begin = tmp.begin;
  end = tmp.end;
}

/////////////////// begin MODEL::INITIALIZE /////////////////////////

atom_index model::sz_to_atom_index(sz i) const {
  if (i < grid_atoms.size())
    return atom_index(i, true);
  else
    return atom_index(i - grid_atoms.size(), false);
}

distance_type model::distance_type_between(const distance_type_matrix& mobility,
    const atom_index& i, const atom_index& j) const {
  if (i.in_grid && j.in_grid) return DISTANCE_FIXED;
  if (i.in_grid)
    return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
  if (j.in_grid)
    return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
  assert(!i.in_grid);
  assert(!j.in_grid);
  assert(i.i < atoms.size());
  assert(j.i < atoms.size());
  sz a = i.i;
  sz b = j.i;
  if (a == b) return DISTANCE_FIXED;
  return (a < b) ? mobility(a, b) : mobility(b, a);
}

const vec& model::atom_coords(const atom_index& i) const {
  return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}

fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
  return vec_distance_sqr(atom_coords(a), atom_coords(b));
}

struct bond_less { // FIXME rm!?
    bool operator()(const bond& a, const bond& b) const {
      return a.connected_atom_index.i < b.connected_atom_index.i;
    }
};

bool model::atom_exists_between(const distance_type_matrix& mobility,
    const atom_index& a, const atom_index& b, const szv& relevant_atoms) const { // there is an atom closer to both a and b then they are to each other and immobile relative to them
  fl r2 = distance_sqr_between(a, b);
  VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
    sz i = relevant_atoms[relevant_atoms_i];
    atom_index c = sz_to_atom_index(i);
    if (a == c || b == c) continue;
    if (is_hydrogen(get_atom(c).sm)) continue; //don't let hydrogens lock things down
    distance_type ac = distance_type_between(mobility, a, c);
    distance_type bc = distance_type_between(mobility, b, c);
    if (ac != DISTANCE_VARIABLE && bc != DISTANCE_VARIABLE
        && distance_sqr_between(a, c) < r2 && distance_sqr_between(b, c) < r2) {
      return true;
    }
  }
  return false;
}

struct beads {
    fl radius_sqr;
    std::vector<std::pair<vec, szv> > data;
    beads(sz reserve_size, fl radius_sqr_)
        : radius_sqr(radius_sqr_) {
      data.reserve(reserve_size);
    }
    void add(sz index, const vec& coords) {
      VINA_FOR_IN(i, data) {
        if (vec_distance_sqr(coords, data[i].first) < radius_sqr) {
          data[i].second.push_back(index);
          return;
        }
      }
      // not found
      std::pair<vec, szv> tmp;
      tmp.first = coords;
      tmp.second.push_back(index);
      data.push_back(tmp);
    }
};

void model::assign_bonds(const distance_type_matrix& mobility) { // assign bonds based on relative mobility, distance and covalent length
  const fl bond_length_allowance_factor = 1.1;
  sz n = grid_atoms.size() + atoms.size();

  // construct beads
  const fl bead_radius = 15;
  beads beads_instance(n, sqr(bead_radius));
  VINA_FOR(i, n) {
    atom_index i_atom_index = sz_to_atom_index(i);
    beads_instance.add(i, atom_coords(i_atom_index));
  }
// assign bonds
  VINA_FOR(i, n) {
    atom_index i_atom_index = sz_to_atom_index(i);
    const vec& i_atom_coords = atom_coords(i_atom_index);
    atom& i_atom = get_atom(i_atom_index);

    const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
    fl i_atom_covalent_radius = covalent_radius(i_atom.sm);

    //find relevant atoms
    szv relevant_atoms;
    const fl bead_cutoff_sqr = sqr(
        bead_radius
            + bond_length_allowance_factor
                * (i_atom_covalent_radius + max_covalent_r));
    VINA_FOR_IN(b, beads_instance.data) {
      if (vec_distance_sqr(beads_instance.data[b].first, i_atom_coords)
          > bead_cutoff_sqr) continue;
      const szv& bead_elements = beads_instance.data[b].second;
      VINA_FOR_IN(bead_elements_i, bead_elements) {
        sz j = bead_elements[bead_elements_i];
        atom_index j_atom_index = sz_to_atom_index(j);
        atom& j_atom = get_atom(j_atom_index);
        const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
        distance_type dt = distance_type_between(mobility, i_atom_index,
            j_atom_index);
        if (dt != DISTANCE_VARIABLE && i != j) {
          fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
          //if(r2 < sqr(bond_length_allowance_factor * bond_length))
          if (r2
              < sqr(
                  bond_length_allowance_factor
                      * (i_atom_covalent_radius + max_covalent_r)))
            relevant_atoms.push_back(j);
        }
      }
    }
    // find bonded atoms
    VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
      sz j = relevant_atoms[relevant_atoms_i];
      if (j <= i) continue; // already considered
      atom_index j_atom_index = sz_to_atom_index(j);
      atom& j_atom = get_atom(j_atom_index);
      const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
      distance_type dt = distance_type_between(mobility, i_atom_index,
          j_atom_index);
      fl r2 = distance_sqr_between(i_atom_index, j_atom_index);

      if (r2 < sqr(bond_length_allowance_factor * bond_length)
          && !atom_exists_between(mobility, i_atom_index, j_atom_index,
              relevant_atoms)) {
        bool rotatable = (dt == DISTANCE_ROTOR);
        fl length = std::sqrt(r2);
        i_atom.bonds.push_back(bond(j_atom_index, length, rotatable));
        j_atom.bonds.push_back(bond(i_atom_index, length, rotatable));
      }

    }
  }
}

bool model::bonded_to_HD(const atom& a) const {
  VINA_FOR_IN(i, a.bonds) {
    const bond& b = a.bonds[i];
    if (get_atom(b.connected_atom_index).sm == smina_atom_type::PolarHydrogen)
      return true;
  }
  return false;
}

bool model::bonded_to_heteroatom(const atom& a) const {
  VINA_FOR_IN(i, a.bonds) {
    const bond& b = a.bonds[i];
    if (get_atom(b.connected_atom_index).is_heteroatom()) return true;
  }
  return false;
}

//dkoes - modify smina types as necessary to take into account bonding information
void model::assign_types() {
  VINA_CHECK(!hydrogens_stripped); //type assignment requires hydrogens

  VINA_FOR(i, grid_atoms.size() + atoms.size()) {
    const atom_index ai = sz_to_atom_index(i);
    atom& a = get_atom(ai);

    //this is where the X-scale atom types diverge from autodock
    a.sm = adjust_smina_type(a.sm, bonded_to_HD(a), bonded_to_heteroatom(a));
  }
}

sz model::find_ligand(sz a) const {
  VINA_FOR_IN(i, ligands)
    if (a >= ligands[i].begin && a < ligands[i].end) return i;
  return ligands.size();
}

void model::bonded_to(sz a, sz n, szv& out) const {
  if (!has(out, a)) { // not found
    out.push_back(a);
    if (n > 0)
      VINA_FOR_IN(i, atoms[a].bonds) {
        const bond& b = atoms[a].bonds[i];
        if (!b.connected_atom_index.in_grid)
          bonded_to(b.connected_atom_index.i, n - 1, out);
      }
  }
}

szv model::bonded_to(sz a, sz n) const {
  szv tmp;
  bonded_to(a, n, tmp);
  return tmp;
}

void model::initialize_pairs(const distance_type_matrix& mobility) {
  VINA_FOR_IN(i, atoms) {
    sz i_lig = find_ligand(i);
    szv bonded_atoms = bonded_to(i, 3);
    sz n = num_atom_types();
    VINA_RANGE(j, i+1, atoms.size()) {
      if (i >= m_num_movable_atoms && j >= m_num_movable_atoms) continue; // exclude inflex-inflex
      if (mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
        smt t1 = atoms[i].get();
        smt t2 = atoms[j].get();
        if (t1 < n && t2 < n && !is_hydrogen(t1) && !is_hydrogen(t2)) { //exclude, say, Hydrogens
          interacting_pair ip(t1, t2, i, j);
          if (i_lig < ligands.size() && find_ligand(j) == i_lig)
            ligands[i_lig].pairs.push_back(ip);
          else
            other_pairs.push_back(ip);
        }
      }
    }
  }
}

void model::initialize(const distance_type_matrix& mobility) {
  VINA_FOR_IN(i, ligands)
    ligands[i].set_range();
  assign_bonds(mobility);
  assign_types();
  initialize_pairs(mobility);
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////

sz model::num_internal_pairs() const {
  sz tmp = 0;
  VINA_FOR_IN(i, ligands)
    tmp += ligands[i].pairs.size();
  return tmp;
}

void model::get_movable_atom_types(std::vector<smt>& movingtypes) const {
  sz n = num_atom_types();
  movingtypes.clear();
  movingtypes.reserve(n);
  VINA_FOR(i, m_num_movable_atoms) {
    const atom& a = atoms[i];
    smt t = a.get();
    if (t < n && !has(movingtypes, t)) movingtypes.push_back(t);
  }
}

conf_size model::get_size() const {
  conf_size tmp;
  tmp.ligands = ligands.count_torsions();
  tmp.flex = flex.count_torsions();
  return tmp;
}

conf model::get_initial_conf(bool enable_receptor) const { // torsions = 0, orientations = identity, ligand positions = current
  conf_size cs = get_size();
  conf tmp(cs, enable_receptor);
  tmp.set_to_null();
  VINA_FOR_IN(i, ligands)
    tmp.ligands[i].rigid.position = ligands[i].node.get_origin();
  return tmp;
}

grid_dims model::movable_atoms_box(fl add_to_each_dimension,
    fl granularity) const {
  vec corner1(0, 0, 0), corner2(0, 0, 0);
  VINA_FOR(i, num_movable_atoms()) {
    const vec& v = movable_coords(i);
    VINA_FOR_IN(j, v) {
      if (i == 0 || v[j] < corner1[j]) corner1[j] = v[j];
      if (i == 0 || v[j] > corner2[j]) corner2[j] = v[j];
    }
  }
  corner1 -= add_to_each_dimension;
  corner2 += add_to_each_dimension;

  grid_dims gd;
  { // always doing this now FIXME ?
    vec center;
    center = 0.5 * (corner2 + corner1);
    VINA_FOR_IN(i, gd) {
      gd[i].n = sz(std::ceil((corner2[i] - corner1[i]) / granularity));
      fl real_span = granularity * gd[i].n;
      gd[i].begin = center[i] - real_span / 2;
      gd[i].end = gd[i].begin + real_span;
    }
  }
  return gd;
}

void string_write_coord(sz i, fl x, std::string& str) {
  VINA_CHECK(i > 0);
  --i;
  std::ostringstream out;
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << std::setw(8) << std::setprecision(3) << x;
  VINA_CHECK(out.str().size() == 8);
  VINA_CHECK(str.size() > i + 8);
  VINA_FOR(j, 8)
    str[i + j] = out.str()[j];
}
std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) {
  std::string tmp(str);
  string_write_coord(31, coords[0], tmp);
  string_write_coord(39, coords[1], tmp);
  string_write_coord(47, coords[2], tmp);
  return tmp;
}

void context::writePDBQT(const vecv& coords, std::ostream& out) const {
  VINA_FOR_IN(i, pdbqttext) {
    const std::string& str = pdbqttext[i].first;
    if (pdbqttext[i].second) {
      out << coords_to_pdbqt_string(coords[pdbqttext[i].second.get()], str)
          << '\n';
    } else
      if (boost::starts_with(str, "BEGIN_RES")
          || boost::starts_with(str, "END_RES")) {
        //dkoes - openbabel thinks these denote separate molecules
        //so we lose all but the first residue if we leave them in
      } else
        out << str << '\n';
  }
}

void sdfcontext::dump(std::ostream& out) const {
  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    const sdfatom& atom = atoms[i];
    out << atom.index << " " << atom.elem << "\n";
  }
}

//output sdf format to out
void sdfcontext::write(const vecv& coords, sz nummove,
    std::ostream& out) const {
  const unsigned bsize = 1024;
  char buff[bsize]; //since sprintf is just so much easier to use
  //name followed by two blank lines
  out << name << "\n\n\n";

  //cnts and version line
  snprintf(buff, bsize, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
      (int) atoms.size(), (int) bonds.size());
  out << buff;

  //atom block
  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    const sdfatom& atom = atoms[i];
    sz idx = atom.index;
    if (atom.inflex) idx += nummove; //rigids are after movable
    const vec& c = coords[idx];
    assert(idx < coords.size());
    snprintf(buff, bsize,
        "%10.4f%10.4f%10.4f %-3.2s 0  0  0  0  0  0  0  0  0  0  0  0\n", c[0],
        c[1], c[2], atom.elem);
    out << buff;
  }

  //bond block
  for (unsigned i = 0, n = bonds.size(); i < n; i++) {
    const sdfbond& bond = bonds[i];
    out << std::setw(3) << bond.a + 1; //indexed from one
    out << std::setw(3) << bond.b + 1;
    out << std::setw(3) << (int) bond.type;
    out << "  0  0  0\n";
  }

  //properties
  for (unsigned i = 0, n = properties.size(); i < n; i++) {
    const sdfprop& prop = properties[i];
    if (prop.type == 'c') //M CHG
        {
      out << "M  CHG 1 " << std::setw(3) << prop.atom + 1 << std::setw(4)
          << (int) prop.value << "\n";
    } else
      if (prop.type == 'i') //M  ISO
          {
        out << "M  ISO 1 " << std::setw(3) << prop.atom + 1 << std::setw(4)
            << (int) prop.value << "\n";
      }
  }

  //end, but leave room for sddata
  out << "M  END\n";
}

void model::write_context(const context& c, std::ostream& out) const {
  verify_bond_lengths();
  c.writePDBQT(coords, out);
}

//more for debugging, dump fixed atoms as xyz to out
//applies transformation
void model::write_rigid_xyz(std::ostream& out, const vec& center) const {
  out << grid_atoms.size() << "\n\n";
  gfloat3 c = gfloat3(center[0], center[1], center[2]);
  gfloat3 t = gfloat3(rec_conf.position[0], rec_conf.position[1],
      rec_conf.position[2]);
  VINA_FOR_IN(i, grid_atoms) {
    const atom& a = grid_atoms[i];
    out << smina_type_to_element_name(a.sm) << " ";
    gfloat3 pt = rec_conf.orientation.transform(a.coords[0], a.coords[1],
        a.coords[2], c, t);
    out << pt.x << " " << pt.y << " " << pt.z << "\n";
  }
}

void model::seti(const conf& c) {
  /* TODO */
  assert(0);
  /* ligands.set_conf(atoms, internal_coords, c.ligands); */
}

void model::sete(const conf& c) {
  VINA_FOR_IN(i, ligands)
    c.ligands[i].rigid.apply(internal_coords, coords, ligands[i].begin,
        ligands[i].end);
  /* TODO */
  flex.set_conf(atoms, coords, c.flex);
}

void model::set(const conf& c) {
  ligands.set_conf(atoms, coords, c.ligands);
  flex.set_conf(atoms, coords, c.flex);
  //for cnn, we do not change the receptor coordinates here
  //instead the cnn layer applies the rigid body transformation, which will
  //apply the inverse of to the ligand when we are done
  rec_conf = c.receptor;
}

//dkoes - return the string corresponding to i'th ligand atoms pdb information
//which is serial+name
std::string model::ligand_atom_str(sz i, sz lig) const {
  if (atoms.size() != coords.size()) abort();
  assert(i < atoms.size());
  std::string pdbline;
  const context& cont = ligands[lig].cont;
  for (sz c = 0, nc = cont.pdbqtsize(); c < nc; c++) {
    /* TODO: nvcc flags error here:
     ../../../src/lib/model.cpp:862: warning: integer conversion
     resulted in a change of sign
     */
    if (cont.pdbqttext[c].second.get_value_or(-1) == i) {
      pdbline = ligands[lig].cont.pdbqttext[c].first;
      break;
    }
  }

  if (pdbline.size() > 6)
    return pdbline.substr(6, 9);
  else
    return "";
}

fl model::gyration_radius(sz ligand_number) const {
  VINA_CHECK(ligand_number < ligands.size());
  const ligand& lig = ligands[ligand_number];
  fl acc = 0;
  unsigned counter = 0;
  VINA_RANGE(i, lig.begin, lig.end) {
    if (!atoms[i].is_hydrogen()) { // only heavy atoms are used
      acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
      ++counter;
    }
  }
  return (counter > 0) ? std::sqrt(acc / counter) : 0;
}

fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
  sz n = x.m_num_movable_atoms;
  VINA_CHECK(n == y.m_num_movable_atoms);
  fl sum = 0;
  unsigned counter = 0;
  VINA_FOR(i, n) {
    const atom& a = x.atoms[i];
    if (!a.is_hydrogen()) {
      fl r2 = max_fl;
      VINA_FOR(j, n) {
        const atom& b = y.atoms[j];
        if (a.same_element(b) && !b.is_hydrogen()) {
          fl this_r2 = vec_distance_sqr(x.coords[i], y.coords[j]);
          if (this_r2 < r2) r2 = this_r2;
        }
      }
      assert(not_max(r2));
      sum += r2;
      ++counter;
    }
  }
  return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_lower_bound(const model& m) const {
  return (std::max)(rmsd_lower_bound_asymmetric(*this, m),
      rmsd_lower_bound_asymmetric(m, *this));
}

fl model::rmsd_upper_bound(const model& m) const {
  VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
  fl sum = 0;
  unsigned counter = 0;
  VINA_FOR(i, m_num_movable_atoms) {
    const atom& a = atoms[i];
    const atom& b = m.atoms[i];
    assert(a.sm == b.sm);
    if (!a.is_hydrogen()) {
      sum += vec_distance_sqr(coords[i], m.coords[i]);
      ++counter;
    }
  }
  return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

fl model::rmsd_ligands_upper_bound(const model& m) const {
  VINA_CHECK(ligands.size() == m.ligands.size());
  fl sum = 0;
  unsigned counter = 0;
  VINA_FOR_IN(ligand_i, ligands) {
    const ligand& lig = ligands[ligand_i];
    const ligand& m_lig = m.ligands[ligand_i];
    VINA_CHECK(lig.begin == m_lig.begin);
    VINA_CHECK(lig.end == m_lig.end);
    VINA_RANGE(i, lig.begin, lig.end) {
      const atom& a = atoms[i];
      const atom& b = m.atoms[i];
      assert(a.sm == b.sm);
      if (!a.is_hydrogen()) {
        sum += vec_distance_sqr(coords[i], m.coords[i]);
        ++counter;
      }
    }
  }
  return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

void model::verify_bond_lengths() const {
  VINA_FOR(i, grid_atoms.size() + atoms.size()) {
    const atom_index ai = sz_to_atom_index(i);
    const atom& a = get_atom(ai);
    VINA_FOR_IN(j, a.bonds) {
      const bond& b = a.bonds[j];
      fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
      bool ok = eq(d, b.length);
      if (!ok) {
        VINA_SHOW(d);
        VINA_SHOW(b.length);
      }
      VINA_CHECK(ok);
    }
  }
}

void model::check_internal_pairs() const {
  VINA_FOR_IN(i, ligands) {
    const ligand& lig = ligands[i];
    const interacting_pairs& pairs = lig.pairs;
    VINA_FOR_IN(j, pairs) {
      const interacting_pair& ip = pairs[j];
      VINA_CHECK(ip.a >= lig.begin);
      VINA_CHECK(ip.b < lig.end);
    }
  }
}

void model::about() const {
  VINA_SHOW(num_movable_atoms());
  VINA_SHOW(num_internal_pairs());
  VINA_SHOW(num_other_pairs());
  VINA_SHOW(num_ligands());
  VINA_SHOW(num_flex());
}

void model::print_stuff() const {
  std::cout << "coords:\n";
  VINA_FOR_IN(i, coords)
    printnl(coords[i]);

  std::cout << "internal_coords:\n";
  VINA_FOR_IN(i, internal_coords)
    printnl(internal_coords[i]);

  std::cout << "atoms:\n";
  VINA_FOR_IN(i, atoms) {
    const atom& a = atoms[i];
    std::cout << a.sm << "    " << a.charge << '\n';
    std::cout << a.bonds.size() << "  ";
    printnl(a.coords);
  }

  std::cout << "grid_atoms:\n";
  VINA_FOR_IN(i, grid_atoms) {
    const atom& a = grid_atoms[i];
    std::cout << a.sm << "    " << a.charge << '\n';
    std::cout << a.bonds.size() << "  ";
    printnl(a.coords);
  }
  about();
}

void model::print_counts(unsigned nrec_atoms) const {
  std::cout << "ligands :" << ligands.size() << "\n" << "torsions: "
      << ligands.count_torsions()[0] << "\n" << "lig atoms: " << coords.size()
      << "\n" << "rec atoms(pruned): " << nrec_atoms << "\n";
}

fl pairwise_clash_penalty(fl r, fl covalent_r) {
  // r = 0          -> max_penalty
  // r = covalent_r -> 1
  // elsewhere      -> hyperbolic function
  assert(r >= 0);
  assert(covalent_r > epsilon_fl);
  const fl x = r / covalent_r;
  if (x > 2) return 0;
  return 1 - x * x / 4;
}

fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
  fl e = 0;
  VINA_FOR_IN(i, pairs) {
    const interacting_pair& ip = pairs[i];
    const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
    const fl covalent_r = atoms[ip.a].covalent_radius()
        + atoms[ip.b].covalent_radius();
    e += pairwise_clash_penalty(r, covalent_r);
  }
  return e;
}

fl model::clash_penalty() const {
  fl e = 0;
  VINA_FOR_IN(i, ligands)
    e += clash_penalty_aux(ligands[i].pairs);
  e += clash_penalty_aux(other_pairs);
  return e;
}
