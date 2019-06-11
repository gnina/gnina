/*
 * PDBQTUtilities.cpp
 *
 *  Created on: Jun 4, 2014
 *      Author: dkoes
 *
 *  This is a bunch of stuff copied from pdqtformat.cpp from OpenBabel.
 */

/**********************************************************************
 Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
 Some portions Copyright (C) 2003-2006 Geoffrey R. Hutchison
 Some portions Copyright (C) 2004 by Chris Morley
 Some portions Copyright (C) 2014 by David Koes and the University of Pittsburgh

 Original Copyright refers to the pdbformat.cpp file, for reading and
 writing pdb format files.
 Extensively modified 2010 Stuart Armstrong (Source Science/InhibOx)
 for the purpose of reading and writing pdbqt format files.
 Some portions Copyright (C) 2010 by Stuart Armstrong of Source Science/
 InhibOx

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 ***********************************************************************/

#include "PDBQTUtilities.h"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>

using namespace std;
using namespace OpenBabel;

unsigned int FindFragments(OpenBabel::OBMol mol,
    std::vector<std::vector<int> >& rigid_fragments) {
  std::vector<int> nr;
  return FindFragments(mol, rigid_fragments, 0, nr);
}

//return the best root atom and fill in the rigid fragments
//if desired_root is set and rotateAroundRoot is true, let bonds rotate about it
unsigned int FindFragments(OBMol mol, vector<vector<int> >& rigid_fragments,
    unsigned desired_root, const vector<int>& norotate) {
  unsigned int best_root_atom = 1;
  unsigned int shortest_maximal_remaining_subgraph = mol.NumAtoms();
  for (unsigned int i = 1; i <= mol.NumAtoms(); i++)
  //finds the root atom by copying the molecule, deleting each atom in turn, and finding the sizes of the resulting pieces
      {
    OBMol mol_pieces = mol;
    OBAtom * atom_to_del = mol_pieces.GetAtom(i);
    vector<vector<int> > frag_list;

    mol_pieces.DeleteAtom(atom_to_del, true);
    mol_pieces.ContigFragList(frag_list);
    unsigned int smrsi = 0;
    for (unsigned int j = 0; j < frag_list.size(); j++) {
      smrsi = smrsi > frag_list.at(j).size() ? smrsi : frag_list.at(j).size();
    }
    if (smrsi < shortest_maximal_remaining_subgraph) {
      shortest_maximal_remaining_subgraph = smrsi;
      best_root_atom = i;
    }
  }

  boost::unordered_set<int> norot(norotate.begin(), norotate.end());
  vector<unsigned int> bonds_to_delete;
  {
    OBMol mol_pieces = mol;
    for (OBBondIterator it = mol_pieces.BeginBonds();
        it != mol_pieces.EndBonds(); it++) {
      OBBond *bond = *it;
      int src = bond->GetBeginAtomIdx();
      int dst = bond->GetEndAtomIdx();

      if (norot.count(src) || norot.count(dst)) {
        //not rotatable
      } else
        if (IsRotBond_PDBQT(bond, desired_root)) {
          bonds_to_delete.push_back((*it)->GetIdx());
        }
    }

    if (bonds_to_delete.size() != 0) //checks there is something to delete
        {
      vector<unsigned int>::iterator itb = bonds_to_delete.end();
      itb--;
      for (OBBondIterator it = mol_pieces.EndBonds(); true;) {
        it--;
        if ((*it)->GetIdx() == (*itb)) {
          mol_pieces.DeleteBond((*it), true);
          if (itb == bonds_to_delete.begin()) {
            break;
          } else {
            itb--;
          }
        }
      }
    }
    mol_pieces.ContigFragList(rigid_fragments);
  }
  return best_root_atom;
}

bool IsRotBond_PDBQT(OBBond * the_bond, unsigned desired_root)
//identifies a bond as rotatable if it is a single bond, not amide, not in a ring,
//and if both atoms it connects have at least one other atom bounded to them
//will also allow bonds to the desired root to rotate
    {
  if (the_bond->GetBondOrder() != 1 || the_bond->IsAmide()
      || the_bond->IsInRing()) {
    return false;
  }
  if ((GET_HVY(the_bond->GetBeginAtom()) == 1)
      || (GET_HVY(the_bond->GetEndAtom()) == 1)) {
    if (the_bond->GetBeginAtomIdx() == desired_root)
      return true;
    else
      if (the_bond->GetEndAtomIdx() == desired_root) return true;

    return false;
  }
  return true;
}

bool IsIn(const vector<int>& vec, const int num) //checks whether a vector of int contains a specific int
    {
  for (vector<int>::const_iterator itv = vec.begin(); itv != vec.end(); itv++) {
    if ((*itv) == num) {
      return true;
    }
  }
  return false;
}

void ConstructTree(map<unsigned int, obbranch>& tree,
    vector<vector<int> > rigid_fragments, unsigned int root_piece,
    const OBMol& mol, bool flexible) {
  unsigned int first_atom = 0;
  unsigned int second_atom = 0;
  unsigned int first_atom_rank = 0;
  unsigned int second_atom_rank = 0;

  obbranch sprog;

  sprog.atoms = rigid_fragments.at(root_piece);
  sprog.rigid_with.insert(0);

  tree.insert(pair<unsigned int, obbranch>(0, sprog));

  rigid_fragments.erase(rigid_fragments.begin() + root_piece);

  unsigned int position = 0;
  unsigned int atoms_moved = 0;
  bool fecund;
  while (!tree[0].done) {
    fecund = !tree[position].done;
    if (fecund) {
      bool sterile = true;
      for (unsigned int i = 0; i < rigid_fragments.size(); i++) {
        if (FindBondedPiece((*tree.find(position)).second.atoms,
            rigid_fragments.at(i), first_atom, second_atom, first_atom_rank,
            second_atom_rank, mol, atoms_moved)) {
          sprog.connecting_atom_parent = first_atom;
          sprog.connecting_atom_branch = second_atom;
          sprog.how_many_atoms_moved = atoms_moved;
          sprog.atoms = rigid_fragments.at(i);

          sprog.depth = (*tree.find(position)).second.depth + 1;
          sprog.parents = (*tree.find(position)).second.parents; //all parents of the parent are parents too
          sprog.parents.push_back(tree.size()); //a branch is its own parent
          sprog.index = tree.size(); //the index is simply the number of precursors
          sprog.rigid_with.clear();
          sprog.rigid_with.insert(sprog.index);

          tree[position].children.insert(tree.size()); //tells the current parent it has an extra child
          tree.insert(pair<unsigned int, obbranch>(tree.size(), sprog)); //adds the current branch to the tree

          rigid_fragments.erase(rigid_fragments.begin() + i);
          sterile = false;
          position = tree.size() - 1;
          break;
        }
      }
      if (sterile) {
        tree[position].done = true;
      }
    } else {
      position--;
    }
  }
}

bool FindBondedPiece(const vector<int>& root, const vector<int>& branched,
    unsigned int& root_atom, unsigned int& branch_atom,
    unsigned int& root_atom_rank, unsigned int& branch_atom_rank,
    const OBMol& mol, unsigned int & atoms_moved) {
  OBBond* the_bond;
  for (unsigned int i = 0; i < root.size(); i++) {
    for (unsigned int j = 0; j < branched.size(); j++) {
      the_bond = mol.GetBond(mol.GetAtom(root.at(i)),
          mol.GetAtom(branched.at(j)));
      if (the_bond != NULL) {
        root_atom = root.at(i);
        branch_atom = branched.at(j);
        root_atom_rank = i;
        branch_atom_rank = j;
        OBMol mol_copy = mol;
        the_bond = mol_copy.GetBond(mol_copy.GetAtom(root.at(i)),
            mol_copy.GetAtom(branched.at(j)));
        mol_copy.DeleteBond(the_bond, true);

        vector<vector<int> > two_pieces;
        mol_copy.ContigFragList(two_pieces);
        atoms_moved = two_pieces.at(1).size();
        return true;
      }
    }
  }
  return false;
}

bool DeleteHydrogens(OBMol & mol) {
  for (OBAtomIterator it = mol.BeginAtoms(); it != mol.EndAtoms(); it++) {
    if ((*it)->IsNonPolarHydrogen()) {
      OBBondIterator voider;
      double charger = (*it)->GetPartialCharge();
      charger += ((*it)->BeginNbrAtom(voider))->GetPartialCharge();
      ((*it)->BeginNbrAtom(voider))->SetPartialCharge(charger);
    }
  }
  return mol.DeleteNonPolarHydrogens();
}

//parsed_atoms store an index refering each atom to a position in atoms
//we need to construct an sdf context with the appropriate atoms in these positions
void createSDFContext(OBMol& mol, vector<OBAtom*> atoms, sdfcontext& sc) {
  sc.atoms.clear();
  sc.bonds.clear();
  sc.properties.clear();
  sc.name = mol.GetTitle();

  //setup mapping between atom indices (getIdx) and position in atoms
  boost::unordered_map<unsigned, unsigned> idx2atompos;
  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    idx2atompos[atoms[i]->GetIdx()] = i;
  }

  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    OBAtom *atom = atoms[i];
    const char *element_name = GET_SYMBOL(atom->GetAtomicNum());
    sc.atoms.push_back(sdfcontext::sdfatom(element_name));

    ///check for special properties
    if (atom->GetFormalCharge()) {
      sdfcontext::sdfprop prop(i, 'c', atom->GetFormalCharge());
      sc.properties.push_back(prop);
    }
    if (atom->GetIsotope()) {
      sdfcontext::sdfprop prop(i, 'i', atom->GetIsotope());
      sc.properties.push_back(prop);
    }
  }

  //now bonds
  for (OBMolBondIter bitr(mol); bitr; ++bitr) {
    unsigned first = idx2atompos[bitr->GetBeginAtomIdx()];
    unsigned second = idx2atompos[bitr->GetEndAtomIdx()];
    sc.bonds.push_back(sdfcontext::sdfbond(first, second, bitr->GetBondOrder()));
  }
}

static void OutputAtom(OBAtom* atom, context& lines, vector<OBAtom*>& atomorder,
    parsing_struct& p, const unsigned int index, unsigned immobile_num) {
  char buffer[BUFF_SIZE];
  char type_name[10], padded_name[10];
  char the_res[10];
  char the_chain = ' ';
  const char *element_name;
  string element_name_string;
  int res_num;
  bool het = false;
  stringstream ofs;

  OBResidue *res;
  strncpy(type_name, GET_SYMBOL(atom->GetAtomicNum()), sizeof(type_name));
  type_name[sizeof(type_name) - 1] = '\0';
  //two char. elements are on position 13 and 14 one char. start at 14

  if (strlen(type_name) > 1)
    type_name[1] = toupper(type_name[1]);
  else {
    char tmp[10];
    strncpy(tmp, type_name, 10);
    snprintf(type_name, sizeof(type_name), " %-3s", tmp);
  }

  if ((res = atom->GetResidue()) != 0) {
//			het = res->IsHetAtom(atom);
    snprintf(the_res, 4, "%s", (char*) res->GetName().c_str());
    snprintf(type_name, 5, "%s", (char*) res->GetAtomID(atom).c_str());
    the_chain = res->GetChain();

    //two char. elements are on position 13 and 14 one char. start at 14
    if (strlen(GET_SYMBOL(atom->GetAtomicNum())) == 1)
    {
      if (strlen(type_name) < 4)
      {
        char tmp[10];
        strncpy(tmp, type_name, 10);
        snprintf(padded_name, sizeof(padded_name), " %-3s", tmp);
        strncpy(type_name, padded_name, 4);
        type_name[4] = '\0';
      }
      else
      {
        type_name[4] = '\0';
      }
    }
    res_num = res->GetNum();
  } else {
    strcpy(the_res, "UNK");
    snprintf(padded_name, sizeof(padded_name), "%s", type_name);
    strncpy(type_name, padded_name, 4);
    type_name[4] = '\0';
    res_num = 1;
  }

  element_name = GET_SYMBOL(atom->GetAtomicNum());
  char element_name_final[3];
  element_name_final[2] = '\0';

  if (atom->GetAtomicNum() == 1) {
    element_name_final[0] = 'H';
    element_name_final[1] = 'D';
  } else
    if ((atom->GetAtomicNum() == 6) && (atom->IsAromatic())) {
      element_name_final[0] = 'A';
      element_name_final[1] = '\0';
    } else
      if (atom->GetAtomicNum() == 8) {
        element_name_final[0] = 'O';
        element_name_final[1] = 'A';
      } else
        if ((atom->GetAtomicNum() == 7) && (atom->IsHbondAcceptor())) {
          element_name_final[0] = 'N';
          element_name_final[1] = 'A';
        } else
          if ((atom->GetAtomicNum() == 16) && (atom->IsHbondAcceptor())) {
            element_name_final[0] = 'S';
            element_name_final[1] = 'A';
          } else {
            if (!isalnum(element_name[0])) {
              element_name_final[0] = '\0';
            } else
              element_name_final[0] = element_name[0];
            if (!isalnum(element_name[1])) {
              element_name_final[1] = '\0'; //null terminate
            } else
              element_name_final[1] = element_name[1];
          }

  double charge = atom->GetPartialCharge();
  snprintf(buffer, BUFF_SIZE,
      "%s%5d %-4s %-3s %c%3d     %8.3f%8.3f%8.3f  0.00  0.00    %+5.3f %.2s",
      het ? "HETATM" : "ATOM  ", index, type_name, the_res, the_chain, res_num,
      atom->GetX(), atom->GetY(), atom->GetZ(), charge, element_name_final);
  ofs << buffer;

  smt sm = string_to_smina_type(element_name_final);
  assert(sm < smina_atom_type::NumTypes);
  parsed_atom patom(sm, charge, vec(atom->GetX(), atom->GetY(), atom->GetZ()),
      index);
  //add_pdbqt_context(lines, ofs.str());
  if (patom.number == immobile_num) p.immobile_atom = p.atoms.size();
  p.add(patom, lines, atomorder.size());

  atomorder.push_back(atom);

}

static void OutputGroup(OBMol& mol, context& lines, vector<OBAtom*>& atomorder,
    parsing_struct& p, unsigned immobile_num, const vector<int>& group,
    map<unsigned int, unsigned int> new_indexes) {
  for (vector<int>::const_iterator it = group.begin(); it != group.end();
      it++) {
    OutputAtom(mol.GetAtom((*it)), lines, atomorder, p,
        new_indexes.find(*it)->second, immobile_num);
  }
}

//this we actually modify - we want the pdbqt output for the context, but
//also want to build up the parallel smina data structure for output
bool OutputTree(OBMol& mol, context& lines, parsing_struct& p,
    map<unsigned int, obbranch> & tree, unsigned int depth) {
  if (tree.size() == 0) {
    return false;
  }

  set_fixed_rotable_hydrogens(true);

  if (depth >= tree.size() - 1) {
    depth = tree.size() - 1;
  }

  vector<OBAtom*> atomorder;
  map<unsigned int, unsigned int> new_order; //gives the new ordering of the indexes of atoms, so that they are in increasing order from 1 in the output

  //generates the new ordering

  unsigned int current_atom_index = 1; //the index of the current atom
  for (unsigned int i = 0; i < tree.size(); i++) {
    assert(tree.count(i));
    for (set<unsigned int>::iterator it = tree[i].rigid_with.begin();
        it != tree[i].rigid_with.end(); it++) {
      vector<int> atoms = tree[*it].atoms;
      for (unsigned int j = 0; j < atoms.size(); j++) {
        new_order.insert(make_pair(atoms[j], current_atom_index));
        current_atom_index++;
      }
    }
  }

  stack<parsing_struct> pstack;
  pstack.push(parsing_struct());
  stack<pair<unsigned, unsigned> > bnumstack;

  //no legacy pdbqt context for smina output
  //add_pdbqt_context(lines, "ROOT");
  for (set<unsigned int>::iterator it = tree[0].rigid_with.begin();
      it != tree[0].rigid_with.end(); it++) {
    OutputGroup(mol, lines, atomorder, pstack.top(), INT_MAX, tree[*it].atoms,
        new_order);
  }

  //add_pdbqt_context(lines, "ENDROOT");

  for (unsigned int i = 1; i < tree.size(); i++) {
    stringstream ofs;
    ofs << "BRANCH";
    unsigned int parent_atom = tree[i].connecting_atom_parent;
    unsigned int child_atom = tree[i].connecting_atom_branch;

    ofs.width(4);
    unsigned parnum = (new_order.find(parent_atom))->second;
    ofs << parnum;
    ofs.width(4);
    unsigned childnum = (new_order.find(child_atom))->second;
    ofs << childnum;

    pstack.push(parsing_struct());
    bnumstack.push(make_pair(parnum, childnum)); //keep track of the parent number for this fragment
    //add_pdbqt_context(lines, ofs.str());

    for (set<unsigned int>::iterator it = tree[i].rigid_with.begin();
        it != tree[i].rigid_with.end(); it++) {
      OutputGroup(mol, lines, atomorder, pstack.top(), childnum,
          tree[*it].atoms, new_order);
    }

    for (vector<unsigned int>::iterator it = tree[i].parents.end();
        it != tree[i].parents.begin();) {
      it--;
      if ((*it) == 0) {
        break;
      } //do not close the main root; that is closed seperately
      vector<unsigned int>::iterator it_parent = it;
      it_parent--;
      if (tree[*it].children.size() == 0) {
        stringstream ofs;
        ofs << "ENDBRANCH";
        ofs.width(4);
        unsigned int parent_atom = tree[*it].connecting_atom_parent;
        unsigned int child_atom = tree[*it].connecting_atom_branch;
        ofs << new_order[parent_atom];
        ofs.width(4);
        ofs << new_order[child_atom];

        //finished with branch, so pop
        parsing_struct branch = pstack.top();
        pstack.pop();
        unsigned pnum = bnumstack.top().first;
        bnumstack.pop();
        //find parent frag
        unsigned pos = 0;
        for (unsigned n = pstack.top().atoms.size(); pos < n; pos++) {
          if (pstack.top().atoms[pos].a.number == pnum) break;
        }
        assert(pos < pstack.top().atoms.size());

        if (branch.mobile_hydrogens_only())
          pstack.top().mergeInto(branch);
        else
          pstack.top().atoms[pos].ps.push_back(branch);

        //add_pdbqt_context(lines, ofs.str());
        tree[*it_parent].children.erase(*it);
      }
    }
  }

  createSDFContext(mol, atomorder, lines.sdftext);
  assert(pstack.size() == 1);
  assert(bnumstack.size() == 0);
  p = pstack.top();
  return true;
}
