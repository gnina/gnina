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

#ifndef VINA_ATOM_CONSTANTS_H
#define VINA_ATOM_CONSTANTS_H

#include "common.h"
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/babelconfig.h>

#if (OB_VERSION >= OB_VERSION_CHECK(2,4,90))
# include <openbabel/elements.h>
# define GET_SYMBOL OpenBabel::OBElements::GetSymbol
# define GET_HVY(a) a->GetHvyDegree()
#else
# define GET_SYMBOL etab.GetSymbol
# define GET_HVY(a) a->GetHvyValence()
#endif

//SMINA unified atom types - these must represent all possible combinations of autodock and x-scale atom types

namespace smina_atom_type
{
  enum type {
    /* 0 */Hydrogen, // H_H_X,
    /* 1 */PolarHydrogen,//(can donate) H_HD_X,
    /* 2 */AliphaticCarbonXSHydrophobe,// C_C_C_H,  //hydrophobic according to xscale
    /* 3 */AliphaticCarbonXSNonHydrophobe,//C_C_C_P, //not hydrophobic (according to xs)
    /* 4 */AromaticCarbonXSHydrophobe,//C_A_C_H,
    /* 5 */AromaticCarbonXSNonHydrophobe,//C_A_C_P,
    /* 6 */Nitrogen,//N_N_N_P, no hydrogen bonding
    /* 7 */NitrogenXSDonor,//N_N_N_D,
    /* 8 */NitrogenXSDonorAcceptor,//N_NA_N_DA, also an autodock acceptor
    /* 9 */NitrogenXSAcceptor,//N_NA_N_A, also considered an acceptor by autodock
    /* 10 */Oxygen,//O_O_O_P,
    /* 11 */OxygenXSDonor,//O_O_O_D,
    /* 12 */OxygenXSDonorAcceptor,//O_OA_O_DA, also an autodock acceptor
    /* 13 */OxygenXSAcceptor,//O_OA_O_A, also an autodock acceptor
    /* 14 */Sulfur,//S_S_S_P,
    /* 15 */SulfurAcceptor,//S_SA_S_P, XS doesn't do sulfur acceptors
    /* 16 */Phosphorus,//P_P_P_P,
    /* 17 */Fluorine,//F_F_F_H,
    /* 18 */Chlorine,//Cl_Cl_Cl_H,
    /* 19 */Bromine,//Br_Br_Br_H,
    /* 20 */Iodine,//I_I_I_H,
    /* 21 */Magnesium,//Met_Mg_Met_D,
    /* 22 */Manganese,//Met_Mn_Met_D,
    /* 23 */Zinc,// Met_Zn_Met_D,
    /* 24 */Calcium,//Met_Ca_Met_D,
    /* 25 */Iron,//Met_Fe_Met_D,
    /* 26 */GenericMetal,//Met_METAL_Met_D,
    /* 27 */Boron,//there are 160 cmpds in pdbbind (general, not refined) with boron
    NumTypes
  };

//store all the desired properties in smina_type_info
  struct info
  {
    type sm;
    const char* smina_name; //this must be more than 2 chars long
    const char* adname;//this must be no longer than 2 chars
    sz anum;
    fl ad_radius;
    fl ad_depth;
    fl ad_solvation;
    fl ad_volume;
    fl covalent_radius;
    fl xs_radius;
    bool xs_hydrophobe;
    bool xs_donor;
    bool xs_acceptor;
    bool ad_heteroatom;

  };

  extern info data[NumTypes];

//dkoes - eventually I'd like to switch to a single unified atom typing, but
//for now stitch everything together with a smina atom type
  const info default_data[NumTypes] = { //el, ad, xs
    { Hydrogen,"Hydrogen", "H", 1, 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.37, false, false, false, false},
    { PolarHydrogen, "PolarHydrogen", "HD", 1, 1.000000, 0.020000, 0.000510, 0.000000, 0.370000, 0.370000, false, false, false, false},
    //note we typically use the xs_radius, which assumes a heavy atom-only model
    { AliphaticCarbonXSHydrophobe, "AliphaticCarbonXSHydrophobe", "C", 6, 2.000000, 0.150000, -0.001430, 33.510300, 0.770000, 1.900000, true, false, false, false},
    { AliphaticCarbonXSNonHydrophobe, "AliphaticCarbonXSNonHydrophobe", "C", 6, 2.000000, 0.150000, -0.001430, 33.510300, 0.770000, 1.900000, false, false, false, false},
    { AromaticCarbonXSHydrophobe, "AromaticCarbonXSHydrophobe", "A", 6, 2.000000, 0.150000, -0.000520, 33.510300, 0.770000, 1.900000, true, false, false, false},
    { AromaticCarbonXSNonHydrophobe, "AromaticCarbonXSNonHydrophobe", "A", 6, 2.000000, 0.150000, -0.000520, 33.510300, 0.770000, 1.900000, false, false, false, false},
    { Nitrogen, "Nitrogen", "N", 7, 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, false, false, true},
    { NitrogenXSDonor, "NitrogenXSDonor", "N", 7, 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, true, false, true},
    { NitrogenXSDonorAcceptor, "NitrogenXSDonorAcceptor", "NA", 7, 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, true, true, true},
    { NitrogenXSAcceptor, "NitrogenXSAcceptor", "NA", 7, 1.750000, 0.160000, -0.001620, 22.449300, 0.750000, 1.800000, false, false, true, true},
    { Oxygen, "Oxygen", "O", 8, 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, false, false, true},
    { OxygenXSDonor, "OxygenXSDonor", "O", 8, 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, true, false, true},
    { OxygenXSDonorAcceptor, "OxygenXSDonorAcceptor", "OA", 8, 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, true, true, true},
    { OxygenXSAcceptor, "OxygenXSAcceptor", "OA", 8, 1.600000, 0.200000, -0.002510, 17.157300, 0.730000, 1.700000, false, false, true, true},
    { Sulfur, "Sulfur", "S", 16, 2.000000, 0.200000, -0.002140, 33.510300, 1.020000, 2.000000, false, false, false, true},
    { SulfurAcceptor, "SulfurAcceptor", "SA", 16, 2.000000, 0.200000, -0.002140, 33.510300, 1.020000, 2.000000, false, false, false, true},
    { Phosphorus, "Phosphorus", "P", 15, 2.100000, 0.200000, -0.001100, 38.792400, 1.060000, 2.100000, false, false, false, true},
    { Fluorine, "Fluorine", "F", 9, 1.545000, 0.080000, -0.001100, 15.448000, 0.710000, 1.500000, true, false, false, true},
    { Chlorine, "Chlorine", "Cl", 17, 2.045000, 0.276000, -0.001100, 35.823500, 0.990000, 1.800000, true, false, false, true},
    { Bromine, "Bromine", "Br", 35, 2.165000, 0.389000, -0.001100, 42.566100, 1.140000, 2.000000, true, false, false, true},
    { Iodine, "Iodine", "I", 53, 2.360000, 0.550000, -0.001100, 55.058500, 1.330000, 2.200000, true, false, false, true},
    { Magnesium, "Magnesium", "Mg", 12, 0.650000, 0.875000, -0.001100, 1.560000, 1.300000, 1.200000, false, true, false, true},
    { Manganese, "Manganese", "Mn", 25, 0.650000, 0.875000, -0.001100, 2.140000, 1.390000, 1.200000, false, true, false, true},
    { Zinc, "Zinc", "Zn", 30, 0.740000, 0.550000, -0.001100, 1.700000, 1.310000, 1.200000, false, true, false, true},
    { Calcium, "Calcium", "Ca", 20, 0.990000, 0.550000, -0.001100, 2.770000, 1.740000, 1.200000, false, true, false, true},
    { Iron, "Iron", "Fe", 26, 0.650000, 0.010000, -0.001100, 1.840000, 1.250000, 1.200000, false, true, false, true},
    { GenericMetal, "GenericMetal", "M", 0, 1.200000, 0.000000, -0.001100, 22.449300, 1.750000, 1.200000, false, true, false, true},
    //note AD4 doesn't have boron, so copying from carbon
    { Boron, "Boron", "B", 5, 2.04, 0.180000, -0.0011, 12.052, 0.90, 1.920000, true, false, false, false}

  };

}

typedef smina_atom_type::type smt;

struct atom_equivalence {
    std::string name;
    std::string to;
};

const atom_equivalence atom_equivalence_data[] = { { "Se", "S" } };

const sz atom_equivalences_size = sizeof(atom_equivalence_data)
    / sizeof(const atom_equivalence);

inline bool is_hydrogen(smt t) {
  return t == smina_atom_type::Hydrogen || t == smina_atom_type::PolarHydrogen;
}

inline bool is_heteroatom(smt t) {
  assert(t < smina_atom_type::NumTypes);
  return smina_atom_type::data[t].ad_heteroatom;
}

inline fl covalent_radius(const smt t) {
  assert(t < smina_atom_type::NumTypes);
  return smina_atom_type::data[t].covalent_radius;
}

inline fl xs_radius(smina_atom_type::type t) {
  assert(t < smina_atom_type::NumTypes);
  return smina_atom_type::data[t].xs_radius;
}

const std::string non_ad_metal_names[] = { // expand as necessary
    "Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni", "Si" };

inline bool is_non_ad_metal_name(const std::string& name) {
  const sz s = sizeof(non_ad_metal_names) / sizeof(const std::string);
  VINA_FOR(i, s)
    if (non_ad_metal_names[i] == name) return true;
  return false;
}

inline bool xs_is_hydrophobic(smt sm) {
  assert(sm < smina_atom_type::NumTypes);
  return smina_atom_type::data[sm].xs_hydrophobe;
}

inline bool xs_is_acceptor(smt sm) {
  assert(sm < smina_atom_type::NumTypes);
  return smina_atom_type::data[sm].xs_acceptor;
}

inline bool xs_is_donor(smt sm) {
  assert(sm < smina_atom_type::NumTypes);
  return smina_atom_type::data[sm].xs_donor;
}

//only checks one order! should only be used as a helper for
//xs_h_bond_possible
inline bool xs_donor_acceptor(smt t1, smt t2) {
  return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_h_bond_possible(smt t1, smt t2) {
  return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}

//return true if both types are strictly donors or both strictly acceptors
inline bool xs_anti_h_bond(smt t1, smt t2) {
  if (xs_is_donor(t1) && !xs_is_acceptor(t1)) {
    return xs_is_donor(t2) && !xs_is_acceptor(t2);
  }
  if (!xs_is_donor(t1) && xs_is_acceptor(t1)) {
    return !xs_is_donor(t2) && xs_is_acceptor(t2);
  }
  return false;
}

inline const char* smina_type_to_string(smt at) {
  return smina_atom_type::data[at].smina_name;
}

inline std::string smina_type_to_element_name(smt at) {
  //try to figure out from adname, this shouldn't necessarily be relied on, but is useful for debug output (e.g. xyz)
  std::string ret = smina_atom_type::data[at].adname;
  if (ret == "A") {
    return "C";
  } else
    if (ret.back() == 'A' || ret.back() == 'D') {
      ret.pop_back();
    }
  return ret;
}

inline smt string_to_smina_type(const std::string& name) {
  //dkoes - returns NumTypes if can't identify the type
  // if name is 2 chars or less, assume it is an AD4 type,
  //otherwise assume it is a full smina type name
  //I'm assuming this will not be called frequently and so am not using lookup tables
  if (name.length() <= 2) {
    VINA_FOR(i, smina_atom_type::NumTypes)
      if (smina_atom_type::data[i].adname == name)
        return smina_atom_type::data[i].sm;
    VINA_FOR(i, atom_equivalences_size)
      if (atom_equivalence_data[i].name == name)
        return string_to_smina_type(atom_equivalence_data[i].to);
    if (is_non_ad_metal_name(name)) return smina_atom_type::GenericMetal; //generic metal
    return smina_atom_type::GenericMetal; //TODO: implement default catch-all type

  } else {
    VINA_FOR(i, smina_atom_type::NumTypes)
      if (smina_atom_type::data[i].smina_name == name)
        return smina_atom_type::data[i].sm;
    return smina_atom_type::NumTypes;
  }
}

inline fl max_covalent_radius() {
  fl tmp = 0;
  VINA_FOR(i, smina_atom_type::NumTypes)
    if (smina_atom_type::data[i].covalent_radius > tmp) tmp =
        smina_atom_type::data[i].covalent_radius;
  return tmp;
}

inline fl solvation_parameter(smt t) {
  assert(t < smina_atom_type::NumTypes);
  return smina_atom_type::data[t].ad_solvation;
}

inline fl ad_volume(smt t) {
  assert(t < smina_atom_type::NumTypes);
  return smina_atom_type::data[t].ad_volume;
}

inline fl ad_depth(smt t) {
  assert(t < smina_atom_type::NumTypes);
  return smina_atom_type::data[t].ad_depth;
}

//take atom's neighborhood into account when setting atom type
//IMPORTANT: make sure this is conistent with obatom_to_smina_type
inline smt adjust_smina_type(smt t, bool Hbonded, bool heteroBonded) {
  using namespace smina_atom_type;
  switch (t) {
  case AliphaticCarbonXSHydrophobe: // C_C_C_H, //hydrophobic according to xscale
  case AliphaticCarbonXSNonHydrophobe: //C_C_C_P,
    return
        heteroBonded ?
            AliphaticCarbonXSNonHydrophobe : AliphaticCarbonXSHydrophobe;
  case AromaticCarbonXSHydrophobe: //C_A_C_H,
  case AromaticCarbonXSNonHydrophobe: //C_A_C_P,
    return
        heteroBonded ?
            AromaticCarbonXSNonHydrophobe : AromaticCarbonXSHydrophobe;
  case NitrogenXSDonor: //N_N_N_D,
  case Nitrogen: //N_N_N_P, no hydrogen bonding
    return Hbonded ? NitrogenXSDonor : Nitrogen;
  case NitrogenXSDonorAcceptor: //N_NA_N_DA, also an autodock acceptor
  case NitrogenXSAcceptor: //N_NA_N_A, also considered an acceptor by autodock
    return Hbonded ? NitrogenXSDonorAcceptor : NitrogenXSAcceptor;
  case OxygenXSDonor: //O_O_O_D,
  case Oxygen: //O_O_O_P,
    return Hbonded ? OxygenXSDonor : Oxygen;
  case OxygenXSDonorAcceptor: //O_OA_O_DA, also an autodock acceptor
  case OxygenXSAcceptor: //O_OA_O_A, also an autodock acceptor
    return Hbonded ? OxygenXSDonorAcceptor : OxygenXSAcceptor;
  default:
    return t;
  }

}

//return smina atom type of provided atom; this duplicates
//the parsing code and above adjust_smina_type code, but including it here
//let's us to atom typing without a link dependency
//IMPORTANT: make sure this is consistent with adjust_smina_type and pdbqt parsing
inline smt obatom_to_smina_type(OpenBabel::OBAtom& atom) {
  using namespace OpenBabel;
  //from pdbqt format
  const char *element_name = GET_SYMBOL(atom.GetAtomicNum());
  std::string ename(element_name);

  if (atom.GetAtomicNum() == 1)
    ename = "HD";
  else
    if ((atom.GetAtomicNum() == 6) && (atom.IsAromatic()))
      ename = "A";
    else
      if (atom.GetAtomicNum() == 8)
        ename = "OA";
      else
        if ((atom.GetAtomicNum() == 7) && (atom.IsHbondAcceptor()))
          ename = "NA";
        else
          if ((atom.GetAtomicNum() == 16) && (atom.IsHbondAcceptor())) ename =
              "SA";

  smt atype = string_to_smina_type(ename);

  bool hbonded = false;
  bool heteroBonded = false;

  FOR_NBORS_OF_ATOM(neigh, atom){
  if(neigh->GetAtomicNum() == 1)
  hbonded = true;
  else if(neigh->GetAtomicNum() != 6)
  heteroBonded = true; //hetero anything that is not hydrogen and not carbon
}

  return adjust_smina_type(atype, hbonded, heteroBonded);
}

#endif
