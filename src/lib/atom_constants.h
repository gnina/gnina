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

// based on SY_TYPE_* but includes H
const sz EL_TYPE_H    =  0;
const sz EL_TYPE_C    =  1;
const sz EL_TYPE_N    =  2;
const sz EL_TYPE_O    =  3;
const sz EL_TYPE_S    =  4;
const sz EL_TYPE_P    =  5;
const sz EL_TYPE_F    =  6;
const sz EL_TYPE_Cl   =  7;
const sz EL_TYPE_Br   =  8;
const sz EL_TYPE_I    =  9;
const sz EL_TYPE_Met  = 10;
const sz EL_TYPE_SIZE = 11;

// AutoDock4
const sz AD_TYPE_C    =  0;
const sz AD_TYPE_A    =  1;
const sz AD_TYPE_N    =  2;
const sz AD_TYPE_O    =  3;
const sz AD_TYPE_P    =  4;
const sz AD_TYPE_S    =  5;
const sz AD_TYPE_H    =  6; // non-polar hydrogen
const sz AD_TYPE_F    =  7;
const sz AD_TYPE_I    =  8;
const sz AD_TYPE_NA   =  9;
const sz AD_TYPE_OA   = 10;
const sz AD_TYPE_SA   = 11;
const sz AD_TYPE_HD   = 12;
const sz AD_TYPE_Mg   = 13;
const sz AD_TYPE_Mn   = 14;
const sz AD_TYPE_Zn   = 15;
const sz AD_TYPE_Ca   = 16;
const sz AD_TYPE_Fe   = 17;
const sz AD_TYPE_Cl   = 18;
const sz AD_TYPE_Br   = 19;
const sz AD_TYPE_METAL = 20; //generic metal, not actually part of autodock
const sz AD_TYPE_SIZE = 21;

// X-Score
const sz XS_TYPE_C_H   =  0;
const sz XS_TYPE_C_P   =  1;
const sz XS_TYPE_N_P   =  2;
const sz XS_TYPE_N_D   =  3;
const sz XS_TYPE_N_A   =  4;
const sz XS_TYPE_N_DA  =  5;
const sz XS_TYPE_O_P   =  6;
const sz XS_TYPE_O_D   =  7;
const sz XS_TYPE_O_A   =  8;
const sz XS_TYPE_O_DA  =  9;
const sz XS_TYPE_S_P   = 10;
const sz XS_TYPE_P_P   = 11;
const sz XS_TYPE_F_H   = 12;
const sz XS_TYPE_Cl_H  = 13;
const sz XS_TYPE_Br_H  = 14;
const sz XS_TYPE_I_H   = 15;
const sz XS_TYPE_Met_D = 16;
const sz XS_TYPE_SIZE  = 17;



//SMINA unified atom types - these must represent all possible combinations of autodock and x-scale atom types

namespace smina_atom_type
{
enum type {
	Hydrogen, // H_H_X,
	PolarHydrogen, //(can donate) H_HD_X,
	AliphaticCarbonXSHydrophobe, // C_C_C_H, //hydrophobic according to xscale
	AliphaticCarbonXSNonHydrophobe, //C_C_C_P, //not hydrophobic (according to xs)
	AromaticCarbonXSHydrophobe, //C_A_C_H,
	AromaticCarbonXSNonHydrophobe, //C_A_C_P,
	Nitrogen, //N_N_N_P, no hydrogen bonding
	NitrogenXSDonor, //N_N_N_D,
	NitrogenXSDonorAcceptor, //N_NA_N_DA, also an autodock acceptor
	NitrogenXSAcceptor, //N_NA_N_A, also considered an acceptor by autodock
	Oxygen, //O_O_O_P,
	OxygenXSDonor, //O_O_O_D,
	OxygenXSDonorAcceptor, //O_OA_O_DA, also an autodock acceptor
	OxygenXSAcceptor, //O_OA_O_A, also an autodock acceptor
	Sulfur, //S_S_S_P,
	SulfurAcceptor, //S_SA_S_P, XS doesn't do sulfur acceptors
	Phosphorus, //P_P_P_P,
	Fluorine, //F_F_F_H,
	Chlorine, //Cl_Cl_Cl_H,
	Bromine, //Br_Br_Br_H,
	Iodine, //I_I_I_H,
	Magnesium, //Met_Mg_Met_D,
	Manganese, //Met_Mn_Met_D,
	Zinc, // Met_Zn_Met_D,
	Calcium, //Met_Ca_Met_D,
	Iron, //Met_Fe_Met_D,
	GenericMetal, //Met_METAL_Met_D,
	NumTypes
};

//store all the desired properties in smina_type_info
struct info
{
	type sm;
	sz el;
	sz ad;
	sz xs;
	const char* smina_name; //this must be more than 2 chars long
	const char* adname; //this must be no longer than 2 chars
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
		{Hydrogen, EL_TYPE_H, AD_TYPE_H, XS_TYPE_SIZE,"Hydrogen",
				"H",	1.000000,	0.020000,	0.000510,	0.000000,	0.370000,	0.000000,	false,	false,	false,	false},
		{PolarHydrogen, EL_TYPE_H, AD_TYPE_HD, XS_TYPE_SIZE,"PolarHydrogen",
				"HD",	1.000000,	0.020000,	0.000510,	0.000000,	0.370000,	0.000000,	false,	false,	false,	false},
		{AliphaticCarbonXSHydrophobe, EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_H,"AliphaticCarbonXSHydrophobe",
				"C",	2.000000,	0.150000,	-0.001430,	33.510300,	0.770000,	1.900000,	true,	false,	false,	false},
		{AliphaticCarbonXSNonHydrophobe, EL_TYPE_C, AD_TYPE_C, XS_TYPE_C_P,"AliphaticCarbonXSNonHydrophobe",
				"C",	2.000000,	0.150000,	-0.001430,	33.510300,	0.770000,	1.900000,	false,	false,	false,	false},
		{AromaticCarbonXSHydrophobe, EL_TYPE_C, AD_TYPE_A, XS_TYPE_C_H,"AromaticCarbonXSHydrophobe",
				"A",	2.000000,	0.150000,	-0.000520,	33.510300,	0.770000,	1.900000,	true,	false,	false,	false},
		{AromaticCarbonXSNonHydrophobe, EL_TYPE_C, AD_TYPE_A, XS_TYPE_C_P,"AromaticCarbonXSNonHydrophobe",
				"A",	2.000000,	0.150000,	-0.000520,	33.510300,	0.770000,	1.900000,	false,	false,	false,	false},
		{Nitrogen, EL_TYPE_N, AD_TYPE_N, XS_TYPE_N_P,"Nitrogen",
				"N",	1.750000,	0.160000,	-0.001620,	22.449300,	0.750000,	1.800000,	false,	false,	false,	true},
		{NitrogenXSDonor, EL_TYPE_N, AD_TYPE_N, XS_TYPE_N_D,"NitrogenXSDonor",
				"N",	1.750000,	0.160000,	-0.001620,	22.449300,	0.750000,	1.800000,	false,	true,	false,	true},
		{NitrogenXSDonorAcceptor, EL_TYPE_N, AD_TYPE_NA, XS_TYPE_N_DA,"NitrogenXSDonorAcceptor",
				"NA",	1.750000,	0.160000,	-0.001620,	22.449300,	0.750000,	1.800000,	false,	true,	true,	true},
		{NitrogenXSAcceptor, EL_TYPE_N, AD_TYPE_NA, XS_TYPE_N_A,"NitrogenXSAcceptor",
				"NA",	1.750000,	0.160000,	-0.001620,	22.449300,	0.750000,	1.800000,	false,	false,	true,	true},
		{Oxygen, EL_TYPE_O, AD_TYPE_O, XS_TYPE_O_P,"Oxygen",
				"O",	1.600000,	0.200000,	-0.002510,	17.157300,	0.730000,	1.700000,	false,	false,	false,	true},
		{OxygenXSDonor, EL_TYPE_O, AD_TYPE_O, XS_TYPE_O_D,"OxygenXSDonor",
				"O",	1.600000,	0.200000,	-0.002510,	17.157300,	0.730000,	1.700000,	false,	true,	false,	true},
		{OxygenXSDonorAcceptor, EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_DA,"OxygenXSDonorAcceptor",
				"OA",	1.600000,	0.200000,	-0.002510,	17.157300,	0.730000,	1.700000,	false,	true,	true,	true},
		{OxygenXSAcceptor, EL_TYPE_O, AD_TYPE_OA, XS_TYPE_O_A,"OxygenXSAcceptor",
				"OA",	1.600000,	0.200000,	-0.002510,	17.157300,	0.730000,	1.700000,	false,	false,	true,	true},
		{Sulfur, EL_TYPE_S, AD_TYPE_S, XS_TYPE_S_P,"Sulfur",
				"S",	2.000000,	0.200000,	-0.002140,	33.510300,	1.020000,	2.000000,	false,	false,	false,	true},
		{SulfurAcceptor, EL_TYPE_S, AD_TYPE_SA, XS_TYPE_S_P,"SulfurAcceptor",
				"SA",	2.000000,	0.200000,	-0.002140,	33.510300,	1.020000,	2.000000,	false,	false,	false,	true},
		{Phosphorus, EL_TYPE_P, AD_TYPE_P, XS_TYPE_P_P,"Phosphorus",
				"P",	2.100000,	0.200000,	-0.001100,	38.792400,	1.060000,	2.100000,	false,	false,	false,	true},
		{Fluorine, EL_TYPE_F, AD_TYPE_F, XS_TYPE_F_H,"Fluorine",
				"F",	1.545000,	0.080000,	-0.001100,	15.448000,	0.710000,	1.500000,	true,	false,	false,	true},
		{Chlorine, EL_TYPE_Cl, AD_TYPE_Cl, XS_TYPE_Cl_H,"Chlorine",
				"Cl",	2.045000,	0.276000,	-0.001100,	35.823500,	0.990000,	1.800000,	true,	false,	false,	true},
		{Bromine, EL_TYPE_Br, AD_TYPE_Br, XS_TYPE_Br_H,"Bromine",
				"Br",	2.165000,	0.389000,	-0.001100,	42.566100,	1.140000,	2.000000,	true,	false,	false,	true},
		{Iodine, EL_TYPE_I, AD_TYPE_I, XS_TYPE_I_H,"Iodine",
				"I",	2.360000,	0.550000,	-0.001100,	55.058500,	1.330000,	2.200000,	true,	false,	false,	true},
		{Magnesium, EL_TYPE_Met, AD_TYPE_Mg, XS_TYPE_Met_D,"Magnesium",
				"Mg",	0.650000,	0.875000,	-0.001100,	1.560000,	1.300000,	1.200000,	false,	true,	false,	true},
		{Manganese, EL_TYPE_Met, AD_TYPE_Mn, XS_TYPE_Met_D,"Manganese",
				"Mn",	0.650000,	0.875000,	-0.001100,	2.140000,	1.390000,	1.200000,	false,	true,	false,	true},
		{Zinc, EL_TYPE_Met, AD_TYPE_Zn, XS_TYPE_Met_D,"Zinc",
				"Zn",	0.740000,	0.550000,	-0.001100,	1.700000,	1.310000,	1.200000,	false,	true,	false,	true},
		{Calcium, EL_TYPE_Met, AD_TYPE_Ca, XS_TYPE_Met_D,"Calcium",
				"Ca",	0.990000,	0.550000,	-0.001100,	2.770000,	1.740000,	1.200000,	false,	true,	false,	true},
		{Iron, EL_TYPE_Met, AD_TYPE_Fe, XS_TYPE_Met_D,"Iron",
				"Fe",	0.650000,	0.010000,	-0.001100,	1.840000,	1.250000,	1.200000,	false,	true,	false,	true},
		{GenericMetal, EL_TYPE_Met, AD_TYPE_METAL, XS_TYPE_Met_D,"GenericMetal",
				"M",	1.200000,	0.000000,	-0.001100,	22.449300,	1.750000,	1.200000,	false,	true,	false,	true}
};

}

typedef smina_atom_type::type smt;

struct atom_equivalence {
	std::string name;
	std::string to;
};

const atom_equivalence atom_equivalence_data[] = {
	{"Se",  "S"}
};

const sz atom_equivalences_size = sizeof(atom_equivalence_data) / sizeof(const atom_equivalence);


inline bool is_hydrogen(smt t) {
	return t == smina_atom_type::Hydrogen || t == smina_atom_type::PolarHydrogen;
}

inline bool is_heteroatom(smt t) {
	assert(t < smina_atom_type::NumTypes);
	return smina_atom_type::data[t].ad_heteroatom;
}

inline fl covalent_radius(const smt t)
{
	assert(t < smina_atom_type::NumTypes);
	return smina_atom_type::data[t].covalent_radius;
}


inline fl xs_radius(smina_atom_type::type t) {
	assert(t < smina_atom_type::NumTypes);
	return smina_atom_type::data[t].xs_radius;
}

const std::string non_ad_metal_names[] = { // expand as necessary
	"Cu", "Fe", "Na", "K", "Hg", "Co", "U", "Cd", "Ni"
};

inline bool is_non_ad_metal_name(const std::string& name) {
	const sz s = sizeof(non_ad_metal_names) / sizeof(const std::string);
	VINA_FOR(i, s)
		if(non_ad_metal_names[i] == name)
			return true;
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
	if(xs_is_donor(t1) && !xs_is_acceptor(t1)) {
		return xs_is_donor(t2) && !xs_is_acceptor(t2);
	}
	if(!xs_is_donor(t1) && xs_is_acceptor(t1)) {
		return !xs_is_donor(t2) && xs_is_acceptor(t2);
	}
	return false;
}


inline smt string_to_smina_type(const std::string& name)
{
	//dkoes - returns NumTypes if can't identify the type
	// if name is 2 chars or less, assume it is an AD4 type,
	//otherwise assume it is a full smina type name
	//I'm assuming this will not be called frequently and so am using lookup tables
	if (name.length() <= 2)
	{
		VINA_FOR(i, smina_atom_type::NumTypes)
			if (smina_atom_type::data[i].adname == name)
				return smina_atom_type::data[i].sm;
		VINA_FOR(i, atom_equivalences_size)
			if (atom_equivalence_data[i].name == name)
				return string_to_smina_type(atom_equivalence_data[i].to);
		if (is_non_ad_metal_name(name))
			return smina_atom_type::GenericMetal; //generic metal
		return smina_atom_type::NumTypes;
	}
	else
	{
		VINA_FOR(i, smina_atom_type::NumTypes)
			if(smina_atom_type::data[i].smina_name == name)
				return smina_atom_type::data[i].sm;
		return smina_atom_type::NumTypes;
	}
}

inline fl max_covalent_radius() {
	fl tmp = 0;
	VINA_FOR(i, smina_atom_type::NumTypes)
		if(smina_atom_type::data[i].covalent_radius > tmp)
			tmp = smina_atom_type::data[i].covalent_radius;
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




inline smt adjust_smina_type(smt t, bool Hbonded, bool heteroBonded)
{
	using namespace smina_atom_type;
	switch(t)
	{
	case AliphaticCarbonXSHydrophobe: // C_C_C_H, //hydrophobic according to xscale
	case AliphaticCarbonXSNonHydrophobe: //C_C_C_P,
		return heteroBonded ? AliphaticCarbonXSNonHydrophobe : AliphaticCarbonXSHydrophobe;
	case AromaticCarbonXSHydrophobe: //C_A_C_H,
	case AromaticCarbonXSNonHydrophobe: //C_A_C_P,
		return heteroBonded ? AromaticCarbonXSNonHydrophobe : AromaticCarbonXSHydrophobe;
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

#endif
