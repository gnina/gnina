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

#include "everything.h"

fl solvation_parameter(const atom_type& a) {
	if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).solvation;
	else if(a.xs == XS_TYPE_Met_D) return metal_solvation_parameter;
	VINA_CHECK(false);
	return 0; // placating the compiler
}

fl volume(const atom_type& a) {
	if(a.ad < AD_TYPE_SIZE) return ad_type_property(a.ad).volume;
	else if(a.xs < XS_TYPE_SIZE) return 4*pi / 3 * int_pow<3>(xs_radius(a.xs));
	VINA_CHECK(false);
	return 0; // placating the compiler
}

fl smooth_div(fl x, fl y) {
	if(std::abs(x) < epsilon_fl) return 0;
	if(std::abs(y) < epsilon_fl) return ((x*y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
	return x / y;
}

everything::everything() { // enabled according to design.out227
	const unsigned d = 0; // default
	const fl cutoff = 8; //6;

	// FIXME? enable some?
	//// distance_additive
	//add(d, new ad4_solvation(3.6, 0.01097,  true, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff
	//add(d, new ad4_solvation(3.6, 0.01097, false, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff

	//add(d, new electrostatic<1>(100, cutoff)); // cap, cutoff
	//add(d, new electrostatic<2>(100, cutoff)); // cap, cutoff

	//add(d, new gauss(0,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(0.5, 0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1.5, 0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2,   0.3, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2.5, 0.3, cutoff)); // offset, width, cutoff

	add(1, new gauss(0, 0.5, cutoff)); // offset, width, cutoff // WEIGHT: -0.035579
	//add(d, new gauss(1, 0.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.5, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 0.7, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 0.7, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.7, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 0.9, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 0.9, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 1.5, cutoff)); // offset, width, cutoff
	//add(d, new gauss(4, 1.5, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 2.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 2.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 2.0, cutoff)); // offset, width, cutoff
	add(1, new gauss(3, 2.0, cutoff)); // offset, width, cutoff // WEIGHT: -0.005156
	//add(d, new gauss(4, 2.0, cutoff)); // offset, width, cutoff

	//add(d, new gauss(0, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(1, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(2, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(3, 3.0, cutoff)); // offset, width, cutoff
	//add(d, new gauss(4, 3.0, cutoff)); // offset, width, cutoff

	//add(d, new repulsion( 0.4, cutoff)); // offset, cutoff
	//add(d, new repulsion( 0.2, cutoff)); // offset, cutoff
	add(1, new repulsion( 0.0, cutoff)); // offset, cutoff // WEIGHT:  0.840245
	//add(d, new repulsion(-0.2, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.4, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.6, cutoff)); // offset, cutoff
	//add(d, new repulsion(-0.8, cutoff)); // offset, cutoff
	//add(d, new repulsion(-1.0, cutoff)); // offset, cutoff

	//add(d, new hydrophobic(0.5, 1, cutoff)); // good, bad, cutoff
	add(1, new hydrophobic(0.5, 1.5, cutoff)); // good, bad, cutoff // WEIGHT:  -0.035069
	//add(d, new hydrophobic(0.5, 2, cutoff)); // good, bad, cutoff
	//add(d, new hydrophobic(0.5, 3, cutoff)); // good, bad, cutoff

	//add(1, new non_hydrophobic(0.5, 1.5, cutoff));

	//add(d, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff

	add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff // WEIGHT:  -0.587439
	//add(d, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff
	//add(d, new non_dir_h_bond(-0.7, 0.2, cutoff)); // good, bad, cutoff
	//add(d, new non_dir_h_bond(-0.7, 0.4, cutoff)); // good, bad, cutoff
	// additive

	// conf-independent
	//add(d, new num_ligands());

	add(1, new num_tors_div()); // WEIGHT: 1.923 -- FIXME too close to limit?
	//add(d, new num_heavy_atoms_div());
	//add(d, new num_heavy_atoms());
	//add(1, new num_tors_add());
	//add(d, new num_tors_sqr());
	//add(d, new num_tors_sqrt());
	//add(d, new num_hydrophobic_atoms());
	///add(1, new ligand_length());

	//add(d, new num_tors(100, 100, false)); // cap, past_cap, heavy_only
	//add(1, new num_tors(100, 100,  true)); // cap, past_cap, heavy_only
	//add(d, new num_tors(  2,   1,  true)); // cap, past_cap, heavy_only
	//add(d, new num_heavy_atoms());
	//add(d, new ligand_max_num_h_bonds());
	//add(1, new num_ligands());
}


dkoes_terms::dkoes_terms() { // chosen through linear regression and cross validation
	const unsigned d = 0; // default
	const fl cutoff = 8; //6;



	add(1, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff
	add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff // WEIGHT:  -0.587439

	add(1, new ad4_solvation(3.6, 0.01097,  true, cutoff)); // desolvation_sigma, solvation_q, charge_dependent, cutoff

	add(1, new num_tors_sqr());
	add(1, new constant_term());
}

old_dkoes_terms::old_dkoes_terms() { // chosen through linear regression and cross validation
	const unsigned d = 0; // default
	const fl cutoff = 8;

	add(1, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff
	add(1, new non_dir_h_bond(-0.7, 0, cutoff)); // good, bad, cutoff
	add(1, new num_tors_sqr());
	add(1, new constant_term());
}

fast_dkoes_terms::fast_dkoes_terms() { // chosen through linear regression and cross validation
	const unsigned d = 0; // default
	const fl cutoff = 8;
	add(1, new vdw<4,  8>(   0, 100, cutoff)); // smoothing, cap, cutoff
	add(1, new num_tors_sqr());
	add(1, new constant_term());
}
