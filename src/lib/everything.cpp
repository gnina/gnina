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


fl ad4_solvation::eval(const atom_base& a, const atom_base& b, fl r) const {
	fl q1 = a.charge;
	fl q2 = b.charge;

	VINA_CHECK(not_max(q1));
	VINA_CHECK(not_max(q2));

	sz t1 = a.ad;
	sz t2 = b.ad;

	fl solv1 = solvation_parameter(a);
	fl solv2 = solvation_parameter(b);

	fl volume1 = volume(a);
	fl volume2 = volume(b);

	fl my_solv = charge_dependent ? solvation_q : 0;

	fl tmp = ((solv1 + my_solv * std::abs(q1)) * volume2 +
		    (solv2 + my_solv * std::abs(q2)) * volume1) * std::exp(-sqr(r/(2*desolvation_sigma)));

	VINA_CHECK(not_max(tmp));
	return tmp;
}
