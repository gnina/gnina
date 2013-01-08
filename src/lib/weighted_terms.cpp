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

#include "weighted_terms.h"

//dkoes - FIX: terms and weights must be in a specific order (usable, da, const)
weighted_terms::weighted_terms(const terms* t, const flv& weights) : t(t), weights(weights), cutoff_(0),conf_indep_start(0) { // does not own t
	VINA_CHECK(t->         additive_terms.num_enabled() == 0);
	VINA_CHECK(t->   intermolecular_terms.num_enabled() == 0);
	VINA_FOR_IN(i, t->usable_terms)
		if(t->usable_terms.enabled[i]) {
			enabled_usable_terms.push_back(i);
			cutoff_ = (std::max)(cutoff_, t->usable_terms[i].cutoff);
		}

	VINA_FOR_IN(i, t->distance_additive_terms)
		if(t->distance_additive_terms.enabled[i]) {
			enabled_distance_additive_terms.push_back(i);
			cutoff_ = (std::max)(cutoff_, t->distance_additive_terms[i].cutoff);
		}

	conf_indep_start = enabled_usable_terms.size() + enabled_distance_additive_terms.size();
}

//dkoes - evaluate usable (atom type) terms only
fl weighted_terms::eval_fast(smt t1, smt t2, fl r) const { // intentionally not checking for cutoff
	fl acc = 0;
	VINA_FOR_IN(i, enabled_usable_terms)
		acc += weights[i] * t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r);
	return acc;
}

//dkoes - evalute da terms here
fl weighted_terms::eval_slow(const atom_base& a, const atom_base& b, fl r) const {
	fl acc = 0;
	sz offset = enabled_usable_terms.size();
	VINA_FOR_IN(i, enabled_distance_additive_terms)
		acc += weights[offset+i] * t->distance_additive_terms[enabled_distance_additive_terms[i]].eval(a, b, r);
	return acc;
}

fl weighted_terms::conf_independent(const model& m, fl e) const {
	flv::const_iterator it = weights.begin() + conf_indep_start;
	conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?
	fl tmp = t->eval_conf_independent(in, e, it);
	assert(it == weights.end());
	return tmp;
}
