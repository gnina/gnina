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

#ifndef VINA_WEIGHTED_TERMS_H
#define VINA_WEIGHTED_TERMS_H

#include "terms.h"

struct weighted_terms: public scoring_function
{
	weighted_terms(const terms* t, const flv& weights); // does not own t
	virtual ~weighted_terms()
	{
	}
	fl cutoff() const
	{
		return cutoff_;
	}
	result_components eval_fast(smt t1, smt t2, fl r) const; // intentionally not checking for cutoff
	fl eval_slow(const atom_base& a, const atom_base& b, fl r) const; //dkoes - da terms

	fl conf_independent(const model& m, fl e) const;

	//dkoes - return true if has slow terms that can't be precalculated
	bool has_slow() const
	{
		return enabled_distance_additive_terms.size() > 0
				|| enabled_distance_additive_terms.size() > 0; //this aren't working yet
	}

	//returns index after last used component
	sz num_used_components() const
	{
		if (enabled_charge_dependent_terms.size() > 0)
			return result_components::size();
		else
			return 1;
	}

	const terms* unweighted_terms() const
	{
		return t;
	}

	fl weight(unsigned i) const {
		return weights[i];
	}

	sz size() const {
		return weights.size();
	}
private:
	weighted_terms() :
			t(NULL), cutoff_(0), conf_indep_start(0)
	{
	}
	const terms* t;
	flv weights;
	fl cutoff_;
	sz conf_indep_start; //dkoes - where conf independent terms are in weights
	szv enabled_charge_independent_terms;
	szv enabled_charge_dependent_terms;
	szv enabled_distance_additive_terms; //additive currently aren't supported
};

#endif
