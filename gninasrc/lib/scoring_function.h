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

#ifndef VINA_SCORING_FUNCTION_H
#define VINA_SCORING_FUNCTION_H

#include "atom_type.h"
#include "result_components.h"

struct model; // forward declaration

struct scoring_function {
	virtual fl cutoff() const = 0;
	virtual result_components eval_fast(smt t1, smt t2, fl r) const = 0;
	//dkoes - split eval between precalculatable scoring and not
	virtual bool has_slow() const = 0;
	virtual sz num_used_components() const = 0;
	virtual fl eval_slow(const atom_base& a, const atom_base& b, fl r) const = 0;
	virtual fl conf_independent(const model& m, fl e) const = 0;
	virtual ~scoring_function() {}
};

#endif
