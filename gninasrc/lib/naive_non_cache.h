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

#ifndef VINA_NAIVE_NON_CACHE_H
#define VINA_NAIVE_NON_CACHE_H

#include "igrid.h"
#include "model.h"

struct naive_non_cache : public igrid {
    naive_non_cache(const precalculate* p_);
    virtual fl eval(model& m, fl v) const; // needs m.coords
    virtual fl eval_deriv(model& m, fl v, const grid& user_grid) const {
      VINA_CHECK(false);
      return 0;
    } // unused'
    //virtual fl eval_usergrid(const model& m, grid& user_grid);
    void eval_atoms(const model& m, std::vector<flv>& per_atom_values) const;
  private:
    const precalculate* p;
};

#endif
