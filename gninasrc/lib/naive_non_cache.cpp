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

#include "naive_non_cache.h"
#include "curl.h"

naive_non_cache::naive_non_cache(const precalculate* p_)
    : p(p_) {
}
fl naive_non_cache::eval(model& m, fl v) const { // needs m.coords
  fl e = 0;
  const fl cutoff_sqr = p->cutoff_sqr();

  sz n = num_atom_types();

  VINA_FOR(i, m.num_movable_atoms()) {
    fl this_e = 0;
    const atom& a = m.atoms[i];
    smt t1 = a.get();
    if (t1 >= n || is_hydrogen(t1)) continue;
    const vec& a_coords = m.coords[i];

    VINA_FOR_IN(j, m.grid_atoms) {
      const atom& b = m.grid_atoms[j];
      smt t2 = b.get();
      if (t2 >= n || is_hydrogen(t2)) continue;
      vec r_ba;
      r_ba = a_coords - b.coords;
      fl r2 = sqr(r_ba);
      if (r2 < cutoff_sqr) {
        this_e += p->eval(a, b, r2);
      }
    }
    curl(this_e, v);
    e += this_e;
  }
  return e;
}

