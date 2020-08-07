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

#include "non_cache.h"
#include "curl.h"
#include "loop_timer.h"

non_cache::non_cache(szv_grid_cache& gcache, const grid_dims& gd_,
    const precalculate* p_, fl slope_)
    : sgrid(gcache, gd_), gd(gd_), p(p_), slope(slope_) {
}

fl non_cache::check_bounds(const grid_dims& dims, const vec& a_coords,
    vec& adjusted_a_coords) const {
  fl out_of_bounds_penalty = 0;
  adjusted_a_coords = a_coords;
  VINA_FOR_IN(j, gd) {
    if (dims[j].n > 0) {
      if (a_coords[j] < dims[j].begin) {
        adjusted_a_coords[j] = dims[j].begin;
        out_of_bounds_penalty += std::abs(a_coords[j] - dims[j].begin);
      } else
        if (a_coords[j] > dims[j].end) {
          adjusted_a_coords[j] = dims[j].end;
          out_of_bounds_penalty += std::abs(a_coords[j] - dims[j].end);
        }
    }
  }
  out_of_bounds_penalty *= slope;
  return out_of_bounds_penalty;
}

fl non_cache::eval(const model& m, fl v) const { // clean up
  fl e = 0;
  const fl cutoff_sqr = p->cutoff_sqr();
  sz n = num_atom_types();
  VINA_FOR(i, m.num_movable_atoms()) {
    const atom& a = m.atoms[i];
    smt t1 = a.get();
    if (t1 >= n || is_hydrogen(t1)) continue;

    const vec& a_coords = m.coords[i];
    vec adjusted_a_coords;
    fl out_of_bounds_penalty = check_bounds(gd, a_coords, adjusted_a_coords);
    fl this_e = 0;
    const szv& possibilities = sgrid.possibilities(adjusted_a_coords);
    VINA_FOR_IN(possibilities_j, possibilities) {
      const sz j = possibilities[possibilities_j];
      const atom& b = m.grid_atoms[j];
      smt t2 = b.get();
      vec r_ba;
      r_ba = adjusted_a_coords - b.coords; // FIXME why b-a and not a-b ?
      fl r2 = sqr(r_ba);
      if (r2 < cutoff_sqr) {
        //jac241 - Use adjusted_a_coords or just a_coords?
        //also how to verify they're ligand coordinates (table lookup?)
        this_e += p->eval(a, b, r2); // + user_grid.evaluate_user(adjusted_a_coords, slope, NULL);
      }
    }
    curl(this_e, v);
    e += this_e + out_of_bounds_penalty;
  }
  return e;
}

bool non_cache::within(const model& m, fl margin) const {
  return gd_within(gd, m, margin);
}

bool non_cache::gd_within(const grid_dims& dims, const model& m,
    fl margin) const {
  VINA_FOR(i, m.num_movable_atoms()) {
    if (m.atoms[i].is_hydrogen()) continue;
    const vec& a_coords = m.coords[i];
    VINA_FOR_IN(j, dims)
      if (dims[j].n > 0)
        if (a_coords[j] < dims[j].begin - margin
            || a_coords[j] > dims[j].end + margin) return false;
  }
  return true;
}

fl non_cache::check_bounds_deriv(const grid_dims& dims, const vec& a_coords,
    vec& adjusted_a_coords, vec& out_of_bounds_deriv) const {
  fl out_of_bounds_penalty = 0;
  adjusted_a_coords = a_coords;
  VINA_FOR_IN(j, gd) {
    if (dims[j].n > 0) {
      if (a_coords[j] < dims[j].begin) {
        adjusted_a_coords[j] = dims[j].begin;
        out_of_bounds_deriv[j] = -1;
        out_of_bounds_penalty += std::abs(a_coords[j] - dims[j].begin);
      } else
        if (a_coords[j] > dims[j].end) {
          adjusted_a_coords[j] = dims[j].end;
          out_of_bounds_deriv[j] = 1;
          out_of_bounds_penalty += std::abs(a_coords[j] - dims[j].end);
        }
    }
  }
  out_of_bounds_penalty *= slope;
  out_of_bounds_deriv *= slope;
  return out_of_bounds_penalty;
}

fl non_cache::eval_deriv(model& m, fl v, const grid& user_grid) const { // clean up
  fl e = 0;
  const fl cutoff_sqr = p->cutoff_sqr();

  sz n = num_atom_types();

  VINA_FOR(i, m.num_movable_atoms()) {
    const atom& a = m.atoms[i];
    smt t1 = a.get();
    if (t1 >= n || is_hydrogen(t1)) {
      m.minus_forces[i].assign(0);
      continue;
    }

    const vec& a_coords = m.coords[i];
    vec adjusted_a_coords;
    vec out_of_bounds_deriv(0, 0, 0);
    fl out_of_bounds_penalty = check_bounds_deriv(gd, a_coords,
        adjusted_a_coords, out_of_bounds_deriv);

    fl this_e = 0;
    vec deriv(0, 0, 0);
    const szv& possibilities = sgrid.possibilities(adjusted_a_coords);
    VINA_FOR_IN(possibilities_j, possibilities) {
      const sz j = possibilities[possibilities_j];
      const atom& b = m.grid_atoms[j];
      smt t2 = b.get();
      vec r_ba;
      r_ba = adjusted_a_coords - b.coords;
      fl r2 = sqr(r_ba);
      if (r2 < cutoff_sqr) {
        if (r2 < epsilon_fl) {
          throw std::runtime_error(
              "Ligand atom exactly overlaps receptor atom.  I can't deal with this.");
        }
        //dkoes - the "derivative" value returned by eval_deriv
        //is normalized by r (dor = derivative over r?)
        pr e_dor = p->eval_deriv(a, b, r2);
        this_e += e_dor.first;
        deriv += e_dor.second * r_ba;
        vec out_deriv = e_dor.second * r_ba;
      }
    }
    if (user_grid.initialized()) {
      vec ug_deriv(0, 0, 0);
      fl uge = user_grid.evaluate_user(a_coords, slope, &ug_deriv);
      this_e += uge;
      deriv += ug_deriv;
    }
    curl(this_e, deriv, v);
    m.minus_forces[i] = deriv + out_of_bounds_deriv;
    e += this_e + out_of_bounds_penalty;
  }
  return e;
}

