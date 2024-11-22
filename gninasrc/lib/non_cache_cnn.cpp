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

#include "non_cache_cnn.h"
#include "curl.h"
#include "loop_timer.h"

non_cache_cnn::non_cache_cnn(szv_grid_cache& gcache, const grid_dims& gd_,
    const precalculate* p_, fl slope_, DLScorer& dl_scorer_)
    : non_cache(gcache, gd_, p_, slope_), dl_scorer(dl_scorer_) {
}

//return the LOSS plus any out of bounds penalties
fl non_cache_cnn::eval(model& m, fl v) const {
  fl e = 0;
  sz n = num_atom_types();
  VINA_FOR(i, m.num_movable_atoms()) {
    const atom_base& a = m.movable_atom(i);
    smt t1 = a.get();
    if (t1 >= n || is_hydrogen(t1)) continue;

    const vec& a_coords = m.movable_coords(i);
    vec adjusted_a_coords;
    fl out_of_bounds_penalty = check_bounds(gd, a_coords, adjusted_a_coords);
    out_of_bounds_penalty += check_bounds(cnn_gd, a_coords, adjusted_a_coords);
    e += out_of_bounds_penalty;
  }
  fl aff = 0;
  fl loss = 0;
  fl variance = 0;
  dl_scorer.score(m, false, aff, loss, variance);
  e += loss;

  return e;
}

//reset center; will apply inverse of receptor transformation to ligand as well
void non_cache_cnn::adjust_center(model& m) {
  //the cnn is only defined over a cubic region, set
  //a second out_of_bound_box to this region

  //set center
  dl_scorer.set_center_from_model(m);

  //apply oob penalty to cnn grid
  dl_scorer.set_bounding_box(cnn_gd);

}

vec non_cache_cnn::get_center() const {
  return dl_scorer.get_center();
}

//check cnn box
bool non_cache_cnn::within(const model& m, fl margin) const {
  return gd_within(cnn_gd, m, margin) || non_cache::within(m, margin);
}

//return the cnn loss plus any out of bounds penalties
fl non_cache_cnn::eval_deriv(model& m, fl v, const grid& user_grid) const {
  fl e = 0;
  sz n = num_atom_types();
  fl aff = 0;
  fl loss = 0;
  fl variance = 0;
  const fl cutoff_sqr = p->cutoff_sqr();
  //this is what compute cnn minus_forces
  dl_scorer.score(m, true, aff, loss, variance);
  e += loss;

  //out of bonds forces
  VINA_FOR(i, m.num_movable_atoms()) {
    const atom_base& a = m.movable_atom(i);
    smt t1 = a.get();
    if (t1 >= n || is_hydrogen(t1)) {
      m.movable_minus_forces(i).assign(0);
      continue;
    }
  
    const vec& a_coords = m.movable_coords(i);
    vec adjusted_a_coords;
    //seperate adjusted coords for empirical bounds
    vec adjusted_emp_a_coords;
    vec out_of_bounds_deriv(0, 0, 0);
    vec cnn_out_of_bounds_deriv(0, 0, 0);
    fl out_of_bounds_penalty = check_bounds_deriv(gd, a_coords,
        adjusted_emp_a_coords, out_of_bounds_deriv);
    out_of_bounds_penalty += check_bounds_deriv(cnn_gd, a_coords,
        adjusted_a_coords, cnn_out_of_bounds_deriv);
    fl this_e = 0;
    vec deriv(0, 0, 0);
    //empirical e and deriv values
    fl emp_e =0;
    vec emp_deriv(0, 0, 0);
    
    if (dl_scorer.options().mix_emp_force){ 
    	
    	const szv& possibilities = sgrid.possibilities(adjusted_emp_a_coords);
        VINA_FOR_IN(possibilities_j, possibilities) {
            const sz j = possibilities[possibilities_j];
            const atom& b = m.grid_atoms[j];
            smt t2 = b.get();
            
            vec r_ba;
            r_ba = adjusted_emp_a_coords - b.coords;
            fl r2 = sqr(r_ba);
            if (r2 < cutoff_sqr) {
                if (r2 < epsilon_fl) {
                    throw std::runtime_error(
                            "Ligand atom exactly overlaps receptor atom.  I can't deal with this.");
                }
                //dkoes - the "derivative" value returned by eval_deriv
                //is normalized by r (dor = derivative over r?)
                pr e_dor = p->eval_deriv(a, b, r2);
                emp_e += e_dor.first;
                emp_deriv += e_dor.second * r_ba;
                vec out_deriv = e_dor.second * r_ba;
            }
        }
    }
    //VINA_FOR_IN(...) { per-atom cnn score would be here }
    if (user_grid.initialized()) {
      vec ug_deriv(0, 0, 0);
      fl uge = user_grid.evaluate_user(a_coords, slope, &ug_deriv);
      this_e += uge;
      deriv += ug_deriv;
      if (dl_scorer.options().mix_emp_force){
          emp_e += uge;
          emp_deriv += ug_deriv;
      }
    }
    curl(this_e, deriv, v);
    m.movable_minus_forces(i) += deriv + out_of_bounds_deriv
        + cnn_out_of_bounds_deriv;
    if (dl_scorer.options().mix_emp_force){
        curl(emp_e,emp_deriv,v);
        m.movable_minus_forces(i) += dl_scorer.options().empirical_weight * (emp_deriv + out_of_bounds_deriv);
        //rescale
        m.movable_minus_forces(i) /= (1.0 + dl_scorer.options().empirical_weight);
    }
    e += this_e + out_of_bounds_penalty;
     if (dl_scorer.options().mix_emp_energy){
      e += dl_scorer.options().empirical_weight * emp_e;
     }
    }
 if (dl_scorer.options().mix_emp_energy){
     e /=  (1.0 + dl_scorer.options().empirical_weight);
	}
  return e;
}

