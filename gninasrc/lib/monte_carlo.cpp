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

#include "monte_carlo.h"
#include "coords.h"
#include "mutate.h"
#include "quasi_newton.h"

output_type monte_carlo::operator()(model& m, const precalculate& p, igrid& ig,
    const vec& corner1, const vec& corner2, incrementable* increment_me,
    rng& generator, grid& user_grid, igrid& ig_metropolis) const {
  output_container tmp;
  this->operator()(m, tmp, p, ig, corner1, corner2, increment_me, generator,
      user_grid,ig_metropolis); // call the version that produces the whole container
  VINA_CHECK(!tmp.empty());
  return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng& generator) {
  if (new_f < old_f) return true;
  const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
  return random_fl(0, 1, generator) < acceptance_probability;
}

void update_energy(model& m, output_type& out, const vec& v, igrid* ig){
  ig->adjust_center(m);
  out.e = ig->eval(m, v[1]);
}

void monte_carlo::single_run(model& m, output_type& out, const precalculate& p,
    igrid& ig, rng& generator, grid& user_grid, igrid& ig_metropolis) const {
  conf_size s = m.get_size();
  change g(s, ig.move_receptor());
  vec authentic_v(1000, 1000, 1000);
  out.e = max_fl;
  output_type current(out);
  minimization_params minparms = ssd_par.minparm;
  if (minparms.maxiters == 0) minparms.maxiters = ssd_par.evals;

  quasi_newton quasi_newton_par(minparms);
  VINA_U_FOR(step, num_steps) {
    output_type candidate(current.c, max_fl);
    mutate_conf(candidate.c, m, mutation_amplitude, generator);
    quasi_newton_par(m, p, ig, candidate, g, hunt_cap, user_grid);
    update_energy(m, candidate, authentic_v, &ig_metropolis);
    if (step == 0
        || metropolis_accept(current.e, candidate.e, temperature, generator)) {
      quasi_newton_par(m, p, ig, candidate, g, authentic_v, user_grid);
      current = candidate;
      update_energy(m, current, authentic_v, &ig_metropolis);
      if (current.e < out.e) out = current;
    }
  }
  quasi_newton_par(m, p, ig, out, g, authentic_v, user_grid);
}

void monte_carlo::many_runs(model& m, output_container& out,
    const precalculate& p, igrid& ig, const vec& corner1, const vec& corner2,
    sz num_runs, rng& generator, grid& user_grid,igrid& ig_metropolis) const {
  conf c(m.get_size(), ig.move_receptor());
  VINA_FOR(run, num_runs) {
    output_type tmp(c, 0);
    tmp.c.randomize(corner1, corner2, generator);
    single_run(m, tmp, p, ig, generator, user_grid,ig_metropolis);
    out.push_back(new output_type(tmp));
  }
  out.sort();
}

output_type monte_carlo::many_runs(model& m, const precalculate& p, igrid& ig,
    const vec& corner1, const vec& corner2, sz num_runs, rng& generator,
    grid& user_grid,igrid& ig_metropolis) const {
  output_container tmp;
  many_runs(m, tmp, p, ig, corner1, corner2, num_runs, generator, user_grid,ig_metropolis);
  VINA_CHECK(!tmp.empty());
  return tmp.front();
}

// out is sorted
void monte_carlo::operator()(model& m, output_container& out,
    const precalculate& p, igrid& ig, const vec& corner1, const vec& corner2,
    incrementable* increment_me, rng& generator, grid& user_grid,igrid& ig_metropolis) const {
  vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
  conf_size s = m.get_size();
  change g(s, ig.move_receptor());
  output_type tmp(conf(s, ig.move_receptor()), 0);
  tmp.c.randomize(corner1, corner2, generator);
  fl best_e = max_fl;
  minimization_params minparms = ssd_par.minparm;

  if (minparms.maxiters == 0) minparms.maxiters = ssd_par.evals;
  quasi_newton quasi_newton_par(minparms);
  VINA_U_FOR(step, num_steps) {
    if (increment_me) ++(*increment_me);
    output_type candidate = tmp;
    mutate_conf(candidate.c, m, mutation_amplitude, generator);

    if (minparms.single_min) //use full v to begin with
      quasi_newton_par(m, p, ig, candidate, g, authentic_v, user_grid);
    else
      quasi_newton_par(m, p, ig, candidate, g, hunt_cap, user_grid);
    
    update_energy(m, candidate, authentic_v, &ig_metropolis);
    
    if (step == 0
        || metropolis_accept(tmp.e, candidate.e, temperature, generator)) {
      tmp = candidate;

      m.set(tmp.c); // FIXME? useless?

      // FIXME only for very promising ones
      if (tmp.e < best_e || out.size() < num_saved_mins) {

        if (!minparms.single_min) { //refine with full v
          quasi_newton_par(m, p, ig, tmp, g, authentic_v, user_grid);
          update_energy(m, tmp, authentic_v, &ig_metropolis);
          m.set(tmp.c); // FIXME? useless?
        }
        tmp.coords = m.get_heavy_atom_movable_coords();
        add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
        if (tmp.e < best_e) {
          best_e = tmp.e;
        }
      }
    }
  }
  VINA_CHECK(!out.empty());
  VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
