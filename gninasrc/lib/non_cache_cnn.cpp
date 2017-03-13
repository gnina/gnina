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

non_cache_cnn::non_cache_cnn(szv_grid_cache& gcache,
			     const grid_dims& gd_,
			     const precalculate* p_,
			     fl slope_,
			     CNNScorer& cnn_scorer_) :
	non_cache(gcache, gd_, p_, slope_), cnn_scorer(cnn_scorer_)
{
}

fl non_cache_cnn::eval(model& m, fl v) const
{
	fl e = 0;
	sz n = num_atom_types();
	VINA_FOR(i, m.num_movable_atoms())
	{
		const atom_base& a = m.movable_atom(i);
		smt t1 = a.get();
		if (t1 >= n || is_hydrogen(t1))
			continue;

		const vec& a_coords = m.movable_coords(i);
		vec adjusted_a_coords;
		fl out_of_bounds_penalty = check_bounds(a_coords, adjusted_a_coords);

		e += out_of_bounds_penalty;
	}
	fl aff = 0;
	e += -cnn_scorer.score(m, false, aff);
	return e;
}

fl non_cache_cnn::eval_deriv(model& m, fl v, const grid& user_grid) const
{
	fl e = 0;
	sz n = num_atom_types();
	VINA_FOR(i, m.num_movable_atoms())
	{
		const atom_base& a = m.movable_atom(i);
		smt t1 = a.get();
		if (t1 >= n || is_hydrogen(t1))
		{
			m.movable_minus_forces(i).assign(0);
			continue;
		}

		const vec& a_coords = m.movable_coords(i);
		vec adjusted_a_coords;
		vec out_of_bounds_deriv(0, 0, 0);
		fl out_of_bounds_penalty = check_bounds_deriv(a_coords, adjusted_a_coords, out_of_bounds_deriv);

		fl this_e = 0;
		vec deriv(0, 0, 0);
		//VINA_FOR_IN(...) { per-atom cnn score would be here }
		if (user_grid.initialized())
		{
			vec ug_deriv(0, 0, 0);
			fl uge = user_grid.evaluate_user(a_coords, slope, &ug_deriv);
			this_e += uge;
			deriv += ug_deriv;
		}
		curl(this_e, deriv, v);
		m.movable_minus_forces(i) = deriv + out_of_bounds_deriv;
		e += this_e + out_of_bounds_penalty;
	}
	fl aff = 0;
	e = -cnn_scorer.score(m, true, aff);
	return e;
}

