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
	fl score = cnn_scorer.score(m, false);
	return -score;
}

bool non_cache_cnn::within(const model& m, fl margin) const
{
	VINA_FOR(i, m.num_movable_atoms())
	{
		if (m.movable_atom(i).is_hydrogen())
			continue;
		const vec& a_coords = m.movable_coords(i);
		VINA_FOR_IN(j, gd)
			if (gd[j].n > 0)
				if (a_coords[j] < gd[j].begin - margin
						|| a_coords[j] > gd[j].end + margin)
					return false;
	}
	return true;
}

fl non_cache_cnn::eval_deriv(model& m, fl v, const grid& user_grid) const
{
	fl score = cnn_scorer.score(m, true);
	return -score;
}

