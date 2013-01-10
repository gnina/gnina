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

#include "szv_grid.h"
#include "brick.h"

szv_grid::szv_grid(const model& m_, const grid_dims& gd, fl cutoff_sqr) : m(m_), cutoff_sq(cutoff_sqr), m_data(gd[0].n, gd[1].n, gd[2].n) {
	vec end;
	VINA_FOR_IN(i, gd) {
		m_init[i] = gd[i].begin;
		end   [i] = gd[i].end;
	}
	m_range = end - m_init;

	const sz nat = num_atom_types();

	VINA_FOR_IN(i, m.grid_atoms) {
		const atom& a = m.grid_atoms[i];
		if(a.get() < nat && !a.is_hydrogen() && brick_distance_sqr(m_init, end, a.coords) < cutoff_sq)
			relevant_indexes.push_back(i);
	}
	//don't precompute - this is particularly inefficient for minimization
}

//get the receptor atoms that are probably within the cutoff at index i,j,k
//the result will be computed and memoized on demand
const szv& szv_grid::get(sz x, sz y, sz z) const
{
	szv& possibles = m_data(x,y,z);
	if(possibles.size() == 0) //not yet computed
	{
		VINA_FOR_IN(ri, relevant_indexes) {
			const sz i = relevant_indexes[ri];
			const atom& a = m.grid_atoms[i];
			if(brick_distance_sqr(index_to_coord(x, y, z), index_to_coord(x+1, y+1, z+1), a.coords) < cutoff_sq)
				possibles.push_back(i);
		}
		if(possibles.size() == 0) //mark that it is suppose to be empty
			possibles.push_back(UINT_MAX);
	}
	if(possibles.size() == 1 && possibles[0] == UINT_MAX)
	{
		//already computed, but nothing was there
		return empty;
	}
	return possibles;
}



const szv& szv_grid::possibilities(const vec& coords) const {
	boost::array<sz, 3> index;
	VINA_FOR_IN(i, index) {
		assert(coords[i] + epsilon_fl >= m_init[i]);
		assert(coords[i] <= m_init[i] + m_range[i] + epsilon_fl);
		const fl tmp = (coords[i] - m_init[i]) * m_data.dim(i) / m_range[i];
		index[i] = fl_to_sz(tmp, m_data.dim(i) - 1);
	}
	return get(index[0], index[1], index[2]);
}

vec szv_grid::index_to_coord(sz i, sz j, sz k) const {
	vec index(i, j, k);
	vec tmp;
	VINA_FOR_IN(n, tmp) 
		tmp[n] = m_init[n] + m_range[n] * index[n] / m_data.dim(n);
	return tmp;
}

grid_dims szv_grid_dims(const grid_dims& gd) {
	grid_dims tmp;
	VINA_FOR_IN(i, tmp) {
		tmp[i].begin = gd[i].begin;
		tmp[i].end   = gd[i].end;
		fl n_fl = (gd[i].end - gd[i].begin) / 3; // 3A preferred size
		int n_int = int(n_fl);
		tmp[i].n     = (n_int < 1) ?  1 : sz(n_int);
	}
	return tmp;
}
