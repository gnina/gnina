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

#ifndef VINA_SZV_GRID_H
#define VINA_SZV_GRID_H

#include "model.h"
#include "grid_dim.h"
#include "array3d.h"

//dkoes - this keeps track of what receptor atoms are possibly close enough
//to grid points to matter and caches their indices
struct szv_grid {
	szv_grid(const model& m_, const grid_dims& gd, fl cutoff_sqr);
	const szv& possibilities(const vec& coords) const;
private:
	const model& m;
	fl cutoff_sq;
	szv relevant_indexes; //rec atoms within distance of docking grid
	mutable array3d<szv> m_data; //this is updated as needed
	vec m_init;
	vec m_range;
	szv empty; //referenced when there are no possiblities
	const szv& get(sz i, sz j, sz k) const;
	vec index_to_coord(sz i, sz j, sz k) const;
};

grid_dims szv_grid_dims(const grid_dims& gd);


#endif
