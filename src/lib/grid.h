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

#ifndef VINA_GRID_H
#define VINA_GRID_H

#include "array3d.h"
#include "grid_dim.h"
#include "curl.h"
#include "result_components.h"
#include "atom.h"

class grid
{ // FIXME rm 'm_', consistent with my new style
	vec m_init;
	vec m_range;
	vec m_factor;
	vec m_dim_fl_minus_1;
	vec m_factor_inv;
	array3d<fl> data;
	array3d<fl> chargedata; //needs to be multiplied by atom charge

	friend class cache;
	friend class non_cache;
	public:
	grid() :
			m_init(0, 0, 0), m_range(1, 1, 1), m_factor(1, 1, 1), m_dim_fl_minus_1(
					-1, -1, -1), m_factor_inv(1, 1, 1)
	{
	} // not private
	grid(const grid_dims& gd, bool hascharged)
	{
		init(gd, hascharged);
	}
	void init(const grid_dims& gd, bool hascharged);
    void init(const grid_dims& gd, std::istream& user_in, fl ug_scaling_factor);
	vec index_to_argument(sz x, sz y, sz z) const
	{
		return vec(m_init[0] + m_factor_inv[0] * x,
				m_init[1] + m_factor_inv[1] * y,
				m_init[2] + m_factor_inv[2] * z);
	}
	bool initialized() const
	{
		return data.dim0() > 0 && data.dim1() > 0 && data.dim2() > 0;
	}
	fl evaluate(const atom& a, const vec& location, fl slope, fl c, vec* deriv = NULL) const;
    fl evaluate_user(const vec& location, fl slope, vec* deriv = NULL) const;
private:
	fl evaluate_aux(const array3d<fl>& m_data, const vec& location, fl slope,
			fl v, vec* deriv) const; // sets *deriv if not NULL
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version)
	{
		ar & m_init;
		ar & data;
		ar & chargedata;
		ar & m_range;
		ar & m_factor;
		ar & m_dim_fl_minus_1;
		ar & m_factor_inv;
	}
};

#endif
