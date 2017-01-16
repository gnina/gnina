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

#include <algorithm> // fill, etc
#if 0 // use binary cache
// for some reason, binary archive gives four huge warnings in VC2008
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
typedef boost::archive::binary_iarchive iarchive;
typedef boost::archive::binary_oarchive oarchive;
#else // use text cache
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
typedef boost::archive::text_iarchive iarchive;
typedef boost::archive::text_oarchive oarchive;
#endif 

#include <boost/serialization/split_member.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/static_assert.hpp>
#include "cache.h"
#include "file.h"
#include "szv_grid.h"

cache::cache(const std::string& scoring_function_version_, const grid_dims& gd_,
		fl slope_) :
		scoring_function_version(scoring_function_version_), gd(gd_), slope(
				slope_), grids(num_atom_types())
{
}

fl cache::eval(const model& m, fl v) const
		{ // needs m.coords
	fl e = 0;
	sz nat = num_atom_types();

	VINA_FOR(i, m.num_movable_atoms())
	{
		const atom& a = m.atoms[i];
		smt t = a.get();
		if (t >= nat || is_hydrogen(t))
			continue;
		const grid& g = grids[t];
		assert(g.initialized());
		e += g.evaluate(a, m.coords[i], slope, v);
	}
	return e;
}

fl cache::eval_deriv(model& m, fl v, const grid& user_grid) const
		{ // needs m.coords, sets m.minus_forces
	fl e = 0;
	sz nat = num_atom_types();

	VINA_FOR(i, m.num_movable_atoms())
	{
		const atom& a = m.atoms[i];
		smt t = a.get();
		if (t >= nat || is_hydrogen(t))
		{
			m.minus_forces[i].assign(0);
			continue;
		}
		const grid& g = grids[t];
		assert(g.initialized());
		vec deriv;
		e += g.evaluate(a, m.coords[i], slope, v, &deriv);
		m.minus_forces[i] = deriv;
	}
	return e;
}


template<class Archive>
void cache::save(Archive& ar, const unsigned version) const
		{
	ar & scoring_function_version;
	ar & gd;
	ar & grids;
}

template<class Archive>
void cache::load(Archive& ar, const unsigned version)
{
	std::string name_tmp;
	ar & name_tmp;
	if (name_tmp != scoring_function_version)
		throw energy_mismatch();
	grid_dims gd_tmp;
	ar & gd_tmp;
	if (!eq(gd_tmp, gd))
		throw grid_dims_mismatch();

	ar & grids;
}

void cache::populate(const model& m, const precalculate& p,
		const std::vector<smt>& atom_types_needed, grid& user_grid, bool display_progress)
{
	std::vector<smt> needed;
	bool haschargeterms = p.has_components();

	VINA_FOR_IN(i, atom_types_needed)
	{
		smt t = atom_types_needed[i];
		if (!grids[t].initialized())
		{
			needed.push_back(t);
			grids[t].init(gd, haschargeterms);
		}
	}
	if (needed.empty())
		return;
	flv affinities(needed.size());
	flv chargeaffinities;
	if(haschargeterms)
		chargeaffinities.resize(needed.size());
           
	sz nat = num_atom_types();

	grid& g = grids[needed.front()];

	const fl cutoff_sqr = p.cutoff_sqr();

	szv_grid_cache igcache(m, cutoff_sqr);
	szv_grid ig(igcache, gd);

	VINA_FOR(x, g.data.dim0())
	{
		VINA_FOR(y, g.data.dim1())
		{
			VINA_FOR(z, g.data.dim2())
			{
				std::fill(affinities.begin(), affinities.end(), 0);
				std::fill(chargeaffinities.begin(), chargeaffinities.end(), 0);
				vec probe_coords;
				probe_coords = g.index_to_argument(x, y, z);
				const szv& possibilities = ig.possibilities(probe_coords);
				VINA_FOR_IN(possibilities_i, possibilities)
				{
					const sz i = possibilities[possibilities_i];
					const atom& a = m.grid_atoms[i];
					const smt t1 = a.get();
					const fl r2 = vec_distance_sqr(a.coords, probe_coords);
					if (r2 <= cutoff_sqr)
					{
						VINA_FOR_IN(j, needed)
						{
							const smt t2 = needed[j];
							assert(t2 < nat);
							//t1 is the receptor atom, a
							//t2 is type from the ligand, not corresponding to any
							//particular atom
							result_components val = p.eval_fast(t1, t2, r2);
							if (haschargeterms)
							{
								//affinities contains the terms that are independent of
								//the ligand atom charge

								affinities[j] +=
										val[result_components::TypeDependentOnly]
												+
												val[result_components::AbsAChargeDependent]
														* fabs(a.charge);
								//this component must be multiplied by the ligand atom charge
								chargeaffinities[j] +=
										val[result_components::AbsBChargeDependent] +
										val[result_components::ABChargeDependent]*a.charge; //not abs value
							}
							else
							{
								affinities[j] +=
										val[result_components::TypeDependentOnly];
							}
						}
					}
				}
				VINA_FOR_IN(j, needed)
				{
					sz t = needed[j];
					assert(t < nat);
					grids[t].data(x, y, z) = affinities[j]; //+ user_grid.evaluate_user(vec(x, y, z));
					if(haschargeterms)
						grids[t].chargedata(x, y, z) = chargeaffinities[j];
                    if(user_grid.initialized())
                        grids[t].data(x, y, z) += user_grid.evaluate_user(vec(x, y ,z), slope);
				}
			}
		}
	}
}
