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

#ifndef VINA_CACHE_H
#define VINA_CACHE_H

#include <string>
#include "igrid.h"
#include "grid.h"
#include "model.h"
#include "array3d.h"

struct cache_mismatch {
};
struct rigid_mismatch : public cache_mismatch {
};
struct grid_dims_mismatch : public cache_mismatch {
};
struct energy_mismatch : public cache_mismatch {
};

struct cache : public igrid {
    cache(const std::string& scoring_function_version_, const grid_dims& gd_,
        fl slope_);
    fl eval(model& m, fl v) const; // needs m.coords // clean up
    fl eval_deriv(model& m, fl v, const grid& user_grid) const; // needs m.coords, sets m.minus_forces // clean up

    virtual void populate(const model& m, const precalculate& p,
        const std::vector<smt>& atom_types_needed, grid& user_grid,
        bool display_progress = true);
    virtual ~cache() {
    }
    ;
  private:
    std::string scoring_function_version;
    atomv atoms; // for verification
    grid_dims gd;
    fl slope; // does not get (de-)serialized
    std::vector<grid> grids;
    friend class boost::serialization::access;
    friend class cache_gpu;
    template<class Archive>
    void save(Archive& ar, const unsigned version) const;
    template<class Archive>
    void load(Archive& ar, const unsigned version);
        BOOST_SERIALIZATION_SPLIT_MEMBER()};

#endif
