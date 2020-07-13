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

#ifndef VINA_NON_CACHE_H
#define VINA_NON_CACHE_H

#include "igrid.h"
#include "szv_grid.h"

struct non_cache : public igrid {
    non_cache(szv_grid_cache& gcache, const grid_dims& gd_,
        const precalculate* p_, fl slope_ = 1e3);
    virtual ~non_cache() {
    }
    virtual fl eval(const model& m, fl v) const; // needs m.coords // clean up
    virtual fl eval_deriv(model& m, fl v, const grid& user_grid) const; // needs m.coords, sets m.minus_forces // clean up

    fl check_bounds(const grid_dims& dims, const vec& a_coords,
        vec& adjusted_a_coords) const;
    fl check_bounds_deriv(const grid_dims& dims, const vec& a_coords,
        vec& adjusted_a_coords, vec& out_of_bounds_deriv) const;

    virtual bool within(const model& m, fl margin = 0.0001) const;
    fl getSlope() {
      return slope;
    }
    virtual void setSlope(fl sl) {
      slope = sl;
    }
    virtual const precalculate* get_precalculate() const {
      return p;
    }
    grid_dims get_grid_dims() const {
      return gd;
    }

  protected:
    fl slope;
    szv_grid sgrid;
    grid_dims gd;
    const precalculate* p;

    bool gd_within(const grid_dims& dims, const model& m, fl margin) const;
};

#endif

