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

#ifndef VINA_NON_CACHE_CNN_H
#define VINA_NON_CACHE_CNN_H

#include "igrid.h"
#include "szv_grid.h"
#include "non_cache.h"
#include "dl_scorer.h"

struct parallel_mc_aux;
struct non_cache_cnn : public non_cache {
    non_cache_cnn(szv_grid_cache& gcache, const grid_dims& gd_,
        const precalculate* p_, fl slope_, DLScorer& dl_scorer_);
    virtual ~non_cache_cnn() {
    }
    virtual fl eval(model& m, fl v) const; // needs m.coords
    virtual fl eval_deriv(model& m, fl v, const grid& user_grid) const; // needs m.coords, sets m.minus_forces
    virtual bool within(const model& m, fl margin = 0.0001) const;
    void setSlope(fl sl) {
      slope = sl;
    }
    fl getSlope() const {
      return slope;
    }
    virtual bool skip_interacting_pairs() const {
      return true;
    }
    virtual void adjust_center(model& m);
    virtual vec get_center() const;
    const DLScorer& get_scorer() const {
      return dl_scorer;
    }
    virtual bool move_receptor() {
      return false; //stopped supporting this with torch
    }
  protected:
    DLScorer& dl_scorer;
    grid_dims cnn_gd;
};

#endif

