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

#ifndef VINA_IGRID_H
#define VINA_IGRID_H

#include "common.h"
#include "grid.h"

struct model;
// forward declaration

struct igrid { // grids interface (that cache, etc. conform to)
    virtual fl eval(model& m, fl v) const = 0; // needs m.coords // clean up
    virtual fl eval_deriv(model& m, fl v, const grid& user_grid) const = 0; // needs m.coords, sets m.minus_forces // clean up
    virtual bool skip_interacting_pairs() const {
      return false;
    } //if true, evaluates the entire model, not just PL interactions
    virtual void adjust_center(model& m) {
    } //for cnn
    virtual vec get_center() const {
      return vec(0, 0, 0);
    } //current center
    virtual bool move_receptor() {
      return false;
    } //for cnn, if we are moving receptor
};

#endif
