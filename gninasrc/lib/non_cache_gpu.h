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

#ifndef VINA_NON_CACHE_GPU_H
#define VINA_NON_CACHE_GPU_H

#include "non_cache.h"
#include "precalculate_gpu.h"

struct non_cache_gpu : public non_cache {
    //TODO: info was private, but I needed to access it in test_runner - reprotect?
    GPUNonCacheInfo info; //host struct of device pointers
    non_cache_gpu(szv_grid_cache& gcache, const grid_dims& gd_,
        const precalculate_gpu* p_, fl slope_);
    virtual void setSlope(fl sl);
    virtual ~non_cache_gpu();
    fl eval(const model& m, fl v) const;
    //evaluate the model on the gpu, v is the curl amount
    //sets m.minus_forces and returns total energy

    const GPUNonCacheInfo& get_info() const {
      return info;
    }
};

#endif

