#pragma once

#include "cache.h"
#include "gpu_math.h"
#include "grid_gpu.h"
#include "gpucode.h"
#include "precalculate_gpu.h"

struct cache_gpu : public cache {
    cache_gpu(const std::string& scoring_function_version_,
        const grid_dims& gd_, fl slope_, precalculate_gpu* prec)
        : cache(scoring_function_version_, gd_, slope_) {
      info.splineInfo = prec->getDeviceData();
      info.cutoff_sq = prec->cutoff_sqr();
    }

    virtual ~cache_gpu() {
    }
    virtual void populate(const model& m, const precalculate& p,
        const std::vector<smt>& atom_types_needed, grid& user_grid,
        bool display_progress = true);
    const GPUCacheInfo& get_info() const {
      return info;
    }

  private:
    GPUCacheInfo info; //host struct of device pointers
};
