#pragma once

#include "cache.h"
#include "gpu_math.h"
#include "grid_gpu.h"
#include "gpucode.h"

struct cache_gpu : public cache{
	cache_gpu(const std::string& scoring_function_version_, const grid_dims& gd_, fl slope_) : 
        cache(scoring_function_version_, gd_, slope);
    virtual ~cache_gpu() {}
    const GPUCacheInfo& get_info() const { return info; }

    private:
    GPUCacheInfo info; //host struct of device pointers
};
