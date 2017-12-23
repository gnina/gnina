#pragma once

#include "cache.h"
#include "gpu_math.h"
#include "grid.h"

struct cache_gpu {
    cache_gpu(const std::string& scoring_function_version_, 
            const grid_dims& gd_, fl slope_);
    __device__ fl eval (const gpu_data& gdata, fl v) const;
    __device__ fl eval_deriv(gpu_data& g, fl v, const grid& user_grid) const;
    void copy_to_gpu(cache c) const;
    ~cache_gpu() {}
private:
    char* scoring_function_version;
    atom_params* atoms;
    sz natoms;
    grid_dims gd;
    fl slope;
    grid* grids;
    sz ngrids;
};
