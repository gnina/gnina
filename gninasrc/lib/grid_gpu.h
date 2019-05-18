#pragma once
#include "grid.h"
#include "gpu_math.h"

struct atom_params;
struct force_energy_tup;

struct grid_gpu {
    gfloat3 m_init;
    gfloat3 m_range;
    gfloat3 m_factor;
    gfloat3 m_dim_fl_minus_1;
    gfloat3 m_factor_inv;
    array3d_gpu<fl, fl> data;
    array3d_gpu<fl, fl> chargedata; //needs to be multiplied by atom charge

    grid_gpu(const grid& g)
        : m_init(g.m_init), m_range(g.m_range), m_factor(g.m_factor),
            m_dim_fl_minus_1(g.m_dim_fl_minus_1), m_factor_inv(g.m_factor_inv),
            data(g.data), chargedata(g.chargedata) {
    }

    __device__ void evaluate(const atom_params& a, float slope, float v,
        force_energy_tup& deriv) const;

    __device__ void evaluate_aux(const array3d_gpu<fl, fl>& data,
        const atom_params& a, float slope, float v,
        force_energy_tup& deriv) const;

    __device__ void evaluate_user(const atom_params& a, float slope,
        force_energy_tup& deriv) const;
};
