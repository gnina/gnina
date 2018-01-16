#pragma once
#include "grid.h"
#include "gpu_math.h"

struct grid_gpu
{
	vec m_init;
	vec m_range;
	vec m_factor;
	vec m_dim_fl_minus_1;
	vec m_factor_inv;
    array3d_gpu<fl> data;
	array3d_gpu<fl> chargedata; //needs to be multiplied by atom charge

    grid_gpu(grid g) : m_init(g.m_init), m_range(g.m_range), 
                       m_factor(g.m_factor), m_dim_fl_minus_1(g.m_dim_fl_minus_1), 
                       m_factor_inv(g.m_factor_inv), data(g.data), 
                       chargedata(g.chargedata) {}
}
