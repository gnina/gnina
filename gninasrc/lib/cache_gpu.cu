#include "cache_gpu.h"
#include "device_buffer.h"

cache_gpu::initialize(model& m)
{
    info.gridbegins = float3(gd[0].begin, gd[1].begin, gd[2].begin);
    info.gridends = float3(gd[0].end, gd[1].end, gd[2].end);
    info.slope = slope;
    info.num_movable_atoms = m.num_movable_atoms();

    //set up grids
    ngrids = grids.size();
    std::vector<grid_gpu> gpu_grids;
    for (auto& g : grids) {
        gpu_grids.push_back(gpu_grid(g));
    }
    thread_buffer.alloc(&info.grids, sizeof(grid_gpu) * gpu_grids.size());
}

__device__ void cache_gpu::eval(const gpu_data& gdata, fl v, 
        const GPUNonCacheInfo& dinfo) const {
	fl e = 0;
	sz nat = num_atom_types();

    unsigned idx = threadIdx.x;
    if (idx < dinfo.num_movable_atoms) {
        const atom_params& a = gdata.coords[idx];
        smt t = a.get();
        if (t < nat && !is_hydrogen(t)) {
            const grid& g = grids[t];
            assert(g.initialized());
            e += g.evaluate(a, gdata.coords[idx], slope, v);
        }
    }
    fl this_e = block_sum<fl>(e);
    if (threadIdx.x == 0)
        gdata.minus_forces[0].energy = this_e;
}

__device__ void eval_deriv(gpu_data& g, fl v, const grid& user_grid) const {
	fl e = 0;
	sz nat = num_atom_types();

    unsigned idx = threadIdx.x;
    if (idx < dinfo.num_movable_atoms) {
		const atom_params& a = gdata.coords[idx];
		smt t = a.get();
		if (t >= nat || is_hydrogen(t)) {
			gdata.minus_forces[i].energy = (fl)0;
            gdata.minus_forces[i].minus_force = zero_vec;
        }
        else {
		    const grid& g = grids[t];
		    assert(g.initialized());
		    vec deriv;
		    e += g.evaluate(a, gdata.coords[idx], slope, v, &deriv);
		    gdata.minus_forces[idx] = deriv;
        }
	}
    fl this_e = block_sum<fl>(e);
    if (threadIdx.x == 0)
        gdata.minus_forces[0].energy = this_e;
}
