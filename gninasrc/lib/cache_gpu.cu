#include "cache_gpu.h"
#include "device_buffer.h"

cache_gpu::initialize(model& m)
{
    info.gridbegins = float3(gd[0].begin, gd[1].begin, gd[2].begin);
    info.gridends = float3(gd[0].end, gd[1].end, gd[2].end);
    info.slope = slope;
    info.num_movable_atoms = m.num_movable_atoms();
    thread_buffer.alloc(&info.types, sizeof(unsigned[num_movable_atoms]));

    //set up grids
    ngrids = grids.size();
    std::vector<grid_gpu> gpu_grids;
    for (auto& g : grids) {
        gpu_grids.push_back(gpu_grid(g));
    }
    thread_buffer.alloc(&info.grids, sizeof(grid_gpu) * gpu_grids.size());
}

