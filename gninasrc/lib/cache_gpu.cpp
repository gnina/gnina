#include "cache_gpu.h"
#include "device_buffer.h"

void cache_gpu::populate(const model& m, const precalculate& p, const
        std::vector<smt>& atom_types_needed, grid& user_grid, bool
        display_progress) 
{
    const precalculate_splines* ps = dynamic_cast<const precalculate_splines*>(&p);
    cache::populate_deriv(m, *ps, atom_types_needed, user_grid, deriv_grids, display_progress);
    info.gridbegins = float3(gd[0].begin, gd[1].begin, gd[2].begin);
    info.gridends = float3(gd[0].end, gd[1].end, gd[2].end);
    info.slope = slope;
    info.num_movable_atoms = m.num_movable_atoms();
    thread_buffer.alloc(&info.types, sizeof(unsigned[info.num_movable_atoms]));
    std::vector<unsigned> movingtypes(info.num_movable_atoms);

    VINA_FOR(i, info.num_movable_atoms)
    {
        movingtypes[i] = m.atoms[i].get();
    }

    definitelyPinnedMemcpy(info.types, &movingtypes[0], sizeof(unsigned[info.num_movable_atoms]), cudaMemcpyHostToDevice);

    //set up grids
    info.ngrids = grids.size();
    std::vector<grid_gpu> gpu_grids;
    for (auto& g : grids) {
        gpu_grids.push_back(grid_gpu(g));
    }
    thread_buffer.alloc(&info.grids, sizeof(grid_gpu) * gpu_grids.size());
    definitelyPinnedMemcpy(info.grids, &gpu_grids[0], sizeof(grid_gpu) * gpu_grids.size(), cudaMemcpyHostToDevice);

    std::vector<grid_gpu> gpu_deriv_grids;
    for (auto& g : deriv_grids) {
        gpu_deriv_grids.push_back(grid_gpu(g));
    }
    thread_buffer.alloc(&info.deriv_grids, sizeof(grid_gpu) * gpu_deriv_grids.size());
    definitelyPinnedMemcpy(info.deriv_grids, &gpu_deriv_grids[0], sizeof(grid_gpu) * gpu_deriv_grids.size(), cudaMemcpyHostToDevice);
}

