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

#include "parallel.h"
#include "parallel_mc.h"
#include "coords.h"
#include "parallel_progress.h"
#include "gpucode.h"
#include "device_buffer.h"
#include "non_cache_cnn.h"
#include "user_opts.h"

struct parallel_mc_task {
    model m;
    output_container out;
    rng generator;
    parallel_mc_task(const model& m_, int seed)
        : m(m_), generator(static_cast<rng::result_type>(seed)) {
      if (m_.gpu_initialized() && m.gdata.device_on) {
        //TODO: need to ensure that worker threads using these copies can't
        //deallocate GPU memory - race condition in
        //parallel_mc_aux::operator(), plus inefficiency of having to copy
        //buffers that won't change but could be freed
        m.gdata.coords = m_.gdata.coords;
        m.gdata.atom_coords = m_.gdata.atom_coords;
        m.gdata.minus_forces = m_.gdata.minus_forces;
        m.gdata.treegpu = m_.gdata.treegpu;
        m.gdata.interacting_pairs = m_.gdata.interacting_pairs;
        m.gdata.other_pairs = m_.gdata.other_pairs;
        m.gdata.dfs_order_bfs_indices = m_.gdata.dfs_order_bfs_indices;
        m.gdata.bfs_order_dfs_indices = m_.gdata.bfs_order_dfs_indices;
        m.gdata.coords_size = m_.gdata.coords_size;
        m.gdata.atom_coords_size = m_.gdata.atom_coords_size;
        m.gdata.forces_size = m_.gdata.forces_size;
        m.gdata.pairs_size = m_.gdata.pairs_size;
        m.gdata.other_pairs_size = m_.gdata.other_pairs_size;
        m.gdata.nlig_roots = m_.gdata.nlig_roots;
      }
    }
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

struct parallel_mc_aux {
    const monte_carlo* mc;
    const precalculate* p;
    igrid* ig;
    const vec* corner1;
    const vec* corner2;
    parallel_progress* pg;
    grid* user_grid;
    non_cache* nc;
    parallel_mc_aux(const monte_carlo* mc_, const precalculate* p_, igrid* ig_,
        const vec* corner1_, const vec* corner2_, parallel_progress* pg_,
        grid* user_grid_, non_cache* nc_)
        : mc(mc_), p(p_), ig(ig_), corner1(corner1_), corner2(corner2_),
            pg(pg_), user_grid(user_grid_), nc(nc_) {
    }

    void operator()(parallel_mc_task& t) const {
      //TODO: remove when the CNN is using the device buffer
      const non_cache_cnn* cnn = dynamic_cast<const non_cache_cnn*>(nc);
      if (t.m.gpu_initialized() && t.m.gdata.device_on && !cnn) { //only triggered with non-cnn gpu docking
        thread_buffer.reinitialize();
        //update our copy of gpu_data to have local buffers; N.B. this
        //theoretically could be restricted to fields we intend to modify,
        //but currently we also need copies of the things that are
        //deallocated in the model destructor
        atom_params *coords;
        thread_buffer.alloc(&coords,
            sizeof(atom_params[t.m.gdata.coords_size]));
        definitelyPinnedMemcpy(coords, t.m.gdata.coords,
            sizeof(atom_params[t.m.gdata.coords_size]),
            cudaMemcpyDeviceToDevice);
        t.m.gdata.coords = coords;

        vec *atom_coords;
        thread_buffer.alloc(&atom_coords,
            sizeof(vec[t.m.gdata.atom_coords_size]));
        definitelyPinnedMemcpy(atom_coords, t.m.gdata.atom_coords,
            sizeof(vec[t.m.gdata.atom_coords_size]), cudaMemcpyDeviceToDevice);
        t.m.gdata.atom_coords = atom_coords;

        force_energy_tup *minus_forces;
        thread_buffer.alloc(&minus_forces,
            sizeof(force_energy_tup[t.m.gdata.forces_size]));
        definitelyPinnedMemcpy(minus_forces, t.m.gdata.minus_forces,
            sizeof(force_energy_tup[t.m.gdata.forces_size]),
            cudaMemcpyDeviceToDevice);
        t.m.gdata.minus_forces = minus_forces;

        segment_node *device_nodes;
        gfloat4p *force_torques;

        //TODO: we really just want to copy data device-to-device
        tree_gpu old_tree;
        definitelyPinnedMemcpy(&old_tree, t.m.gdata.treegpu, sizeof(tree_gpu),
            cudaMemcpyDeviceToHost);
        unsigned num_nodes = old_tree.num_nodes;
        thread_buffer.alloc(&device_nodes, sizeof(segment_node[num_nodes]));
        thread_buffer.alloc(&force_torques, sizeof(gfloat4p[num_nodes]));
        definitelyPinnedMemcpy(device_nodes, old_tree.device_nodes,
            sizeof(segment_node[num_nodes]), cudaMemcpyDeviceToDevice);
        definitelyPinnedMemcpy(force_torques, old_tree.force_torques,
            sizeof(gfloat4p[num_nodes]), cudaMemcpyDeviceToDevice);
        tree_gpu* new_tree;
        thread_buffer.alloc(&new_tree, sizeof(tree_gpu));
        old_tree.device_nodes = device_nodes;
        old_tree.force_torques = force_torques;
        definitelyPinnedMemcpy(new_tree, &old_tree, sizeof(tree_gpu),
            cudaMemcpyHostToDevice);
        t.m.gdata.treegpu = new_tree;

        thread_buffer.alloc(&t.m.gdata.scratch, sizeof(float));

        size_t* dfs_order_bfs_indices = new size_t[num_nodes];
        memcpy(dfs_order_bfs_indices, t.m.gdata.dfs_order_bfs_indices,
            sizeof(size_t[num_nodes]));
        t.m.gdata.dfs_order_bfs_indices = dfs_order_bfs_indices;

        size_t* bfs_order_dfs_indices = new size_t[num_nodes];
        memcpy(bfs_order_dfs_indices, t.m.gdata.bfs_order_dfs_indices,
            sizeof(size_t[num_nodes]));
        t.m.gdata.bfs_order_dfs_indices = bfs_order_dfs_indices;
      }
      if (cnn && cnn->get_scorer().options().cnn_scoring>CNNrefinement) {
        CNNScorer cnn_scorer(cnn->get_scorer().options());
        const precalculate* p = cnn->get_precalculate();
        szv_grid_cache gridcache(t.m, p->cutoff_sqr());
        non_cache_cnn new_cnn(gridcache, cnn->get_grid_dims(), p,
            cnn->getSlope(), cnn_scorer);
        if (cnn_scorer.options().cnn_scoring==CNNmetropolisrefine||cnn_scorer.options().cnn_scoring==CNNmetropolisrescore)
        {
            (*mc)(t.m, t.out, *p, *ig, *corner1, *corner2, pg, t.generator,
                *user_grid, new_cnn);
        }
        else{//cnn_scoring=CNNall
            (*mc)(t.m, t.out, *p, new_cnn, *corner1, *corner2, pg, t.generator,
                *user_grid,new_cnn);
        }
      } else
        (*mc)(t.m, t.out, *p, *ig, *corner1, *corner2, pg, t.generator,
            *user_grid, *ig);
    }
};

//TODO: null model.gdata pointers at task exit

void merge_output_containers(const output_container& in, output_container& out,
    fl min_rmsd, sz max_size) {
  VINA_FOR_IN(i, in)
    add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_mc_task_container& many,
    output_container& out, fl min_rmsd, sz max_size) {
  min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
  VINA_FOR_IN(i, many) {
    merge_output_containers(many[i].out, out, min_rmsd, max_size);
  }
  out.sort();
}

void parallel_mc::operator()(const model& m, output_container& out,
    const precalculate& p, igrid& ig, const vec& corner1, const vec& corner2,
    rng& generator, grid& user_grid, non_cache& nc) const {
  parallel_progress pp;
  parallel_mc_aux parallel_mc_aux_instance(&mc, &p, &ig, &corner1, &corner2,
      (display_progress ? (&pp) : NULL), &user_grid, &nc);
  parallel_mc_task_container task_container;
  VINA_FOR(i, num_tasks)
    task_container.push_back(
        new parallel_mc_task(m, random_int(0, 1000000, generator)));
  if (display_progress) pp.init(num_tasks * mc.num_steps);

  auto thread_init = [&]()
  {
    initializeCUDA(m.gdata.device_id); //harmless to do if there isn't a gpu
    if (m.gdata.device_on)
    {
      const non_cache_cnn* cnn = dynamic_cast<const non_cache_cnn*>(&nc);
      if (!cnn)
        thread_buffer.init(available_mem(num_threads));
    }
  };

  parallel_iter<parallel_mc_aux,
  parallel_mc_task_container, parallel_mc_task,
      decltype(thread_init), true> parallel_iter_instance(
      &parallel_mc_aux_instance, num_threads, thread_init);
  parallel_iter_instance.run(task_container);

  merge_output_containers(task_container, out, mc.min_rmsd, mc.num_saved_mins);

}
