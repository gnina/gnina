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

struct parallel_mc_task
{
	model m;
	output_container out;
	rng generator;
	parallel_mc_task(const model& m_, int seed) :
			m(m_), generator(static_cast<rng::result_type>(seed))
	{
        if (run_on_gpu) {
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
        }
	}
};

typedef boost::ptr_vector<parallel_mc_task> parallel_mc_task_container;

struct parallel_mc_aux
{
	const monte_carlo* mc;
	const precalculate* p;
	const igrid* ig;
	const vec* corner1;
	const vec* corner2;
	parallel_progress* pg;
	grid* user_grid;
	parallel_mc_aux(const monte_carlo* mc_, const precalculate* p_,
			const igrid* ig_, const vec* corner1_, const vec* corner2_,
			parallel_progress* pg_, grid* user_grid_)
	:
			mc(mc_), p(p_), ig(ig_), corner1(corner1_), corner2(corner2_), pg(
					pg_), user_grid(user_grid_)
	{
	}
	void operator()(parallel_mc_task& t) const
	{
        if (run_on_gpu) {
            thread_buffer.reinitialize();
            //update our copy of gpu_data to have local buffers for the fields
            //we intend to modify
            atom_params *coords;
            thread_buffer.alloc(&coords, sizeof(t.m.gdata.coords[t.m.gdata.coords_size]));
            definitelyPinnedMemcpy(coords, t.m.gdata.coords, sizeof(t.m.gdata.coords[t.m.gdata.coords_size]), cudaMemcpyDeviceToDevice);
            t.m.gdata.coords = coords;

            force_energy_tup *minus_forces;
            thread_buffer.alloc(&minus_forces, sizeof(t.m.gdata.minus_forces[t.m.gdata.forces_size]));
            definitelyPinnedMemcpy(minus_forces, t.m.gdata.minus_forces, sizeof(t.m.gdata.minus_forces[t.m.gdata.forces_size]), cudaMemcpyDeviceToDevice);
            t.m.gdata.minus_forces = minus_forces;

            segment_node *device_nodes;
            gfloat4p *force_torques;

            //TODO: we really just want to copy data device-to-device
            tree_gpu old_tree;
            definitelyPinnedMemcpy(&old_tree, t.m.gdata.treegpu, sizeof(tree_gpu), cudaMemcpyDeviceToHost);
            unsigned num_nodes = old_tree.num_nodes;
            thread_buffer.alloc(&device_nodes, sizeof(old_tree.device_nodes[num_nodes]));
            thread_buffer.alloc(&force_torques, sizeof(old_tree.force_torques[num_nodes]));
            definitelyPinnedMemcpy(device_nodes, old_tree.device_nodes, sizeof(old_tree.device_nodes[num_nodes]), cudaMemcpyDeviceToDevice);
            definitelyPinnedMemcpy(force_torques, old_tree.force_torques, sizeof(old_tree.force_torques[num_nodes]), cudaMemcpyDeviceToDevice);
            tree_gpu* new_tree;
            thread_buffer.alloc(&new_tree, sizeof(tree_gpu));
            old_tree.device_nodes = device_nodes;
            old_tree.force_torques = force_torques;
            definitelyPinnedMemcpy(new_tree, &old_tree, sizeof(tree_gpu), cudaMemcpyHostToDevice);
            t.m.gdata.treegpu = new_tree;

            thread_buffer.alloc(&t.m.gdata.scratch, sizeof(float));
        }

		(*mc)(t.m, t.out, *p, *ig, *corner1, *corner2, pg, t.generator, *user_grid);
	}
};

void merge_output_containers(const output_container& in, output_container& out,
		fl min_rmsd, sz max_size)
{
	VINA_FOR_IN(i, in)
	add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_mc_task_container& many,
		output_container& out, fl min_rmsd, sz max_size)
{
	min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
	VINA_FOR_IN(i, many)
	merge_output_containers(many[i].out, out, min_rmsd, max_size);
	out.sort();
}

void parallel_mc::operator()(const model& m, output_container& out,
		const precalculate& p, const igrid& ig, const vec& corner1,
		const vec& corner2, rng& generator, grid& user_grid) const
{
	parallel_progress pp;
	parallel_mc_aux parallel_mc_aux_instance(&mc, &p, &ig, &corner1, &corner2,
			(display_progress ? (&pp) : NULL), &user_grid);
	parallel_mc_task_container task_container;
	VINA_FOR(i, num_tasks)
		task_container.push_back(new parallel_mc_task(m, random_int(0, 1000000, generator)));
	if (display_progress)
		pp.init(num_tasks * mc.num_steps);
	parallel_iter<parallel_mc_aux, parallel_mc_task_container, parallel_mc_task,
			true> parallel_iter_instance(&parallel_mc_aux_instance,
			num_threads);
	parallel_iter_instance.run(task_container);
	merge_output_containers(task_container, out, mc.min_rmsd,
			mc.num_saved_mins);
}
