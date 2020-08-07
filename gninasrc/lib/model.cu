/*
 model.cu

 Contains the GPU-specific methods of model.
 Currently also have all eval methods (cpu or gpu) here for easy reference.

 */

#include "model.h"
#include "common.h"
#include "file.h"
#include "curl.h"
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include "non_cache_gpu.h"
#include "loop_timer.h"
#include "gpu_debug.h"
#include "device_buffer.h"

#define MAX_THREADS 1024

fl model::eval_interacting_pairs(const precalculate& p, fl v,
    const interacting_pairs& pairs, const vecv& coords) const { // clean up
  const fl cutoff_sqr = p.cutoff_sqr();
  fl e = 0;
  VINA_FOR_IN(i, pairs) {
    const interacting_pair& ip = pairs[i];
    fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);
    if (r2 < cutoff_sqr) {
      fl tmp = p.eval(atoms[ip.a], atoms[ip.b], r2);
      curl(tmp, v);
      e += tmp;
    }
  }
  return e;
}

fl model::eval_interacting_pairs_deriv(const precalculate& p, fl v,
    const interacting_pairs& pairs, const vecv& coords, vecv& forces) const { // adds to forces  // clean up
  const fl cutoff_sqr = p.cutoff_sqr();
  fl e = 0;
  VINA_FOR_IN(i, pairs) {
    const interacting_pair& ip = pairs[i];
    vec r;
    r = coords[ip.b] - coords[ip.a]; // a -> b
    fl r2 = sqr(r);
    if (r2 < cutoff_sqr) {
      pr tmp = p.eval_deriv(atoms[ip.a], atoms[ip.b], r2);
      vec force;
      force = tmp.second * r;
      curl(tmp.first, force, v);
      e += tmp.first;

      // FIXME inefficient, if using hard curl
      forces[ip.a] -= force; // we could omit forces on inflex here
      forces[ip.b] += force;
    }
  }
  return e;
}

//evaluates interacting pairs (which is all of them) on the gpu
template<typename infoT>
__host__  __device__ fl gpu_data::eval_interacting_pairs_deriv_gpu(
    const infoT& info, fl v, interacting_pair* pairs, unsigned pairs_sz) const { // adds to forces  // clean up

  float e = 0;
  // If there aren't any pairs, just return.
  if (pairs_sz == 0) {
    return e;
  }

  const fl cutoff_sqr = info.cutoff_sq;
#ifdef __CUDA_ARCH__
  scratch[0] = 0;
#else
  cudaMemsetAsync(scratch, 0, sizeof(float), cudaStreamPerThread);
#endif

  //TODO: this should not be a dynamic launch...
  if (pairs_sz < CUDA_THREADS_PER_BLOCK) {
    eval_intra_kernel<<<1, pairs_sz>>>(info.splineInfo, coords, pairs, pairs_sz,
        cutoff_sqr, v, minus_forces, scratch);

  } else {
    eval_intra_kernel<<<CUDA_GET_BLOCKS(pairs_sz, 1024),
    CUDA_THREADS_PER_BLOCK>>>(info.splineInfo, coords, pairs, pairs_sz,
        cutoff_sqr, v, minus_forces, scratch);
  }
#ifdef __CUDA_ARCH__ 
  cudaDeviceSynchronize();
  return scratch[0];
#else
  float out;
  definitelyPinnedMemcpy(&out, scratch, sizeof(float), cudaMemcpyDeviceToHost);
  return out;
#endif
}

fl model::evali(const precalculate& p, const vec& v) const { // clean up

  assert(0);
  /* TODO */
  /* fl e = 0; */
  /* VINA_FOR_IN(i, ligands) */
  /* 	e += eval_interacting_pairs(p, v[0], ligands[i].pairs, internal_coords); */
  //probably might as well use coords here
  /* return e; */
  return 0;
}

fl model::evale(const precalculate& p, const igrid& ig, const vec& v) const { // clean up
  fl e = ig.eval(*this, v[1]);
  e += eval_interacting_pairs(p, v[2], other_pairs, coords);
  return e;
}

fl model::eval(const precalculate& p, const igrid& ig, const vec& v,
    const conf& c, const grid& user_grid) { // clean up
  set(c);
  fl e = evale(p, ig, v);
  VINA_FOR_IN(i, ligands)
    e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords
  //std::cout << "smina_contribution: " << e << "\n";
  if (user_grid.initialized()) {
    fl uge = 0;
    vecv l_coords = this->get_ligand_coords();
    VINA_FOR_IN(i, l_coords) {
      fl tmp = user_grid.evaluate_user(l_coords[i], (fl) 1000);
      uge += tmp;
      e += tmp;
    }
    std::cout << "User_grid_contribution: " << uge << "\n";
  }

  return e;
}

static __global__
void derivatives_kernel(tree_gpu *t, const vec * coords, const vec* forces,
    change_gpu c) {
  t->derivative(coords, forces, &c);
}

static __global__
void set_conf_kernel(tree_gpu *t, const vec *atom_coords, vec *coords,
    const conf_gpu c) {
  t->set_conf(atom_coords, coords, &c);
}

template<typename infoT>
__device__ fl gpu_data::eval_deriv_gpu(const infoT& info, const vec& v,
    const conf_gpu& c, change_gpu& g) {
  // static loop_timer t;
  // t.resume();
  fl e, ie = 0;
  if (threadIdx.x == 0) {
    set_conf_kernel<<<1, treegpu->num_atoms>>>(treegpu, atom_coords,
        (vec*) coords, c);
    memset(minus_forces, 0, sizeof(force_energy_tup) * (forces_size));
    cudaDeviceSynchronize();
    e = single_point_calc(info, coords, minus_forces, v[1]);
    cudaDeviceSynchronize();
    if (other_pairs_size)
      e += eval_interacting_pairs_deriv_gpu(info, v[2], other_pairs,
          other_pairs_size);
    if (pairs_size)
      ie = eval_interacting_pairs_deriv_gpu(info, v[0], interacting_pairs,
          pairs_size); // adds to minus_forces

    e += ie;
    derivatives_kernel<<<1, treegpu->num_atoms>>>(treegpu, (vec*) coords,
        (vec*) minus_forces, g);

    cudaDeviceSynchronize();
  }

  // t.stop();

  /* flex.derivative(coords, minus_forces, g.flex); // inflex forces are ignored */
  return e;
}

//NB: the next two are only used for score_only, but they also compute the forces
//unnecessarily. this is fine as long as you don't care about using score_only
//for performance metrics, just correctness, which is what you should do. 
fl gpu_data::eval(const GPUNonCacheInfo& info, const float v) {
  fl e;
  cudaMemsetAsync(minus_forces, 0,
      sizeof(force_energy_tup) * info.num_movable_atoms, cudaStreamPerThread);
  e = single_point_calc(info, coords, minus_forces, v);
  return e;
}

fl gpu_data::eval_intramolecular(const GPUNonCacheInfo& info, const float v) {
  fl ie;
  ie = eval_interacting_pairs_deriv_gpu(info, v, interacting_pairs, pairs_size);
  return ie;
}

fl model::eval_deriv(const precalculate& p, const igrid& ig, const vec& v,
    const conf& c, change& g, const grid& user_grid) { // clean up
  static loop_timer t;
  t.resume();

  set(c);

  fl e = ig.eval_deriv(*this, v[1], user_grid); // sets minus_forces, except inflex
  fl ie = 0;

  if (!ig.skip_interacting_pairs()) {
    ie += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords,
        minus_forces); // adds to minus_forces

    VINA_FOR_IN(i, ligands)
      ie += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords,
          minus_forces); // adds to minus_forces
    e += ie;
  }

  // calculate derivatives
  ligands.derivative(coords, minus_forces, g.ligands);
  flex.derivative(coords, minus_forces, g.flex); // inflex forces are ignored
  g.receptor = rec_change; //for cnn
  t.stop();
  return e;
}

fl model::eval_intra(const precalculate& p, const vec& v) {
  fl ie = 0;
  VINA_FOR_IN(i, ligands)
    ie += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords,
        minus_forces); // adds to minus_forces
  return ie;
}

void model::clear_minus_forces() {
  minus_forces.clear();
  minus_forces.reserve(m_num_movable_atoms);
  VINA_FOR(i, m_num_movable_atoms) {
    vec force;
    force.data[0] = 0.0;
    force.data[1] = 0.0;
    force.data[2] = 0.0;
    minus_forces.push_back(force);
  }
}

void model::add_minus_forces(const std::vector<gfloat3>& forces) {
  assert(forces.size() <= m_num_movable_atoms);
  unsigned j = 0;
  VINA_FOR(i, m_num_movable_atoms) {
    if (!atoms[i].is_hydrogen()) // no hydrogen forces
    {
      minus_forces[i].data[0] += forces[j].x;
      minus_forces[i].data[1] += forces[j].y;
      minus_forces[i].data[2] += forces[j].z;
      j += 1;
    }
  }
}

void model::sub_minus_forces(const std::vector<gfloat3>& forces) {
  assert(forces.size() <= m_num_movable_atoms);
  unsigned j = 0;
  VINA_FOR(i, m_num_movable_atoms) {
    if (!atoms[i].is_hydrogen()) // no hydrogen forces
    {
      minus_forces[i].data[0] -= forces[j].x;
      minus_forces[i].data[1] -= forces[j].y;
      minus_forces[i].data[2] -= forces[j].z;
      j += 1;
    }
  }
}

void model::scale_minus_forces(fl scale) {
  unsigned j = 0;
  VINA_FOR(i, m_num_movable_atoms) {
    if (!atoms[i].is_hydrogen()) // no hydrogen forces
    {
      minus_forces[i].data[0] *= scale;
      minus_forces[i].data[1] *= scale;
      minus_forces[i].data[2] *= scale;
      j += 1;
    }
  }
}

/* This is my attempt to make gpus more deterministic.. */
void model::round_minus_forces() {
  VINA_FOR(i, m_num_movable_atoms) {
    if (!atoms[i].is_hydrogen()) // no hydrogen forces
    {
      for(unsigned j = 0; j < 3; j++) {
          double val = minus_forces[i].data[j];
          val *= 10000; 
          val = roundf(val);
          minus_forces[i].data[j] = val/10000;
      }
    }
  }
}

fl model::get_minus_forces_sum_magnitude() const {
  fl x = 0, y = 0, z = 0;
  VINA_FOR(i, m_num_movable_atoms) {
    if (!atoms[i].is_hydrogen()) // no hydrogen forces
    {
      x += minus_forces[i].data[0];
      y += minus_forces[i].data[1];
      z += minus_forces[i].data[2];
    }
  }
  return sqrt(x * x + y * y + z * z);
}

//evaluate interactiongs between all of flex (including rigid) and protein
//will ignore grid_atoms greater than max
fl model::eval_flex(const precalculate& p, const vec& v, const conf& c,
    unsigned maxGridAtom) {
  set(c);
  fl e = 0;
  sz nat = num_atom_types();
  const fl cutoff_sqr = p.cutoff_sqr();

  //ignore atoms after maxGridAtom (presumably part of "unfrag")
  sz gridstop = grid_atoms.size();
  if (maxGridAtom > 0 && maxGridAtom < gridstop) gridstop = maxGridAtom;

  // flex-rigid
  VINA_FOR(i, atoms.size()) {
    if (find_ligand(i) < ligands.size()) continue; // we only want flex-rigid interaction
    const atom& a = atoms[i];
    smt t1 = a.get();
    if (t1 >= nat) continue;
    VINA_FOR_IN(j, grid_atoms) {
      if (j >= gridstop) break;
      const atom& b = grid_atoms[j];
      smt t2 = b.get();
      if (t2 >= nat) continue;
      fl r2 = vec_distance_sqr(coords[i], b.coords);
      if (r2 < cutoff_sqr) {
        fl this_e = p.eval(a, b, r2);
        curl(this_e, v[1]);
        e += this_e;
      }
    }
  }

  return e;
}

fl model::eval_intramolecular(const precalculate& p, const vec& v,
    const conf& c) {
  set(c);
  fl e = 0;

  // internal for each ligand
  VINA_FOR_IN(i, ligands)
    e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

  sz nat = num_atom_types();
  const fl cutoff_sqr = p.cutoff_sqr();

  // flex-rigid
  VINA_FOR(i, num_movable_atoms()) {
    if (find_ligand(i) < ligands.size()) continue; // we only want flex-rigid interaction
    const atom& a = atoms[i];
    smt t1 = a.get();
    if (t1 >= nat || is_hydrogen(t1)) continue;
    VINA_FOR_IN(j, grid_atoms) {
      const atom& b = grid_atoms[j];
      smt t2 = b.get();
      if (t2 >= nat || is_hydrogen(t2)) continue;
      fl r2 = vec_distance_sqr(coords[i], b.coords);
      if (r2 < cutoff_sqr) {
        fl this_e = p.eval(a, b, r2);
        curl(this_e, v[1]);
        e += this_e;
      }
    }
  }

// flex-flex
  VINA_FOR_IN(i, other_pairs) {
    const interacting_pair& pair = other_pairs[i];
    if (find_ligand(pair.a) < ligands.size()
        || find_ligand(pair.b) < ligands.size()) continue; // we only need flex-flex
    fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
    if (r2 < cutoff_sqr) {
      fl this_e = p.eval(atoms[pair.a], atoms[pair.b], r2);
      curl(this_e, v[2]);
      e += this_e;
    }
  }
  return e;
}

fl model::eval_adjusted(const scoring_function& sf, const precalculate& p,
    const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy,
    const grid& user_grid) {
  fl e = eval(p, ig, v, c, user_grid); // sets c
  return sf.conf_independent(*this, e - intramolecular_energy);
}

void model::initialize_gpu() {
  //TODO: only re-malloc if need larger size
  deallocate_gpu();

  CUDA_CHECK_GNINA(device_malloc(&gdata.coords, sizeof(vec) * coords.size()));
  CUDA_CHECK_GNINA(
      device_malloc(&gdata.atom_coords, sizeof(vec) * atoms.size()));
  CUDA_CHECK_GNINA(
      device_malloc(&gdata.minus_forces, sizeof(vec) * minus_forces.size()));
  CUDA_CHECK_GNINA(device_malloc(&gdata.scratch, sizeof(float)));

  gdata.coords_size = coords.size();
  gdata.atom_coords_size = atoms.size();
  gdata.forces_size = minus_forces.size();

  //ligand internal pairs 
  if (ligands.size()) {
    std::vector<interacting_pair> ligand_pairs(ligands[0].pairs);
    for (int i = 1; i < ligands.size(); i++)
      ligand_pairs.insert(ligand_pairs.end(), ligands[i].pairs.begin(),
          ligands[i].pairs.end());
    gdata.pairs_size = ligand_pairs.size();

    CUDA_CHECK_GNINA(
        device_malloc(&gdata.interacting_pairs,
            sizeof(interacting_pair) * ligand_pairs.size()));
    CUDA_CHECK_GNINA(
        definitelyPinnedMemcpy(gdata.interacting_pairs, &ligand_pairs[0],
            sizeof(interacting_pair) * ligand_pairs.size(),
            cudaMemcpyHostToDevice));
  }

  //all flexible pairs, but not intra
  if (other_pairs.size()) {
    CUDA_CHECK_GNINA(
        device_malloc(&gdata.other_pairs,
            sizeof(interacting_pair) * other_pairs.size()));
    CUDA_CHECK_GNINA(
        definitelyPinnedMemcpy(gdata.other_pairs, &other_pairs[0],
            sizeof(interacting_pair) * other_pairs.size(),
            cudaMemcpyHostToDevice));
    gdata.other_pairs_size = other_pairs.size();
  }

  //input atom coords do not change
  std::vector<vec> acoords(atoms.size());
  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    acoords[i] = atoms[i].coords;
  }

  //set up tree. Writes padding to mark every atom in acoords with its owner.
  //TODO: quite intrusive
  tree_gpu tg(ligands, flex, &acoords[0], gdata.dfs_order_bfs_indices,
      gdata.bfs_order_dfs_indices);
  gdata.nlig_roots = tg.nlig_roots;
  CUDA_CHECK_GNINA(device_malloc(&gdata.treegpu, sizeof(tree_gpu)));
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(gdata.treegpu, &tg, sizeof(tree_gpu),
          cudaMemcpyHostToDevice));

  //this contains the marked_coords for all the atoms and therefore all the
  //nodes in both trees
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(gdata.atom_coords, &acoords[0],
          sizeof(vec) * atoms.size(), cudaMemcpyHostToDevice));

  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(gdata.coords, &coords[0],
          coords.size() * sizeof(atom_params), cudaMemcpyHostToDevice));
}

void model::deallocate_gpu() {
  if (gdata.coords) {
    CUDA_CHECK_GNINA(device_free(gdata.coords));
    gdata.coords = NULL;
  }
  if (gdata.atom_coords) {
    CUDA_CHECK_GNINA(device_free(gdata.atom_coords));
    gdata.atom_coords = NULL;
  }
  if (gdata.minus_forces) {
    CUDA_CHECK_GNINA(device_free(gdata.minus_forces));
    gdata.minus_forces = NULL;
  }
  if (gdata.treegpu) {
    tree_gpu::deallocate(gdata.treegpu);
    gdata.treegpu = NULL;
  }
  if (gdata.scratch) {
    CUDA_CHECK_GNINA(device_free(gdata.scratch));
    gdata.scratch = NULL;
  }
  if (gdata.dfs_order_bfs_indices) {
    delete[] gdata.dfs_order_bfs_indices;
    gdata.dfs_order_bfs_indices = NULL;
  }
  if (gdata.bfs_order_dfs_indices) {
    delete[] gdata.bfs_order_dfs_indices;
    gdata.bfs_order_dfs_indices = NULL;
  }

  if (gdata.interacting_pairs) {
    CUDA_CHECK_GNINA(device_free(gdata.interacting_pairs));
    gdata.interacting_pairs = NULL;
  }
  gdata.coords_size = gdata.atom_coords_size = gdata.forces_size =
      gdata.pairs_size = 0;
}

//copy relevant data to gpu buffers
void gpu_data::copy_to_gpu(model& m) {
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(coords, &m.coords[0], coords_size * sizeof(vec),
          cudaMemcpyHostToDevice));

  //minus_forces gets initialized in eval_deriv
  //interacting pairs, atom_coords and ligand tree do not change
}

//copy back relevant data from gpu buffers
void gpu_data::copy_from_gpu(model& m) {
  assert(coords);
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(&m.coords[0], coords, coords_size * sizeof(vec),
          cudaMemcpyDeviceToHost));
}

size_t gpu_data::node_idx_cpu2gpu(size_t cpu_idx) const {
  return dfs_order_bfs_indices[cpu_idx];
}

template
__host__  __device__ fl gpu_data::eval_interacting_pairs_deriv_gpu<
    GPUNonCacheInfo>(const GPUNonCacheInfo& info, fl v, interacting_pair* pairs,
    unsigned pairs_sz) const;

template
__host__  __device__ fl gpu_data::eval_interacting_pairs_deriv_gpu<GPUCacheInfo>(
    const GPUCacheInfo& info, fl v, interacting_pair* pairs,
    unsigned pairs_sz) const;

template __device__ fl gpu_data::eval_deriv_gpu<GPUNonCacheInfo>(
    const GPUNonCacheInfo& info, const vec& v, const conf_gpu& c,
    change_gpu& g);
template __device__ fl gpu_data::eval_deriv_gpu<GPUCacheInfo>(
    const GPUCacheInfo& info, const vec& v, const conf_gpu& c, change_gpu& g);

