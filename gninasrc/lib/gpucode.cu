#include <stdio.h>
#include "gpu_util.h"
#include "gpucode.h"
#include "curl.h"

#define THREADS_PER_BLOCK 1024
#define warpSize 32

__global__
void evaluate_splines(float **splines, float r, float fraction, float cutoff,
    float *vals, float *derivs) {
  unsigned i = blockIdx.x;
  float *spline = splines[i];
  vals[i] = 0;
  derivs[i] = 0;

  if (r >= cutoff || r < 0) {
    return;
  }

  unsigned index = r / fraction; //xval*numpoints/cutoff
  unsigned base = 5 * index;
  float x = spline[base];
  float a = spline[base + 1];
  float b = spline[base + 2];
  float c = spline[base + 3];
  float d = spline[base + 4];

  const float lx = r - x;
  vals[i] = ((a * lx + b) * lx + c) * lx + d;
  derivs[i] = (3 * a * lx + 2 * b) * lx + c;
}

//evaluate a single spline
__device__
float evaluate_spline(float *spline, float r, float fraction, float cutoff,
    float& deriv) {
  float val = 0;
  deriv = 0;
  if (r >= cutoff || r < 0) {
    return 0;
  }

  unsigned index = r / fraction; //xval*numpoints/cutoff
  unsigned base = 5 * index;
  float x = spline[base];
  float a = spline[base + 1];
  float b = spline[base + 2];
  float c = spline[base + 3];
  float d = spline[base + 4];

  const float lx = r - x;
  val = ((a * lx + b) * lx + c) * lx + d;
  deriv = (3 * a * lx + 2 * b) * lx + c;
  return val;
}

void evaluate_splines_host(const GPUSplineInfo& spInfo, float r,
    float *device_vals, float *device_derivs) {
  unsigned n = spInfo.n;
  evaluate_splines<<<n, 1>>>((float**) spInfo.splines, r, spInfo.fraction,
      spInfo.cutoff, device_vals, device_derivs);
}

__device__
float eval_deriv_gpu(const GPUSplineInfo* splineInfo, unsigned t, float charge,
    unsigned rt, float rcharge, float r2, float& dor) {
  float r = sqrt(r2);
  unsigned t1, t2;
  float charge1, charge2;
  if (t < rt) {
    t1 = t;
    t2 = rt;
    charge1 = fabs(charge);
    charge2 = fabs(rcharge);
  } else {
    t1 = rt;
    t2 = t;
    charge1 = fabs(rcharge);
    charge2 = fabs(charge);
  }

  unsigned tindex = t1 + t2 * (t2 + 1) / 2;
  GPUSplineInfo spInfo = splineInfo[tindex];
  unsigned n = spInfo.n; //number of components

  float ret = 0, d = 0;

  //ick, hard code knowledge of components here; need to come up with
  //something mroe elegant
  //TypeDependentOnly,//no need to adjust by charge
  if (n > 0) {
    float fraction = spInfo.fraction;
    float cutoff = spInfo.cutoff;
    float val, deriv;
    val = evaluate_spline(spInfo.splines[0], r, fraction, cutoff, deriv);
    ret += val;
    d += deriv;
    //AbsAChargeDependent,//multiply by the absolute value of a's charge
    if (n > 1) {
      val = evaluate_spline(spInfo.splines[1], r, fraction, cutoff, deriv);
      ret += val * charge1;
      d += deriv * charge1;
      //AbsBChargeDependent,//multiply by abs(b)'s charge
      if (n > 2) {
        val = evaluate_spline(spInfo.splines[2], r, fraction, cutoff, deriv);
        ret += val * charge2;
        d += deriv * charge2;
        //ABChargeDependent,//multiply by a*b
        if (n > 3) {
          val = evaluate_spline(spInfo.splines[3], r, fraction, cutoff, deriv);
          ret += val * charge2 * charge1;
          d += deriv * charge2 * charge1;
        }
      }
    }
  }

  dor = d / r; //divide by distance to normalize vector later
  return ret;
}

template<typename T> T __device__ __host__ zero(void);
template<> gfloat3 zero(void) {
  return gfloat3(0, 0, 0);
}

template<> float zero(void) {
  return 0;
}

/* Can't really return an accurate value. Don't need it anyway. */
__device__
void xadd(force_energy_tup *a, const force_energy_tup &b) {
  atomicAdd(&a->minus_force.x, b.minus_force.x);
  atomicAdd(&a->minus_force.y, b.minus_force.y);
  atomicAdd(&a->minus_force.z, b.minus_force.z);
  atomicAdd(&a->energy, b.energy);
}

//device functions for warp-based reduction using shufl operations
template<class T>
__device__   __forceinline__ T warp_sum(T mySum) {
  for (int offset = warpSize >> 1; offset > 0; offset >>= 1)
    mySum += shuffle_down(mySum, offset);
  return mySum;
}

__device__ __forceinline__
bool isNotDiv32(unsigned int val) {
  return val & 31;
}

__device__
void operator+=(force_energy_tup &a, const force_energy_tup &b) {
  a.minus_force += b.minus_force;
  a.energy += b.energy;
}

/* requires blockDim.x <= 1024, blockDim.y == 1 */
template<class T>
__device__   __forceinline__ T block_sum(T mySum) {
  const unsigned int lane = threadIdx.x & 31;
  const unsigned int wid = threadIdx.x >> 5;

  __shared__ T scratch[32];

  mySum = warp_sum(mySum);
  if (lane == 0) scratch[wid] = mySum;
  __syncthreads();

  if (wid == 0) {
    mySum = (threadIdx.x < blockDim.x >> 5) ? scratch[lane] : zero<T>();
    mySum = warp_sum(mySum);
    if (threadIdx.x == 0 && isNotDiv32(blockDim.x))
      mySum += scratch[blockDim.x >> 5];
  }
  return mySum;
}

__device__
void eval_intra_st(const GPUSplineInfo * spinfo, const atom_params * atoms,
    const interacting_pair* pairs, unsigned npairs, float cutoff_sqr, float v,
    float *st_e) {
  float total = 0.0;
  for (unsigned i = 0; i < npairs; i += 1) {

    const interacting_pair& ip = pairs[i];
    gfloat3 r = atoms[ip.b].coords - atoms[ip.a].coords;
    float r2 = dot(r, r);
    if (r2 < cutoff_sqr) {
      gfloat3 deriv(0, 0, 0);
      float dor;
      unsigned t1 = ip.t1;
      unsigned t2 = ip.t2;
      float energy = eval_deriv_gpu(spinfo, t1, atoms[ip.a].charge, t2,
          atoms[ip.b].charge, r2, dor);
      deriv = r * dor;
      curl(energy, deriv, v);

      // st_out[ip.b].minus_force.x += deriv.x;
      // st_out[ip.b].minus_force.y += deriv.y;
      // st_out[ip.b].minus_force.z += deriv.z;
// 
      // st_out[ip.a].minus_force.x -= deriv.x;
      // st_out[ip.a].minus_force.y -= deriv.y;
      // st_out[ip.a].minus_force.z -= deriv.z;

      total += energy;
    }
  }
  *st_e = total;
}

//calculates the energies of all ligand-prot interactions and combines the results
//into energies and minus forces
//needs enough shared memory for derivatives and energies of single ligand atom
//roffset specifies how far into the receptor atoms we are
template<bool remainder> __global__
void interaction_energy(const GPUNonCacheInfo dinfo, unsigned remainder_offset,
    const atom_params *ligs, force_energy_tup *out) {
  unsigned l = blockIdx.x;
  unsigned r = blockDim.x - threadIdx.x - 1;
  unsigned roffset =
      remainder ? remainder_offset : blockIdx.y * THREADS_PER_BLOCK;
  unsigned ridx = roffset + r;
  //get ligand atom info
  unsigned t = dinfo.types[l];
  float rec_energy = 0;
  gfloat3 rec_deriv(0, 0, 0);

  atom_params lin = ligs[l];
  gfloat3 xyz = lin.coords;

  //TODO: remove hydrogen atoms completely
  if (t > 1) { //ignore hydrogens
    //now consider interaction with every possible receptor atom
    if (ridx < dinfo.nrec_atoms) {
      //compute squared difference
      atom_params rin = dinfo.rec_atoms[ridx];
      gfloat3 diff = xyz - rin.coords;

      float rSq = 0;
      for (unsigned j = 0; j < 3; j++) {
        float d = diff[j];
        rSq += d * d;
      }

      if (rSq < dinfo.cutoff_sq) {
        //dkoes - the "derivative" value returned by eval_deriv
        //is normalized by r (dor = derivative over r?)
        float dor;
        rec_energy = eval_deriv_gpu(dinfo.splineInfo, t, lin.charge,
            dinfo.rectypes[ridx], rin.charge, rSq, dor);
        rec_deriv = diff * dor;
      }
    }
  }
  float this_e = block_sum<float>(rec_energy);
  gfloat3 deriv = block_sum<gfloat3>(rec_deriv);
  if (threadIdx.x == 0) {
    if (dinfo.nrec_atoms > 1024)
      pseudoAtomicAdd(&out[l], force_energy_tup(deriv, this_e));
    else
      out[l] += force_energy_tup(deriv, this_e);
  }
}

__device__ void reduce_energy(force_energy_tup *result, float energy) {
  unsigned idx = threadIdx.x;
  float e = block_sum<float>(energy);
  if (idx == 0) result[0].energy = e;
}

/*Cached version of GPU energy/deriv eval*/
__device__
void interaction_energy(const GPUCacheInfo& dinfo, const atom_params *ligs,
    force_energy_tup *out, float v) {
  force_energy_tup val = force_energy_tup(0, 0, 0, 0);

  unsigned idx = threadIdx.x;
  if (idx < dinfo.num_movable_atoms) {
    const atom_params& a = ligs[idx];
    unsigned t = dinfo.types[idx];
    if (t > 1) {
      const grid_gpu& g = dinfo.grids[t];
      g.evaluate(a, dinfo.slope, v, val);
      out[idx] = val;
    }
  }
  reduce_energy(out, val.energy);
}

__global__ void reduce_energy(force_energy_tup *result, int n, float v,
    gfloat3 gridbegins, gfloat3 gridends, float slope, const atom_params* ligs) {
  int l = threadIdx.x;
  force_energy_tup val = force_energy_tup(0, 0, 0, 0);
  if (l < n) {
    val = result[l];
    curl(val.energy, val.minus_force, v);
    gfloat3 xyz = ligs[l].coords;
    gfloat3 out_of_bounds_deriv(0, 0, 0);
    float out_of_bounds_penalty = 0;
    //evaluate for out of boundsness
    for (unsigned i = 0; i < 3; i++) {
      float min = gridbegins[i];
      float max = gridends[i];
      if (xyz[i] < min) {
        out_of_bounds_deriv[i] = -1;
        out_of_bounds_penalty += fabs(min - xyz[i]);
        xyz[i] = min;
      } else
        if (xyz[i] > max) {
          out_of_bounds_deriv[i] = 1;
          out_of_bounds_penalty += fabs(max - xyz[i]);
          xyz[i] = max;
        }
      out_of_bounds_deriv[i] *= slope;
    }

    out_of_bounds_penalty *= slope;
    val.energy = val.energy + out_of_bounds_penalty;
    val.minus_force = val.minus_force + out_of_bounds_deriv;
    result[l] = val;
  }
  float e = block_sum<float>(val.energy);
  if (l == 0) {
    result[0].energy = e;
  }
}

/* global */
/* void apply_penalties(float3 *coords, float3 gridBegins, float3 gridEnds, */
/*                      float3 *e_penalties, float3 *deriv_penalties){ */

/* } */

__host__ __device__
float single_point_calc(const GPUNonCacheInfo &info, atom_params *ligs,
    force_energy_tup *out, float v) {
  /* Assumed by warp_sum */
  assert(THREADS_PER_BLOCK <= 1024);
  unsigned num_movable_atoms = info.num_movable_atoms;
  unsigned nrec_atoms = info.nrec_atoms;

  //this will calculate the per-atom energies and forces.
  //TODO: there could be one execution stream for the blocks with
  //a full complement of threads and a separate stream
  //for the blocks that have the remaining threads
  unsigned nfull_blocks = nrec_atoms / THREADS_PER_BLOCK;
  unsigned nthreads_remain = nrec_atoms % THREADS_PER_BLOCK;

  if (nfull_blocks)
    interaction_energy<0> <<<dim3(num_movable_atoms, nfull_blocks),
    THREADS_PER_BLOCK>>>(info, 0, ligs, out);
  if (nthreads_remain)
    interaction_energy<1> <<<num_movable_atoms, ROUND_TO_WARP(nthreads_remain)>>>(
        info, nrec_atoms - nthreads_remain, ligs, out);

  //TODO: reduce energy only launches one block, thus enforcing the
  //hardware restriction on the number of threads per block for the number of
  //movable atoms. generalize this to remove this unnecessary constraint
  assert(num_movable_atoms <= 1024);
  reduce_energy<<<1, ROUND_TO_WARP(num_movable_atoms)>>>(out, num_movable_atoms,
      v, info.gridbegins, info.gridends, info.slope, ligs);
  abort_on_gpu_err();
  cudaDeviceSynchronize();
#ifdef __CUDA_ARCH__
  return out->energy;
#else
  float cpu_out;
  definitelyPinnedMemcpy(&cpu_out, &out->energy, sizeof(float), cudaMemcpyDeviceToHost);
  return cpu_out;
#endif
}

__global__
void cache_gpu_kernel(const GPUCacheInfo info, atom_params *ligs,
    force_energy_tup *out, float v) {
  interaction_energy(info, ligs, out, v);
}

__host__ __device__
float single_point_calc(const GPUCacheInfo &info, atom_params *ligs,
    force_energy_tup *out, float v) {
#ifdef __CUDA_ARCH__
  /* Assumed by warp_sum */
  assert(THREADS_PER_BLOCK <= 1024);
  interaction_energy(info, ligs, out, v);
  return out->energy;
#else
  /*If we're on the CPU we need to launch a kernel to do this eval*/
  cache_gpu_kernel<<<1, ROUND_TO_WARP(info.num_movable_atoms)>>>(info, ligs, out, v);
  abort_on_gpu_err();
  float cpu_out;
  definitelyPinnedMemcpy(&cpu_out, &out->energy, sizeof(float), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  return cpu_out;
#endif
}

/* evaluate contribution of interacting pairs, add to forces and place total */
/* energy in e (which must be zero initialized). */

__global__
void eval_intra_kernel(const GPUSplineInfo * spinfo, const atom_params * atoms,
    const interacting_pair* pairs, unsigned npairs, float cutoff_sqr, float v,
    force_energy_tup *out, float *e) {

  for (unsigned i = blockDim.x * blockIdx.x + threadIdx.x; i < npairs; i +=
  CUDA_THREADS_PER_BLOCK * gridDim.x) {

    const interacting_pair& ip = pairs[i];
    gfloat3 r = atoms[ip.b].coords - atoms[ip.a].coords;
    float r2 = dot(r, r);
    if (r2 < cutoff_sqr) {
      gfloat3 deriv(0, 0, 0);
      float dor;
      unsigned t1 = ip.t1;
      unsigned t2 = ip.t2;
      float energy = eval_deriv_gpu(spinfo, t1, atoms[ip.a].charge, t2,
          atoms[ip.b].charge, r2, dor);
      deriv = r * dor;
      curl(energy, deriv, v);

      atomicAdd(&out[ip.b].minus_force.x, deriv.x);
      atomicAdd(&out[ip.b].minus_force.y, deriv.y);
      atomicAdd(&out[ip.b].minus_force.z, deriv.z);

      atomicAdd(&out[ip.a].minus_force.x, -deriv.x);
      atomicAdd(&out[ip.a].minus_force.y, -deriv.y);
      atomicAdd(&out[ip.a].minus_force.z, -deriv.z);

      atomicAdd(e, energy); //can't use blocksum unless we have multiple of 32 threads
      /*float this_e = block_sum<float>(energy);
       if (threadIdx.x == 0)
       {
       printf("Intra %f\n",this_e);
       *e += this_e;
       }*/
    }
  }
}

void *getHostMem() {
  void *r = nullptr;
  CUDA_CHECK_GNINA(cudaHostAlloc(&r, 40960 * 1024, cudaHostAllocDefault));
  return r;
}

cudaError definitelyPinnedMemcpy(void* dst, const void *src, size_t n,
    cudaMemcpyKind k) {
  assert(k != cudaMemcpyDefault);
  if (k == cudaMemcpyDeviceToDevice || k == cudaMemcpyHostToHost)
    return cudaMemcpy(dst, src, n, k);
  thread_local void *buf = getHostMem();
  assert(n < 40960 * 1024);
  if (k == cudaMemcpyHostToDevice) {
    memcpy(buf, src, n);
    CUDA_CHECK_GNINA(cudaMemcpyAsync(dst, buf, n, k, cudaStreamPerThread));
    CUDA_CHECK_GNINA(cudaStreamSynchronize(cudaStreamPerThread));
  } else {
    assert(k == cudaMemcpyDeviceToHost);

    CUDA_CHECK_GNINA(cudaMemcpyAsync(buf, src, n, k, cudaStreamPerThread));
    CUDA_CHECK_GNINA(cudaStreamSynchronize(cudaStreamPerThread));
    memcpy(dst, buf, n);
  }
  return cudaSuccess;
}

