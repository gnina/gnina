/*
 * Header file for gpu code.
 */

#ifndef __GPUCODE_H
#define __GPUCODE_H
// CUDA runtime (Metal build uses compatibility stubs instead)
#ifdef USE_METAL
  #include "cuda_metal_compat.h"
#else
  #include <cuda_runtime.h>
#endif
#include <vector>
#include "interacting_pairs.h"
#include "grid_gpu.h"

#ifdef USE_METAL
// Metal-friendly mirror of GPUSplineInfo — no pointer fields.
// Used by noncache_pair_energy / noncache_postprocess shaders.
// Must match SplineInfoMSL in gnina_kernels.metal (16 bytes).
struct GPUSplineInfo_MSL {
    unsigned n;        // number of active components (1..4)
    float    fraction; // knot bin width
    float    cutoff;   // cutoff distance
    unsigned _pad;
};
#endif // USE_METAL

struct GPUSplineInfo {
    unsigned n; //number of components
    float **splines; //pointer to spline data in device memory, size is number of components
    float fraction; //how spline is binned
    float cutoff; //where to stop

    GPUSplineInfo()
        : n(0), splines(NULL), fraction(0), cutoff(0) {
    }
};

/* float3 reads/writes can't be coalesced into a single load/store. But
 float4 ops can. Rather than padding 3-dimensional data to exploit this,
 pack it in with a relevant piece of 1-dimensional data. NB: without
 __align__, the compiler can't do this coalescing. */
struct __align__(sizeof(float4)) atom_params {
    gfloat3 coords;
    float charge;
};

struct __align__(sizeof(float4)) force_energy_tup {
    gfloat3 minus_force;
    float energy;

    __host__ __device__ force_energy_tup(void)
        : minus_force(0, 0, 0), energy(0) {
    }
    __host__ __device__ force_energy_tup(gfloat3 f, float e)
        : minus_force(f), energy(e) {
    }
    ;
    __host__ __device__ force_energy_tup(float f1, float f2, float f3, float f4)
        : minus_force(f1, f2, f3), energy(f4) {
    }
    ;

    __host__  __device__
   const fl& operator[](sz i) const {
      return i == 0 ? minus_force.x : i == 1 ? minus_force.y :
             i == 2 ? minus_force.z : energy;
    }
    __host__  __device__ fl& operator[](sz i) {
      return i == 0 ? minus_force.x : i == 1 ? minus_force.y :
             i == 2 ? minus_force.z : energy;
    }
};

inline __host__  __device__ force_energy_tup operator+(force_energy_tup& a,
    force_energy_tup& b) {
  return force_energy_tup(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]);
}
inline __host__  __device__ force_energy_tup& operator+=(force_energy_tup& a,
    force_energy_tup& b) {
  a = a + b;
  return a;
}

struct GPUNonCacheInfo {
    unsigned num_movable_atoms, nrec_atoms;
    float cutoff_sq;
    float slope;

    //device pointers for grid data
    gfloat3 gridends; //max range of grid
    gfloat3 gridbegins; //min range of grid

    //device pointers for ligand data
    force_energy_tup *lig_penalties;
    unsigned *types; //n

    //device pointers for receptor data
    atom_params *rec_atoms;
    unsigned *rectypes;

    //triangular matrix of spline data, indexed by type, device pointer
    unsigned ntypes; //number of atom types; also, dimension of triangular splineInfo
    GPUSplineInfo *splineInfo;

#ifdef USE_METAL
    // Standalone Metal buffers for spline data (owned by precalculate_gpu).
    // Stored as void* to avoid pulling in metal_context.h here.
    void*    metal_flat_buf;    // float[]          — flat knot data (5 floats/knot)
    void*    metal_info_buf;    // GPUSplineInfo_MSL[] — one per type-pair (tindex order)
    void*    metal_offset_buf;  // unsigned[]        — offset table [tindex*numc+c]
    unsigned metal_numc;        // spline components per type-pair

    // 1-float atomic energy accumulator in thread_buffer (unified memory).
    void*    metal_total_buf;   // MTLBufferHandle of thread_buffer slab
    size_t   metal_total_off;   // byte offset of the float within that buffer
    float*   metal_total_scratch; // CPU pointer (same location, for zeroing/reading)
#endif
};

struct GPUCacheInfo {
    gfloat3 gridends;
    gfloat3 gridbegins;
    fl slope;
    float cutoff_sq;
    unsigned num_movable_atoms;

    //lig atom types
    unsigned *types;
    //grids used to interpolate atom energies
    grid_gpu* grids;
    unsigned ngrids;
    GPUSplineInfo *splineInfo;
};

void evaluate_splines_host(const GPUSplineInfo& spInfo, float r,
    float *device_vals, float *device_derivs);

__host__ __device__
float single_point_calc(const GPUNonCacheInfo &dinfo, atom_params *lig,
    force_energy_tup *out, float v);

__host__ __device__
float single_point_calc(const GPUCacheInfo &dinfo, atom_params *lig,
    force_energy_tup *out, float v);

#ifndef USE_METAL
__global__
void eval_intra_kernel(const GPUSplineInfo * spinfo, const atom_params * atoms,
    const interacting_pair* pairs, unsigned npairs, float cutoff_sqr, float v,
    force_energy_tup *out, float *e);
#endif // USE_METAL

#ifdef USE_METAL
// Under Metal, eval_deriv_gpu is a plain C++ function (defined in gpucode.cu).
// Declaration here allows model.cu to call it directly for CPU pair evaluation.
float eval_deriv_gpu(const GPUSplineInfo* splineInfo, unsigned t, float charge,
    unsigned rt, float rcharge, float r2, float& dor);
#endif // USE_METAL

#endif
