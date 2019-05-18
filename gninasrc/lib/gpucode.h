/*
 * Header file for gpu code.
 */

#ifndef __GPUCODE_H
#define __GPUCODE_H
// CUDA runtime
#include <cuda_runtime.h>
#include <vector>
#include "interacting_pairs.h"
#include "grid_gpu.h"

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

__global__
void eval_intra_kernel(const GPUSplineInfo * spinfo, const atom_params * atoms,
    const interacting_pair* pairs, unsigned npairs, float cutoff_sqr, float v,
    force_energy_tup *out, float *e);

#endif
