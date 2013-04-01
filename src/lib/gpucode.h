/*
 * Header file for gpu code.
 */

#ifndef __GPUCODE_H
#define __GPUCODE_H
#ifdef SMINA_GPU
//everything is guarded by SMINA_GPU
// CUDA runtime
#include <cuda5/cuda_runtime.h>
#include <vector>

struct GPUSplineInfo
{
	unsigned n; //number of components
	float **splines; //pointer to spline data in device memory, size is number of components
	float fraction; //how spline is binned
	float cutoff; //where to stop

	GPUSplineInfo(): n(0), splines(NULL), fraction(0), cutoff(0) {}
};


struct GPUNonCacheInfo
{
	unsigned natoms, nrecatoms;
	float cutoff_sq;
	//device pointers for ligand data
	float *coords; //3*n - this is the only thing that actually changes
	float *charges; //n
	unsigned *types; //n
	float *minus_forces; //3*n
	float *energies; //per atom energy

	//device pointers for grid data
	float *gridends; //max range of grid
	float *gridbegins; //min range of grid

	//device pointers for receptor data
	float *recoords;
	float *reccharges;
	unsigned *rectypes;

	//triangular matrix of spline data, indexed by type, device pointer
	unsigned ntypes; //number of atom types; also, dimension of triangular splineInfo
	GPUSplineInfo *splineInfo;
};

void evaluate_splines_host(const GPUSplineInfo& spInfo, float r, std::vector<float>& vals, std::vector<float>& derivs);
float single_point_calc(GPUNonCacheInfo *dinfo, float *energies, float slope, unsigned natoms, unsigned nrecatoms, float v);


#endif //SMINA_GPU
#endif
