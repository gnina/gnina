/* dkoes
 * This file contains all the standalone gpu kernels.  There is (hopefully)
 * a nicer way to organize this, but I'm currently slightly flummoxed as to
 * how to cleaning mix object-oriented cpu and gpu code.
 */
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <stdio.h>
#include "gpu_util.h"
#include "gpucode.h"

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

//TODO: buy compute 3.0 or greater card and implement dynamic paralellism
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
			val = evaluate_spline(spInfo.splines[1], r, fraction, cutoff,
					deriv);
			ret += val * charge1;
			d += deriv * charge1;
			//AbsBChargeDependent,//multiply by abs(b)'s charge
			if (n > 2) {
				val = evaluate_spline(spInfo.splines[2], r, fraction, cutoff,
						deriv);
				ret += val * charge2;
				d += deriv * charge2;
				//ABChargeDependent,//multiply by a*b
				if (n > 3) {
					val = evaluate_spline(spInfo.splines[3], r, fraction,
							cutoff, deriv);
					ret += val * charge2 * charge1;
					d += deriv * charge2 * charge1;
				}
			}
		}
	}

	dor = d / r; //divide by distance to normalize vector later
	return ret;
}

//curl function to scale back positive energies and match vina calculations
//assume v is reasonable
__device__
void curl(float& e, float *deriv, float v) {
	if (e > 0) {
		float tmp = (v / (v + e));
		e *= tmp;
		tmp *= tmp;
		for (unsigned i = 0; i < 3; i++)
			deriv[i] *= tmp;
	}
}

template<typename T> T __device__ __host__ zero(void);
template<> float3 zero(void) {
	return float3(0, 0, 0);
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
__device__  __forceinline__ T warp_sum(T mySum) {
	for (int offset = warpSize >> 1; offset > 0; offset >>= 1)
		mySum += __shfl_down(mySum, offset);
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
__device__  __forceinline__ T block_sum(T mySum) {
	const unsigned int lane = threadIdx.x & 31;
	const unsigned int wid = threadIdx.x >> 5;

	__shared__ T scratch[32];

	mySum = warp_sum(mySum);
	if (lane == 0)
		scratch[wid] = mySum;
	__syncthreads();

	if (wid == 0) {
		mySum = (threadIdx.x < blockDim.x >> 5) ? scratch[lane] : zero<T>();
		mySum = warp_sum(mySum);
		if (threadIdx.x == 0 && isNotDiv32(blockDim.x))
			mySum += scratch[blockDim.x >> 5];
	}
	return mySum;
}

//calculates the energies of all ligand-prot interactions and combines the results
//into energies and minus forces
//needs enough shared memory for derivatives and energies of single ligand atom
//roffset specifies how far into the receptor atoms we are
template<bool remainder> __global__
void interaction_energy(
		const GPUNonCacheInfo dinfo, //intentionally copying from host to device
		unsigned remainder_offset, float slope, float v,
		const atom_params *ligs, force_energy_tup *out) {
	unsigned l = blockIdx.x;
	unsigned r = blockDim.x - threadIdx.x - 1;
	unsigned roffset =
			remainder ? remainder_offset : blockIdx.y * THREADS_PER_BLOCK;
	unsigned ridx = roffset + r;
	//get ligand atom info
	unsigned t = dinfo.types[l];
	float rec_energy = 0;
	float3 rec_deriv = float3(0, 0, 0);
	float3 out_of_bounds_deriv = float3(0, 0, 0);
	float out_of_bounds_penalty = 0;

	//TODO: remove hydrogen atoms completely
	if (t > 1) { //ignore hydrogens
		//now consider interaction with every possible receptor atom
		//TODO: parallelize
		atom_params lin = ligs[l];
		float3 xyz = lin.coords;

		//evaluate for out of boundsness
		for (unsigned i = 0; i < 3; i++) {
			float min = dinfo.gridbegins[i];
			float max = dinfo.gridends[i];
			if (xyz[i] < min) {
				out_of_bounds_deriv[i] = -1;
				out_of_bounds_penalty += fabs(min - xyz[i]);
				xyz[i] = min;
			} else if (xyz[i] > max) {
				out_of_bounds_deriv[i] = 1;
				out_of_bounds_penalty += fabs(max - xyz[i]);
				xyz[i] = max;
			}
			out_of_bounds_deriv[i] *= slope;
		}

		out_of_bounds_penalty *= slope;

		//compute squared difference
		atom_params rin = dinfo.rec_atoms[ridx];
		float3 diff = xyz - rin.coords;

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
	float this_e = block_sum<float>(rec_energy);
	float3 deriv = block_sum<float3>(rec_deriv);
	if (threadIdx.x == 0) {
		curl(this_e, (float *) &deriv, v);
		out[l] += force_energy_tup(deriv + out_of_bounds_deriv,
				this_e + out_of_bounds_penalty);
	}
}

__global__ void reduce_energy(force_energy_tup *result, int n) {
	int l = threadIdx.x;
	float val = 0;
	if (l < n)
		val = result[l].energy;
	float e = block_sum<float>(val);
	if (l == 0) {
		result[0].energy = e;
	}
}

/* global */
/* void apply_penalties(float3 *coords, float3 gridBegins, float3 gridEnds, */
/*                      float3 *e_penalties, float3 *deriv_penalties){ */

/* } */

//host side of single point_calculation, energies and coords should already be initialized
float single_point_calc(const GPUNonCacheInfo *info, atom_params *ligs,
		force_energy_tup *out, float slope, unsigned nlig_atoms,
		unsigned nrec_atoms, float v) {
	/* Assumed by warp_sum */
	assert(THREADS_PER_BLOCK <= 1024);

	//this will calculate the per-atom energies and forces.
	//there is one execution stream for the blocks with
	//a full complement of threads and a separate stream
	//for the blocks that have the remaining threads
	unsigned nfull_blocks = nrec_atoms / THREADS_PER_BLOCK;
	unsigned nthreads_remain = nrec_atoms % THREADS_PER_BLOCK;

	if (nfull_blocks)
		interaction_energy<0> <<<dim3(nlig_atoms, nfull_blocks),
				THREADS_PER_BLOCK>>>(*info, 0, slope, v, ligs, out);
	if (nthreads_remain)
		interaction_energy<1> <<<nlig_atoms, nthreads_remain>>>(*info,
				nrec_atoms - nthreads_remain, slope, v, ligs, out);

	//get total energy
	reduce_energy<<<1, ROUND_TO_WARP(nlig_atoms)>>>(out, nlig_atoms);
	cudaThreadSynchronize();
	/* cudaStreamSynchronize(0); */
	abort_on_gpu_err();

	force_energy_tup ret;
	cudaMemcpy(&ret, out, sizeof(force_energy_tup), cudaMemcpyDeviceToHost);
	return ret.energy;
}

/* evaluate contribution of interacting pairs, add to forces and place total energy in e (which must be zero initialized) */
__global__
void eval_intra_kernel(const GPUSplineInfo * spinfo, const atom_params * atoms,
		const interacting_pair* pairs, unsigned npairs, float cutoff_sqr,
		float v, force_energy_tup *out, float *e) {

	float total = 0.0;
	for (unsigned i = 0; i < npairs; i += 1) {

		const interacting_pair& ip = pairs[i];
		float3 r = atoms[ip.b].coords - atoms[ip.a].coords;
		float r2 = dot(r, r);
		if (r2 < cutoff_sqr) {
			float3 deriv = float3(0, 0, 0);
			float dor;
			unsigned t1 = ip.t1;
			unsigned t2 = ip.t2;
			float energy = eval_deriv_gpu(spinfo, t1, atoms[ip.a].charge, t2,
					atoms[ip.b].charge, r2, dor);
			deriv = r * dor;
			curl(energy, (float*) &deriv, v);

			out[ip.b].minus_force.x += deriv.x;
			out[ip.b].minus_force.y += deriv.y;
			out[ip.b].minus_force.z += deriv.z;

			out[ip.a].minus_force.x -= deriv.x;
			out[ip.a].minus_force.y -= deriv.y;
			out[ip.a].minus_force.z -= deriv.z;

			total += energy;
		}
	}
	*e = total;
}

