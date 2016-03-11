/* dkoes
 * This file contains all the standalone gpu kernels.  There is (hopefully)
 * a nicer way to organize this, but I'm currently slightly flummoxed as to
 * how to cleaning mix object-oriented cpu and gpu code.
 */
#include "gpucode.h"
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <stdio.h>
#include "gpu_util.h"

#define THREADS_PER_BLOCK 1024
#define warpSize 32

global
void evaluate_splines(float **splines, float r, float fraction,
                      float cutoff, float *vals, float *derivs)
{
	unsigned i = blockIdx.x;
	float *spline = splines[i];
	vals[i] = 0;
	derivs[i] = 0;

	if (r >= cutoff || r < 0)
	{
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
device
float evaluate_spline(float *spline, float r, float fraction,
                      float cutoff, float& deriv)
{
	float val = 0;
	deriv = 0;
	if (r >= cutoff || r < 0)
	{
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

void evaluate_splines_host(const GPUSplineInfo& spInfo,
                           float r, float *device_vals, float *device_derivs)
{
    unsigned n = spInfo.n;
	evaluate_splines<<<n,1>>>((float**)spInfo.splines, r, spInfo.fraction, spInfo.cutoff,
                              device_vals, device_derivs);
}

device
float eval_deriv_gpu(GPUNonCacheInfo *dinfo, unsigned t,
                     float charge, unsigned rt, float rcharge, float r2, float& dor)
{
	float r = sqrt(r2);
	unsigned t1, t2;
	float charge1, charge2;
	if (t < rt)
	{
		t1 = t;
		t2 = rt;
		charge1 = fabs(charge);
		charge2 = fabs(rcharge);
	}
	else
	{
		t1 = rt;
		t2 = t;
		charge1 = fabs(rcharge);
		charge2 = fabs(charge);
	}

	unsigned tindex = t1 + t2 * (t2 + 1) / 2;
	GPUSplineInfo spInfo = dinfo->splineInfo[tindex];
	unsigned n = spInfo.n; //number of components

	float ret = 0, d = 0;

	//ick, hard code knowledge of components here; need to come up with
	//something mroe elegant
	//TypeDependentOnly,//no need to adjust by charge
	if (n > 0)
	{
		float fraction = spInfo.fraction;
		float cutoff = spInfo.cutoff;
		float val, deriv;
		val = evaluate_spline(spInfo.splines[0], r, fraction, cutoff, deriv);
		ret += val;
		d += deriv;
		//AbsAChargeDependent,//multiply by the absolute value of a's charge
		if (n > 1)
		{
			val = evaluate_spline(spInfo.splines[1], r, fraction, cutoff,
                                  deriv);
			ret += val * charge1;
			d += deriv * charge1;
			//AbsBChargeDependent,//multiply by abs(b)'s charge
			if (n > 2)
			{
				val = evaluate_spline(spInfo.splines[2], r, fraction, cutoff,
                                      deriv);
				ret += val * charge2;
				d += deriv * charge2;
				//ABChargeDependent,//multiply by a*b
				if (n > 3)
				{
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
device
void curl(float& e, float *deriv, float v)
{
	if (e > 0)
	{
		float tmp = (v / (v + e));
		e *= tmp;
		tmp *= tmp;
		for (unsigned i = 0; i < 3; i++)
			deriv[i] *= tmp;
	}
}

template <typename T> T __device__ __host__ zero(void);
template <> float3 zero(void){
    return float3(0,0,0);
}

template <> float zero(void){
    return 0;
}

//device functions for warp-based reduction using shufl operations
template <class T>
device __forceinline__
T warp_sum(T mySum) {
	for (int offset = warpSize>>1; offset > 0; offset>>=1)
        mySum += __shfl_down(mySum, offset);
	return mySum;
}

__device__ __forceinline__ 
bool isNotDiv32(unsigned int val) {
	return val & 31;
}

template <class T>
device __forceinline__
T block_sum(T* sdata, T mySum) {
	const unsigned int lane = threadIdx.x & 31;
	const unsigned int wid = threadIdx.x>>5;

	mySum = warp_sum(mySum);
	if (lane==0)
        sdata[wid] = mySum;
	__syncthreads();

	if (wid == 0) {
		mySum = (threadIdx.x < blockDim.x >> 5) ? sdata[lane] : zero<T>();
		mySum = warp_sum(mySum);
        if (threadIdx.x == 0 && isNotDiv32(blockDim.x))
            mySum += sdata[blockDim.x >> 5];
	}
	return mySum;
}

//calculates the energies of all ligand-prot interactions and combines the results
//into energies and minus forces
//needs enough shared memory for derivatives and energies of single ligand atom
//roffset specifies how far into the receptor atoms we are
template<bool remainder> global
void interaction_energy(GPUNonCacheInfo *dinfo, unsigned remainder_offset,
                        float slope, float v)
{
	unsigned l = blockIdx.x;
	unsigned r = blockDim.x - threadIdx.x - 1;
	unsigned roffset = remainder ? remainder_offset : blockIdx.y * THREADS_PER_BLOCK;
	unsigned ridx = roffset + r;
	//get ligand atom info
	unsigned t = dinfo->types[l];
	//TODO: remove hydrogen atoms completely
	if (t <= 1) //hydrogen ligand atom
		return;
    
	//now consider interaction with every possible receptor atom
	//TODO: parallelize
    float3 xyz = ((float3 *) dinfo->coords)[l];
    float3 out_of_bounds_deriv = float3(0, 0, 0);
    float out_of_bounds_penalty = 0;

    //evaluate for out of boundsness
    for (unsigned i = 0; i < 3; i++)
    {
        float min = dinfo->gridbegins[i];
        float max = dinfo->gridends[i];
        if (get(xyz, i) < min)
        {
            get(out_of_bounds_deriv, i) = -1;
            out_of_bounds_penalty += fabs(min - get(xyz, i));
            get(xyz, i) = min;
        }
        else if (get(xyz, i) > max)
        {
            get(out_of_bounds_deriv, i) = 1;
            out_of_bounds_penalty += fabs(max - get(xyz, i));
            get(xyz, i) = max;
        }
        get(out_of_bounds_deriv, i) *= slope;
    }

    out_of_bounds_penalty *= slope;


	//compute squared difference
	float rSq = 0;
	float3 diff = xyz - ((float3 *) dinfo->recoords)[ridx];
	for (unsigned j = 0; j < 3; j++)
	{
		float d = get(diff, j);
		get(diff, j) = d;
		rSq += d * d;
	}
	
	float rec_energy = 0;
	float3 rec_deriv = make_float3(0,0,0);
	if (rSq < dinfo->cutoff_sq)
	{
		//dkoes - the "derivative" value returned by eval_deriv
		//is normalized by r (dor = derivative over r?)
		float dor;
		rec_energy = eval_deriv_gpu(dinfo, t,
                                    dinfo->charges[l],
                                    dinfo->rectypes[ridx],
                                    dinfo->reccharges[ridx], rSq,
                                    dor);
		rec_deriv = diff * dor;
	}
    
    shared float energies[32];
	shared float3 derivs[32];
	float this_e = block_sum<float>(energies, rec_energy); 
	float3 deriv = block_sum<float3>(derivs, rec_deriv);
	if (threadIdx.x == 0)
	{
		curl(this_e, (float *) &deriv, v);
		
        ((float3 *) dinfo->minus_forces)[l] += deriv + out_of_bounds_deriv;
		dinfo->energies[l] += this_e + out_of_bounds_penalty;
	}
}

global void reduce_energy(GPUNonCacheInfo *dinfo) {
	int l  = threadIdx.x;
	shared float energies[warpSize];
	float my_energy = dinfo->energies[l];
	float e = block_sum<float>(energies, my_energy);
	if ( l == 0 ) {
		dinfo->energies[0] = e;
	}	
}


//host side of single point_calculation, energies and coords should already be initialized
float single_point_calc(GPUNonCacheInfo *dinfo, float *energies,
                        float slope, unsigned natoms,
                        unsigned nrecatoms, float v)
{
    /* Assumed by warp_sum */
    assert(THREADS_PER_BLOCK <= 1024);

	//this will calculate the per-atom energies and forces.
	//there is one execution stream for the blocks with
	//a full complement of threads and a separate stream
	//for the blocks that have the remaining threads
	unsigned nfull_blocks = nrecatoms / THREADS_PER_BLOCK;
	unsigned nthreads_remain = nrecatoms % THREADS_PER_BLOCK;
	
	dim3 fullGrid(natoms,nfull_blocks);
	dim3 remainGrid(natoms);

	if (nfull_blocks)
		interaction_energy<0><<<fullGrid,THREADS_PER_BLOCK,1>>>(dinfo, 0, slope, v);
	if (nthreads_remain) 
		interaction_energy<1><<<remainGrid,nthreads_remain>>>(dinfo, nrecatoms - nthreads_remain, slope, v);

	cudaThreadSynchronize();
	cudaError err = cudaGetLastError();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
				__FILE__, __LINE__, cudaGetErrorString(err));
		exit(-1);
	}
	//get total energy
	reduce_energy<<<1, natoms>>>(dinfo);
	cudaError err2 = cudaGetLastError();
	if (cudaSuccess != err2)
	{
		fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
				__FILE__, __LINE__, cudaGetErrorString(err2));
		exit(-1);
	}
	cudaThreadSynchronize();
	float e;
	cudaMemcpy(&e, &energies[0], sizeof(float), cudaMemcpyDeviceToHost);
	return e;
}
