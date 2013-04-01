/* dkoes
 * This file contains all the standalone gpu kernels.  There is (hopefully)
 * a nicer way to organize this, but I'm currently slightly flummoxed as to
 * how to cleaning mix object-oriented cpu and gpu code.
 */
#include "gpucode.h"
#include <thrust/reduce.h>
#include <stdio.h>

__global__ void evaluate_splines(float **splines, float r, float fraction,
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
__device__ float evaluate_spline(float *spline, float r, float fraction,
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
		float r, std::vector<float>& vals, std::vector<float>& derivs)
{
	unsigned n = spInfo.n;
	vals.resize(n);
	derivs.resize(n);

	float *device_vals, *device_derivs;
	cudaMalloc(&device_vals, sizeof(float) * n);
	cudaMalloc(&device_derivs, sizeof(float) * n);

	evaluate_splines<<<n,1>>>((float**)spInfo.splines, r, spInfo.fraction, spInfo.cutoff,
			device_vals, device_derivs);

	cudaMemcpy(&vals[0], device_vals, n * sizeof(float),
			cudaMemcpyDeviceToHost);
	cudaMemcpy(&derivs[0], device_derivs, n * sizeof(float),
			cudaMemcpyDeviceToHost);
	cudaFree(device_vals);
	cudaFree(device_derivs);

}

__device__ float eval_deriv_gpu(GPUNonCacheInfo *dinfo, unsigned t,
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
__device__ void curl(float& e, float *deriv, float v)
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

//calculates the energies of all ligand-prot interactions and combines the results
//into energies and minus forces
//needs enough shared memory for derivatives and energies of single ligand atom
//roffset specifies how far into the receptor atoms we are
__global__ void interaction_energy(GPUNonCacheInfo *dinfo, unsigned roffset,
		float slope, float v)
{
	unsigned l = blockIdx.x;
	unsigned r = threadIdx.x;
	unsigned ridx = roffset + r;
	//get ligand atom info
	unsigned t = dinfo->types[l];
	//TODO: remove hydrogen atoms completely
	if (t <= 1) //hydrogen ligand atom
		return;
	float out_of_bounds_deriv[3] =
			{ 0, 0, 0 };
	float out_of_bounds_penalty = 0;

	extern __shared__ float mysmem[];
	float *myenergies = mysmem;
	float *derivs = mysmem+blockDim.x;

	//initailize shared memory
	myenergies[r] = 0;
	derivs[3 * r] = 0;
	derivs[3 * r + 1] = 0;
	derivs[3 * r + 2] = 0;

	float charge = dinfo->charges[l];
	float xyz[3] =
			{ dinfo->coords[3 * l], dinfo->coords[3 * l + 1],
					dinfo->coords[3 * l + 2] };

	//evaluate for out of boundsness
	for (unsigned i = 0; i < 3; i++)
	{
		float min = dinfo->gridbegins[i];
		float max = dinfo->gridends[i];
		if (xyz[i] < min)
		{
			out_of_bounds_deriv[i] = -1;
			out_of_bounds_penalty += fabs(min - xyz[i]);
			xyz[i] = min;
		}
		else if (xyz[i] > max)
		{
			out_of_bounds_deriv[i] = 1;
			out_of_bounds_penalty += fabs(max - xyz[i]);
			xyz[i] = max;
		}
		out_of_bounds_deriv[i] *= slope;
	}

	out_of_bounds_penalty *= slope;

	//now consider interaction with every possible receptor atom
	//TODO: parallelize

	float cutoff = dinfo->cutoff_sq;

	unsigned rt;
	float rcharge;
	float rxyz[3];
	float diff[3];

	rt = dinfo->rectypes[ridx];
	rcharge = dinfo->reccharges[ridx];
	rxyz[0] = dinfo->recoords[3 * ridx];
	rxyz[1] = dinfo->recoords[3 * ridx + 1];
	rxyz[2] = dinfo->recoords[3 * ridx + 2];

	//compute squared difference
	float rSq = 0;
	for (unsigned j = 0; j < 3; j++)
	{
		float d = xyz[j] - rxyz[j];
		diff[j] = d;
		rSq += d * d;
	}

	if (rSq < cutoff)
	{
		//dkoes - the "derivative" value returned by eval_deriv
		//is normalized by r (dor = derivative over r?)
		float dor;
		myenergies[r] = eval_deriv_gpu(dinfo, t, charge, rt, rcharge, rSq,
				dor);
		for (unsigned j = 0; j < 3; j++)
		{
			derivs[3 * r + j] = dor * diff[j];
		}
	}

	__syncthreads();
	//horribly inefficient reduction; TODO improve
	if (r == 0)
	{
		float this_e = 0;
		float deriv[3] =
				{ 0, 0, 0 };
		unsigned nr = blockDim.x;
		for (unsigned i = 0; i < nr; i++)
		{
			this_e += myenergies[i];
			deriv[0] += derivs[3 * i];
			deriv[1] += derivs[3 * i + 1];
			deriv[2] += derivs[3 * i + 2];
		}
		curl(this_e, deriv, v);
		
		//update minus forces
		for (unsigned j = 0; j < 3; j++)
		{
			dinfo->minus_forces[3 * l + j] += deriv[j] + out_of_bounds_deriv[j];
		}
		//and energy
		dinfo->energies[l] += this_e + out_of_bounds_penalty;
	}
}
//calculates the energies of a single ligand atom (determined by block id)
__global__ void per_ligand_atom_energy(GPUNonCacheInfo *dinfo, float slope, float v)
{
        unsigned l = blockIdx.x;

        //get ligand atom info
        unsigned t = dinfo->types[l];

        //TODO: remove hydrogen atoms completely
        if (t <= 1)
                return; // hydrogen

        float charge = dinfo->charges[l];
        float xyz[3] =
        { dinfo->coords[3 * l], dinfo->coords[3 * l + 1],
                        dinfo->coords[3 * l + 2] };

        float out_of_bounds_deriv[3] =
        { 0, 0, 0 };
        float out_of_bounds_penalty = 0;

        //evaluate for out of boundsness
        for (unsigned i = 0; i < 3; i++)
        {
                float min = dinfo->gridbegins[i];
                float max = dinfo->gridends[i];
                if (xyz[i] < min)
                {
                        out_of_bounds_deriv[i] = -1;
                        out_of_bounds_penalty += fabs(min - xyz[i]);
                        xyz[i] = min;
                }
                else if (xyz[i] > max)
                {
                        out_of_bounds_deriv[i] = 1;
                        out_of_bounds_penalty += fabs(max - xyz[i]);
                        xyz[i] = max;
                }
                out_of_bounds_deriv[i] *= slope;
        }

        out_of_bounds_penalty *= slope;

        //now consider interaction with every possible receptor atom
        //TODO: parallelize
        
        float cutoff = dinfo->cutoff_sq;
        float this_e = 0;
        float deriv[3] = {0,0,0};
        unsigned n = dinfo->nrecatoms;
        unsigned rt;
        float rcharge;
        float rxyz[3];
        float diff[3];
        for(unsigned r = 0; r < n; r++) {
                rt = dinfo->rectypes[r];
                rcharge = dinfo->reccharges[r];
                rxyz[0] = dinfo->recoords[3*r];
                rxyz[1] = dinfo->recoords[3*r+1];
                rxyz[2] = dinfo->recoords[3*r+2];
                
                //compute squared difference
                float rSq = 0;
                for(unsigned j = 0; j < 3; j++) {
                        float d = xyz[j]-rxyz[j];
                        diff[j] = d;
                        rSq += d*d;
                }

                if(rSq < cutoff)
                {
                        //dkoes - the "derivative" value returned by eval_deriv
                        //is normalized by r (dor = derivative over r?)
                        float dor;
                        float e = eval_deriv_gpu(dinfo, t, charge, rt, rcharge, rSq, dor);
                        this_e += e;
                        for(unsigned j = 0; j < 3; j++) {
                                deriv[j] +=  dor * diff[j];
                        }
                }
        }

        curl(this_e, deriv, v);
        //update minus forces
        for(unsigned j = 0; j < 3; j++) {
                dinfo->minus_forces[3*l+j] = deriv[j]+out_of_bounds_deriv[j];
        }
        //and energy
        dinfo->energies[l] = this_e + out_of_bounds_penalty;
}


//host side of single point_calculation, energies and coords should already be initialized
float single_point_calc(GPUNonCacheInfo *dinfo, float *energies,
		float slope,
		unsigned natoms, unsigned nrecatoms, float v)
{
#if 1
	//this will calculate the per-atom energies and forces
#define THREADS_PER_BLOCK 1024
	for (unsigned off = 0; off < nrecatoms; off += THREADS_PER_BLOCK)
	{
		unsigned nr = nrecatoms - off;
		if (nr > THREADS_PER_BLOCK)
			nr = THREADS_PER_BLOCK;
		interaction_energy<<<natoms,nr, sizeof(float)*nr*4>>>(dinfo, off,slope, v);
		cudaError err = cudaGetLastError();
		if (cudaSuccess != err)
		{
			fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
					__FILE__, __LINE__, cudaGetErrorString(err));
			exit(-1);
		}
		cudaThreadSynchronize();
	}
#else
    //this will calculate the per-atom energies and forces
    per_ligand_atom_energy<<<natoms,1>>>(dinfo, slope, v);
#endif
	//get total energy
	thrust::device_ptr<float> dptr(energies);
	return thrust::reduce(dptr, dptr + natoms);
}