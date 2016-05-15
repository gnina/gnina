#include "nngridder.h"

// GPU routines for nngridder



//return squared distance between pt and (x,y,z)
__device__ float sqDistance(float3 pt, float x, float y, float z)
{
	float ret;
	float tmp = pt.x-x;
	ret = tmp*tmp;
	tmp  = pt.y-y;
	ret += tmp*tmp;
	tmp = pt.z-z;
	ret += tmp*tmp;
	return ret;
}

//origin is grid origin
//dim is dimension of cubic grid
//resolution is grid resolution
//n is number of atoms
//coords are xyz coors
//gridindex is which grid they belong in
//radii are atom radii
//grids are the output and are assumed to be zeroed
template<bool Binary> __global__ void gpu_grid_set(float3 origin, int dim, float resolution, float rmult, int n, float3 *coords, short *gridindex, float *radii, float *grids) 
{
	//figure out what grid point we are 
	unsigned xi = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned yi = threadIdx.y + blockIdx.y*blockDim.y;
	unsigned zi = threadIdx.z + blockIdx.z*blockDim.z;
			
	if(xi >= dim || yi >= dim || zi >= dim)
		return;	//bail if we're off-grid, this should not be common
	
	unsigned gsize = dim*dim*dim;
	//compute x,y,z coordinate of grid point
	float x = xi*resolution+origin.x;
	float y = yi*resolution+origin.y;
	float z = zi*resolution+origin.z;
	
	//TODO: evaluate setting to zero here
	
	//iterate over all atoms
	for(unsigned i = 0; i < n; i++)
	{
		float3 coord = coords[i];
		short which = gridindex[i];
		if(which >= 0) { //because of hydrogens on ligands
			float r = radii[i];
			float rsq  = r*r;
			float d = sqDistance(coord, x,y,z);
			
			if(Binary)
			{
				if(d < rsq)
				{
					//set gridpoint to 1
					unsigned goffset = which*gsize;
					unsigned off = (xi*dim+yi)*dim+zi;
					//printf("%f,%f,%f %d,%d,%d  %d  %d %d\n",x,y,z, xi,yi,zi, which, goffset,off);
					grids[goffset+off] = 1.0;
				}
			}
			else
			{
				//for non binary we want a gaussian were 2 std occurs at the radius
				//after which which switch to a quadratic
				//the quadratic is to fit to have both the same value and first order
				//derivative at the cross over point and a value and derivative of zero
				//at 1.5*radius
				//TODO: figure if we can do the math without sqrt
				double dist = sqrt(d);
				if (dist < r * rmult)
				{
					unsigned goffset = which*gsize;
					unsigned off = (xi*dim+yi)*dim+zi;
					unsigned gpos = goffset+off;
					
					if (dist <= r)
					{
						//return gaussian
						double h = 0.5 * r;
						double ex = -dist * dist / (2 * h * h);
						grids[gpos] += exp(ex);
					}
					else //return quadratic
					{
						double h = 0.5 * r;
						double eval = 1.0 / (M_E * M_E); //e^(-2)
						double q = dist * dist * eval / (h * h) - 6.0 * eval * dist / h
								+ 9.0 * eval;
						grids[gpos] += q;
					}
				}
			}
		}
	}

}


void NNGridder::setAtomsGPU(unsigned natoms, float3 *coords, short *gridindex, float *radii, unsigned ngrids, float *grids)
{
	//each thread is responsible for a grid point location and will handle all atom types
	//each block is 8x8x8=512 threads
	float3 origin = make_float3(dims[0].begin, dims[1].begin, dims[2].begin);
	dim3 threads(8,8,8);
	unsigned dim = dims[0].n+1; //number of grid points
	unsigned blocksperside = ceil(dim/8.0);
	dim3 blocks(blocksperside,blocksperside,blocksperside);
	
	CUDA_CHECK(cudaMemset(grids, 0, ngrids*dim*dim*dim*sizeof(float))); //TODO: see if faster to do in kernel
	if(binary)
	{
		gpu_grid_set<true><<<blocks,threads>>>(origin, dim, resolution, radiusmultiple, natoms, coords, gridindex, radii, grids);
		CUDA_CHECK(cudaPeekAtLastError() );
	}
	else 
	{
		gpu_grid_set<false><<<blocks,threads>>>(origin, dim, resolution, radiusmultiple, natoms, coords, gridindex, radii, grids);
		CUDA_CHECK(cudaPeekAtLastError() );
	}
}
