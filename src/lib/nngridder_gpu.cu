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
__global__ void binary_set(float3 origin, int dim, float resolution, int n, float3 *coords, short *gridindex, float *radii, float *grids) 
{
	//figure out what grid point we are 
	unsigned xi = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned yi = threadIdx.y + blockIdx.y*blockDim.y;
	unsigned zi = threadIdx.z + blockIdx.y*blockDim.z;
			
	printf("%d %d %d %d\n",xi,yi,zi,dim);

	if(xi < dim || yi < dim || zi < dim)
		return;	//bail if we're off-grid, this should not be common
	
	unsigned gsize = dim*dim*dim;
	//compute x,y,z coordinate of grid point
	float x = xi*resolution+origin.x;
	float y = yi*resolution+origin.y;
	float z = zi*resolution+origin.z;
	
	//TODO: evaluate setting to zero here
	
	printf("%f %f %f\n",x,y,z);
	//iterate over all atoms
	for(unsigned i = 0; i < n; i++)
	{
		float3 coord = coords[i];
		short which = gridindex[i];
		if(which >= 0) { //because of hydrogens on ligands
			float r = radii[i];
			r *= r; //square radius
			float d = sqDistance(coord, x,y,z);
			if(d < r)
			{
				//set gridpoint to 1
				unsigned goffset = which*gsize;
				unsigned off = (xi*dim+yi)*dim+zi;
				printf("%d %d\n",goffset,off);
				grids[goffset+off] = 1.0;
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
		binary_set<<<blocks,threads>>>(origin, dim, resolution, natoms, coords, gridindex, radii, grids);
		CUDA_CHECK(cudaPeekAtLastError() );
	}
	else 
	{
		cerr << " non binary not gpu implemented yet" << "\n";
		exit(1);
	}
}
