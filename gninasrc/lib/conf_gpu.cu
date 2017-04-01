/*

   GPU optimized versions for conf and change.

*/

#include "conf_gpu.h"

#define CUDA_KERNEL_LOOP(i, n) \
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; \
       i < (n); \
       i += blockDim.x * gridDim.x)

__global__ void scalar_mult_kernel(float mult, const int n, float *vals) {
	CUDA_KERNEL_LOOP(index, n)
	{
		vals[index] *= mult;
	}
}

change_gpu::change_gpu(const change& src, float_buffer& buffer) :
		values(NULL), n(0) {
	std::vector<float> data;
    for (auto& ligchange : src.ligands) {
	    n += 6; //position + orientation
	    for (unsigned i = 0; i < 3; i++)
	    	data.push_back(ligchange.rigid.position[i]);

	    for (unsigned i = 0; i < 3; i++)
	    	data.push_back(ligchange.rigid.orientation[i]);

	    n += ligchange.torsions.size();
	    for (unsigned i = 0, nn = ligchange.torsions.size(); i < nn; i++) {
	    	data.push_back(ligchange.torsions[i]);
	    }
    }

    flex_offset = n;

    for (auto& reschange : src.flex) {
	    n += reschange.torsions.size();
	    for (unsigned j = 0, m = reschange.torsions.size(); j < m; j++) {
	    	data.push_back(reschange.torsions[j]);
	    }
    }
	//copy to buffer, leaving scratch space for dot
	assert(n == data.size());
    data.push_back(0);
    values = buffer.copy(&data[0], n+1, cudaMemcpyHostToDevice);
}

//allocate and copy
change_gpu::change_gpu(const change_gpu& src, float_buffer& buffer) :
		n(src.n), values(NULL), flex_offset(src.flex_offset) {
    values = buffer.copy(src.values, n+1, cudaMemcpyDeviceToDevice);
}

__device__ 
change_gpu& change_gpu::operator=(const change_gpu& src) {
    assert(values && n == src.n);
#ifndef __CUDA_ARCH__    
    CUDA_CHECK_GNINA(cudaMemcpy(values, src.values, sizeof(float) * n,
            cudaMemcpyDeviceToDevice));
#else
    memcpy(values, src.values, sizeof(float) * n);
#endif
    
	return *this;
}

//dkoes - zeros out all differences
__device__ void change_gpu::clear() {
	memset(values, 0, sizeof(float) * n);
}

//dkoes - multiply by -1
void change_gpu::invert() {
	scalar_mult_kernel<<<1, min(GNINA_CUDA_NUM_THREADS, n)>>>(-1.0, n,
			values);
}

//return dot product
__device__ float change_gpu::dot(const change_gpu& rhs) const {
    //TODO: n is no longer necessarily small
    if (threadIdx.x < WARPSIZE) {
	    int start = threadIdx.x;
	    float val = 0.0;
	    for(int i = start; i < n; i += WARPSIZE)
	    {
	    	val += values[i]*rhs.values[i];
	    }
	    //now warp reduce with shuffle

	    for(uint offset = WARPSIZE>>1; offset > 0; offset >>= 1)
	    	val += __shfl_down(val, offset);

	    if(start == 0)
	    	values[n] = val;
    }
    __syncthreads();
	return values[n];
}

//subtract rhs from this
__device__ void change_gpu::sub(const change_gpu& rhs) {
    int nthreads = blockDim.x < n ? blockDim.x : n;
	for (int index = threadIdx.x; index < n; index += nthreads)
		values[index] -= rhs.values[index];
}

__device__
void change_gpu::minus_mat_vec_product(const flmat_gpu& m, change_gpu& out) const {
    int idx = threadIdx.x;
    fl sum = 0;
    VINA_FOR(j,n)
        sum += m(m.index_permissive(idx,j)) * values[j];
    out.values[idx] = -sum;
}

__host__ __device__
sz change_gpu::num_floats() const {
	return n;
}

conf_gpu::conf_gpu(const conf& src, float_buffer& buffer) :
		values(NULL), n(0) {
	std::vector<float> data;
    for (auto& ligconf : src.ligands) {
	    n += 7; //position + orientation(qt)
	    for (unsigned i = 0; i < 3; i++)
	    	data.push_back(ligconf.rigid.position[i]);

	    data.push_back(ligconf.rigid.orientation.R_component_1());
	    data.push_back(ligconf.rigid.orientation.R_component_2());
	    data.push_back(ligconf.rigid.orientation.R_component_3());
	    data.push_back(ligconf.rigid.orientation.R_component_4());

	    n += ligconf.torsions.size();
	    for (unsigned i = 0, nn = ligconf.torsions.size(); i < nn; i++) {
	    	data.push_back(ligconf.torsions[i]);
	    }
    }

    flex_offset = n;

    for (auto& resconf : src.flex) {
	    n += resconf.torsions.size();
	    for (unsigned j = 0, m = resconf.torsions.size(); j < m; j++) {
	    	data.push_back(resconf.torsions[j]);
	    }
    }

	assert(n == data.size());
    values = buffer.copy(&data[0], n, cudaMemcpyHostToDevice);
}

//set cpu to gpu values, assumes correctly sized
void conf_gpu::set_cpu(conf& dst) const {
	std::vector<float> d;
	get_data(d);
    unsigned n = 0;
    //N.B. this assumes the torsions are in the same order
    //TODO: ugliness
    for (auto& lig : dst.ligands) {
        for (unsigned i = 0; i < 3; i++) {
            lig.rigid.position[i] = d[n];
            n++;
        }
        lig.rigid.orientation.x = d[n];
        n++;
        lig.rigid.orientation.y = d[n];
        n++;
        lig.rigid.orientation.z = d[n];
        n++;
        lig.rigid.orientation.w = d[n];
        n++;
        for (unsigned i = 0; i < lig.torsions.size(); i++) {
            lig.torsions[i] = d[n];
            n++;
        }
    }

    for (auto& res : dst.flex) {
        for (unsigned i = 0; i < res.torsions.size(); i++) {
            res.torsions[i] = d[n];
            n++;
        }
    }
}

//copy within buffer
conf_gpu::conf_gpu(const conf_gpu& src, float_buffer& buffer) :
		n(src.n), values(NULL), flex_offset(src.flex_offset) {
    values = buffer.copy(src.values, n, cudaMemcpyDeviceToDevice);
}

__host__ __device__
conf_gpu& conf_gpu::operator=(const conf_gpu& src) {
    assert(values && n == src.n);
#ifndef __CUDA_ARCH__
	CUDA_CHECK_GNINA(
			cudaMemcpy(values, src.values, sizeof(float) * n,
					cudaMemcpyDeviceToDevice));
#else
    memcpy(values, src.values, sizeof(float) * n);
#endif
    
	return *this;
}

__device__ void conf_gpu::increment(const change_gpu& c, fl factor, unsigned*
        lig_subtree_sizes) {
    unsigned idx = threadIdx.x;
    unsigned nroots = lig_subtree_sizes[0];
    //update rigid with early threads
    if (idx < nroots) {
        //position
        unsigned torsion_offset = idx == 0 ? 0 : lig_subtree_sizes[idx];
        unsigned conf_offset = torsion_offset + (7*idx);
        unsigned change_offset = conf_offset - idx;
        for (int i = 0; i < 3; i++) 
            values[conf_offset + i] += c.values[change_offset + i]*factor;

	    //rotation
	    qt orientation(values[conf_offset + 3],values[conf_offset +
                4],values[conf_offset + 5],values[conf_offset + 6]);
	    vec rotation(factor * c.values[change_offset + 3], factor *
                c.values[change_offset + 4], factor *
                c.values[change_offset + 5]);
	    quaternion_increment(orientation, rotation);
	    values[conf_offset + 3] = orientation.R_component_1();
	    values[conf_offset + 4] = orientation.R_component_2();
	    values[conf_offset + 5] = orientation.R_component_3();
	    values[conf_offset + 6] = orientation.R_component_4();
    }
	//torsions updated by everybody else, with indexing to avoid touching rigid again
    else if (idx < (n - (nroots*6))) {     //n-(7*nroots)+nroots 
        unsigned treeid = nroots;
        for (unsigned i = 1; i < nroots+1; i++) {
            if (idx < (lig_subtree_sizes[i] + nroots))
                treeid = i;
        }
        unsigned conf_offset = idx - nroots + (7 * treeid);
        unsigned change_offset = conf_offset - treeid;
	    values[conf_offset] += normalized_angle(factor*c.values[change_offset]);
	    normalize_angle(values[conf_offset]);
    }
}

void conf_gpu::get_data(std::vector<float>& d) const {
	d.resize(n);
	CUDA_CHECK_GNINA(
			cudaMemcpy(&d[0], values, n * sizeof(float),
					cudaMemcpyDeviceToHost));
    sync_and_errcheck();
}

void conf_gpu::set_data(std::vector<float>& d) const {
	CUDA_CHECK_GNINA(
			cudaMemcpy(values, &d[0], n * sizeof(float),
					cudaMemcpyHostToDevice));
    sync_and_errcheck();
}
