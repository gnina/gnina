/*

 GPU optimized versions for conf and change.

 */

#ifndef VINA_CONF_GPU_H
#define VINA_CONF_GPU_H

#include "conf.h"
#include "matrix.h"
#include "gpu_util.h"
#include "gpu_math.h"
#include "device_buffer.h"

struct change_gpu {
	float *values;
	int n; //size of ligand change_values is 6+torsions; residue is just torsions
    int flex_offset; //index where the residues start

	change_gpu(const change& src, float_buffer& buffer);

	change_gpu(const change_gpu& src, float_buffer& buffer);

    change_gpu(const change_gpu& src) = default;

    __device__
	change_gpu& operator=(const change_gpu& src);

	__device__ void clear();

	void invert();

	__device__ float dot(const change_gpu& rhs) const;

	__device__ void sub(const change_gpu& rhs);

    __device__
	void minus_mat_vec_product(const flmat_gpu& m, change_gpu& out) const;

    __host__ __device__
	sz num_floats() const;
};

struct conf_gpu {

	float *values;
    int n; //size of ligand conf_values is 7+torsions; residue is just torsions
    int flex_offset; //index where the residues start

	conf_gpu(const conf& src, float_buffer& buffer);

	void set_cpu(conf& dst) const;

	conf_gpu(const conf_gpu& src, float_buffer& buffer);

    conf_gpu(const conf_gpu& src) = default;

    __host__ __device__
	conf_gpu& operator=(const conf_gpu& src);

	__device__ void increment(const change_gpu& c, fl factor, unsigned* lig_subtree_sizes);

	void get_data(std::vector<float>& d) const;

	void set_data(std::vector<float>& d) const;
};

#endif
