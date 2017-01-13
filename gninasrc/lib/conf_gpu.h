/*

 GPU optimized versions for conf and change.

 */

#ifndef VINA_CONF_GPU_H
#define VINA_CONF_GPU_H

#include "conf.h"
#include "matrix.h"
#include "gpu_util.h"
#include "gpu_math.h"

/* change is a single GPU allocated array of floats.
 * The first six are the position and orientation and the remaining
 * are torsions.  The class itself is allocated on the CPU (at least for now)
 * and GPU code must be passed raw vector.
 */
struct change_gpu {
	float *change_values;
	int n; //on cpu, size of change_values is 6+torsions

	change_gpu(const change& src);

	change_gpu(const change_gpu& src);

    __host__ __device__
	change_gpu& operator=(const change_gpu& src);

	~change_gpu();

	__host__ __device__ void clear();

	void invert();

	__device__ float dot(const change_gpu& rhs) const;

	__device__ void sub(const change_gpu& rhs);

    __device__
	void minus_mat_vec_product(const flmat_gpu& m, change_gpu& out) const;

    __host__ __device__
	sz num_floats() const;

	void get_data(std::vector<float>& d) const;

	void set_data(std::vector<float>& d) const;

	void print() const;

    void* operator new(size_t count);
    void operator delete(void* ptr) noexcept;
};

union conf_info {
	struct {
		float position[3];
		float orientation[4];
		float torsions[];
	};
	float values[];
};

struct conf_gpu {

	conf_info *cinfo;
	int n; //on cpu, size of conf_values is 7+torsions to include x,y,z and quaternion

	conf_gpu(const conf& src);

	void set_cpu(conf& dst) const;

	conf_gpu(const conf_gpu& src);

    __host__ __device__
	conf_gpu& operator=(const conf_gpu& src);

	~conf_gpu();

	__device__ void increment(const change_gpu& c, fl factor);

	void get_data(std::vector<float>& d) const;

	void set_data(std::vector<float>& d) const;

	void print() const;

    void* operator new(size_t count);
    void operator delete(void* ptr) noexcept;

};

#endif
