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

struct gpu_data;

struct change_gpu {
    fl *values;
    int n; //size of ligand change_values is 6+torsions; residue is just torsions

    change_gpu(const change& src, const gpu_data& d, device_buffer& buffer);

    change_gpu(const change_gpu& src, device_buffer& buffer);

    change_gpu(const change_gpu& src) = default;

    __device__ change_gpu& operator=(const change_gpu& src);

    __device__ void clear();

    void invert();

    __device__ float dot(const change_gpu& rhs) const;

    __device__ void sub(const change_gpu& rhs);

    __device__
    void minus_mat_vec_product(const flmat_gpu& m, change_gpu& out) const;

    __host__  __device__ sz num_floats() const;

    __device__
    void print() const;
  private:
    static size_t idx_cpu2gpu(size_t cpu_node_idx, size_t offset_in_node,
        const gpu_data& d);

};

struct conf_gpu {

    float *values;
    int n; //size of ligand conf_values is 7+torsions; residue is just torsions

    conf_gpu(const conf& src, const gpu_data& d, device_buffer& buffer);

    void set_cpu(conf& dst, const gpu_data& d) const;

    conf_gpu(const conf_gpu& src, device_buffer& buffer);

    conf_gpu(const conf_gpu& src) = default;

    __host__  __device__ conf_gpu& operator=(const conf_gpu& src);

    __device__ void increment(const change_gpu& c, fl factor, gpu_data* gdata);

    void get_data(std::vector<float>& d) const;

    void set_data(std::vector<float>& d) const;

    __device__
    void print() const;
  private:
    static size_t idx_cpu2gpu(size_t cpu_node_idx, size_t offset_in_node,
        const gpu_data& d);
};

#endif
