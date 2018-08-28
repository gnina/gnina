/*

 GPU optimized versions for conf and change.

 */

#include "conf_gpu.h"
#include "model.h"
#include "gpu_debug.h"

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

size_t change_gpu::idx_cpu2gpu(size_t cpu_node_idx, size_t offset_in_node,
    const gpu_data& d) {
  size_t gpu_node_idx = d.node_idx_cpu2gpu(cpu_node_idx);

  constexpr size_t extra_floats_per_lig_root = 5;
  size_t lig_roots_before_node = min(gpu_node_idx, d.nlig_roots);
  size_t gpu_flat_idx = gpu_node_idx
      + extra_floats_per_lig_root * lig_roots_before_node + offset_in_node;

  return gpu_flat_idx;
}

// CPU conf torsions are stored in dfs order, relative to the model's
// trees. GPU conf torsions are in bfs order.
size_t conf_gpu::idx_cpu2gpu(size_t cpu_node_idx, size_t offset_in_node,
    const gpu_data& d) {
  size_t gpu_node_idx = d.node_idx_cpu2gpu(cpu_node_idx);

  constexpr size_t extra_floats_per_lig_root = 6;
  size_t lig_roots_before_node = min(gpu_node_idx, d.nlig_roots);
  size_t gpu_flat_idx = gpu_node_idx
      + extra_floats_per_lig_root * lig_roots_before_node + offset_in_node;

  return gpu_flat_idx;
}

change_gpu::change_gpu(const change& src, const gpu_data& d,
    device_buffer& buffer)
    : n(src.num_floats()) {
  std::unique_ptr<fl[]> data(new fl[n]);

  for (sz i = 0; i < num_floats(); i++) {
    sz cpu_node_idx;
    sz offset_in_node;
    fl cpu_val = src.get_with_node_idx(i, &cpu_node_idx, &offset_in_node);
    assert(offset_in_node < 6);
    assert(cpu_node_idx < n);

    data[change_gpu::idx_cpu2gpu(cpu_node_idx, offset_in_node, d)] = cpu_val;
  }

  values = buffer.copy(data.get(), n, cudaMemcpyHostToDevice);
}

//allocate and copy
change_gpu::change_gpu(const change_gpu& src, device_buffer& buffer)
    : n(src.n), values(NULL) {
  values = buffer.copy(src.values, n + 1, cudaMemcpyDeviceToDevice);
}

__device__ change_gpu& change_gpu::operator=(const change_gpu& src) {
  assert(values && n == src.n);
#ifndef __CUDA_ARCH__
  CUDA_CHECK_GNINA(definitelyPinnedMemcpy(values, src.values, sizeof(float) * n,
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
  scalar_mult_kernel<<<1, min(GNINA_CUDA_NUM_THREADS, n)>>>(-1.0, n, values);
}

//return dot product
__device__ float change_gpu::dot(const change_gpu& rhs) const {
  __shared__ float out;
  //TODO: n is no longer necessarily small
  if (threadIdx.x < WARPSIZE) {
    int start = threadIdx.x;
    float val = 0.0;
    for (int i = start; i < n; i += WARPSIZE)
    {
      val += values[i] * rhs.values[i];
    }
    //now warp reduce with shuffle

    for (uint offset = WARPSIZE >> 1; offset > 0; offset >>= 1)
      val += shuffle_down(val, offset);

    if (start == 0) out = val;
  }
  __syncthreads();
  return out;
}

//subtract rhs from this
__device__ void change_gpu::sub(const change_gpu& rhs) {
  int nthreads = blockDim.x < n ? blockDim.x : n;
  for (int index = threadIdx.x; index < n; index += nthreads)
    values[index] -= rhs.values[index];
}

__device__
void change_gpu::minus_mat_vec_product(const flmat_gpu& m,
    change_gpu& out) const {
  int idx = threadIdx.x;
  fl sum = 0;
  VINA_FOR(j,n)
    sum += m(m.index_permissive(idx, j)) * values[j];
  out.values[idx] = -sum;
}

__host__  __device__ sz change_gpu::num_floats() const {
  return n;
}

inline
bool constructor_valid(const conf_gpu& gpu, const conf& src,
    const gpu_data& d) {
  conf test_dst = src;
  gpu.set_cpu(test_dst, d);
  sz n = src.num_floats();
  assert(n == test_dst.num_floats());
  assert(test_dst == src);
  return true;
}

conf_gpu::conf_gpu(const conf& src, const gpu_data& d, device_buffer& buffer)
    : n(src.num_floats()) {
  std::unique_ptr<fl[]> data(new fl[n]);

  for (sz i = 0; i < n; i++) {
    sz cpu_node_idx;
    sz offset_in_node;
    fl cpu_val = src.get_with_node_idx(i, &cpu_node_idx, &offset_in_node);
    assert(offset_in_node < 7);
    assert(cpu_node_idx < n);

    data[conf_gpu::idx_cpu2gpu(cpu_node_idx, offset_in_node, d)] = cpu_val;
  }

  values = buffer.copy(data.get(), n, cudaMemcpyHostToDevice);

  assert(constructor_valid(*this, src, d));
}

//set cpu to gpu values, assumes correctly sized
void conf_gpu::set_cpu(conf& dst, const gpu_data& d) const {
  std::vector<fl> data;
  get_data(data);

  for (sz i = 0; i < n; i++) {
    // TODO: need get_with_node_idx for node_idx, but need operator()
    // for writeable-ref.
    sz cpu_node_idx;
    sz offset_in_node;
    dst.get_with_node_idx(i, &cpu_node_idx, &offset_in_node);
    assert(offset_in_node < 7);
    assert(cpu_node_idx < n);

    dst.flat_index(i) = data[conf_gpu::idx_cpu2gpu(cpu_node_idx, offset_in_node,
        d)];
  }
}

//copy within buffer
conf_gpu::conf_gpu(const conf_gpu& src, device_buffer& buffer)
    : n(src.n), values(NULL) {
  values = buffer.copy(src.values, n, cudaMemcpyDeviceToDevice);
}

__host__  __device__ conf_gpu& conf_gpu::operator=(const conf_gpu& src) {
  assert(values && n == src.n);
#ifndef __CUDA_ARCH__
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(values, src.values, sizeof(float) * n,
          cudaMemcpyDeviceToDevice));
#else
  memcpy(values, src.values, sizeof(float) * n);
#endif

  return *this;
}

__device__ void conf_gpu::increment(const change_gpu& c, fl factor,
    gpu_data* gdata) {
  unsigned idx = threadIdx.x;
  tree_gpu& tree = *gdata->treegpu;
  unsigned lig_roots = tree.nlig_roots;
  //update rigid with early threads
  if (idx < lig_roots) {
    //position
    unsigned conf_offset = idx * 7;
    unsigned change_offset = idx * 6;
    for (int i = 0; i < 3; i++)
      values[conf_offset + i] += c.values[change_offset + i] * factor;

    //rotation
    qt orientation(values[conf_offset + 3], values[conf_offset + 4],
        values[conf_offset + 5], values[conf_offset + 6]);
    vec rotation(factor * c.values[change_offset + 3],
        factor * c.values[change_offset + 4],
        factor * c.values[change_offset + 5]);
    quaternion_increment(orientation, rotation);
    values[conf_offset + 3] = orientation.R_component_1();
    values[conf_offset + 4] = orientation.R_component_2();
    values[conf_offset + 5] = orientation.R_component_3();
    values[conf_offset + 6] = orientation.R_component_4();
  }
  //torsions updated by everybody else, with indexing to avoid touching rigid again
  else
    if (idx < tree.num_nodes) {
      unsigned conf_offset = idx + (6 * lig_roots);
      unsigned change_offset = idx + (5 * lig_roots);
      values[conf_offset] += normalized_angle(factor * c.values[change_offset]);
      normalize_angle(values[conf_offset]);
    }
}

void conf_gpu::get_data(std::vector<float>& d) const {
  d.resize(n);
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(&d[0], values, n * sizeof(float),
          cudaMemcpyDeviceToHost));
}

void conf_gpu::set_data(std::vector<float>& d) const {
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(values, &d[0], n * sizeof(float),
          cudaMemcpyHostToDevice));
}

__device__
void change_gpu::print() const {
  pretty_print_array(values, n, "change_gpu", "%f");
}

__device__
void conf_gpu::print() const {
  pretty_print_array(values, n, "conf_gpu", "%f");
}
