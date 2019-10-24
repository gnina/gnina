/*
 * Interface for reusing a single memory buffer 
 * to avoid repeated calls to cudaMalloc
 */
#include "device_buffer.h"
#include <cmath>
#include "gpu_util.h"
#include <cassert>
#include <boost/thread/thread.hpp>
#include <cuda.h>

#define align_down_pow2(n, size)            \
    ((decltype (n)) ((uintptr_t) (n) & ~((size) - 1)))

#define align_up_pow2(n, size)                                    \
    ((decltype (n)) align_down_pow2((uintptr_t) (n) + (size) - 1, size))

size_t available_mem(size_t num_cpu_threads) {
  size_t free, total;
  cudaError_t res;
  res = cudaMemGetInfo(&free, &total);
  if (res != cudaSuccess) {
    std::cerr << "cudaMemGetInfo returned status " << res << "\n";
    return 1;
  }
  return (free / num_cpu_threads) * .8;
}

thread_local device_buffer thread_buffer;

device_buffer::device_buffer()
    : capacity(0), begin(nullptr), next_alloc(nullptr) {
}

void device_buffer::init(size_t capacity) {
  this->capacity = capacity;
  CUDA_CHECK_GNINA(cudaMalloc(&begin, capacity));
  next_alloc = begin;
}

bool device_buffer::has_space(size_t n_bytes) {
  size_t bytes_used = next_alloc - begin; //in elements
  return capacity - bytes_used >= n_bytes;
}

cudaError_t device_alloc_bytes(void **alloc, size_t n_bytes) {
  // static here gives us lazy initialization and saves us from
  // constructing a thread_buffer for non-worker threads not included in
  // the available_mem() calculation.
  return thread_buffer.alloc((char **) alloc, n_bytes);
}

cudaError_t device_free(void *buf) {
  return thread_buffer.dealloc(buf);
}

void device_buffer::resize(size_t n_bytes) {
  assert(
      begin == next_alloc
          || !(std::cerr
              << "Device buffer only supports resize when buffer is empty.\n"));
  if (n_bytes > capacity) {
    CUDA_CHECK_GNINA(cudaFree(begin));
    CUDA_CHECK_GNINA(cudaMalloc(&begin, n_bytes));
    capacity = n_bytes;
  }
}

cudaError_t device_buffer::alloc_bytes(void** alloc, size_t n_bytes) {
  //N.B. you need to resize appropriately before starting to copy into the
  //buffer. This function avoids resizing so data structures in the buffer can use
  //pointers, and the only internal protection is the following assert.
  bool can_alloc = has_space(n_bytes);
  assert(can_alloc);
  if (!can_alloc) {
    std::cerr << "Alloc of " << n_bytes << " failed for buffer with "
        << capacity - (next_alloc - begin) << "/" << capacity
        << " bytes free.\n";
    abort();
  }
  *alloc = (void *) next_alloc;
  next_alloc = align_up_pow2(next_alloc + n_bytes, 128);
  return cudaSuccess;
}

void* device_buffer::copy_bytes(void* cpu_object, size_t n_bytes,
    cudaMemcpyKind kind) {
  assert(has_space(n_bytes));
  void *r;
  CUDA_CHECK_GNINA(alloc_bytes(&r, n_bytes));
  CUDA_CHECK_GNINA(definitelyPinnedMemcpy(r, cpu_object, n_bytes, kind));
  return r;
}

device_buffer::~device_buffer() {
  CUDA_CHECK_GNINA(cudaFree(begin));
}
