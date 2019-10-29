/*
 * Interface for reusing a single memory buffer 
 * to avoid repeated calls to cudaMalloc
 */

#ifndef __DEVICE_BUFFER_H
#define __DEVICE_BUFFER_H
#include <cuda_runtime.h>
#include <assert.h>

size_t available_mem(size_t num_cpu_threads);

class device_buffer {
    char* begin; //pointer to the beginning of the buffer
    size_t capacity; //buffer size in bytes is capacity*sizeof(type),
    char* next_alloc;  //pointer to the beginning of the unused region

    void *copy_bytes(void *from, size_t n_bytes, cudaMemcpyKind kind);

    cudaError_t alloc_bytes(void** alloc, size_t n_bytes);

    bool has_space(size_t n_bytes);

  public:
    device_buffer();
    void init(size_t capacity);

    void resize(size_t n_bytes);
    template<typename T>
    T* copy(T* cpu_object, size_t n_requested, cudaMemcpyKind kind); //returns ptr to segment
    void reinitialize() {
      next_alloc = begin;
    }

    ~device_buffer();
    template<typename T>
    cudaError_t alloc(T** alloc, size_t n_bytes);
    template<typename T>
    cudaError_t dealloc(T* alloc);
};

extern thread_local device_buffer thread_buffer;

cudaError_t device_alloc_bytes(void **alloc, size_t n_bytes);

cudaError_t device_free(void *buf);

template<typename T>
cudaError_t device_malloc(T **alloc, size_t n_bytes) {
  return device_alloc_bytes((void **) alloc, n_bytes);
}

template<typename T>
T* device_buffer::copy(T* cpu_object, size_t n_requested, cudaMemcpyKind kind) {
  return (T*) copy_bytes(cpu_object, n_requested * sizeof(T), kind);
}

template<typename T>
cudaError_t device_buffer::alloc(T** alloc, size_t n_bytes) {
  return alloc_bytes((void **) alloc, n_bytes);
}

template<typename T>
cudaError_t device_buffer::dealloc(T* alloc) {
  assert(alloc >= begin && alloc < next_alloc);
  return cudaSuccess;
}

#endif
