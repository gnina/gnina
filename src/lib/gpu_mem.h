#ifndef __GPU_MEM_H
#define __GPU_MEM_H
 
#include <cstddef>
#include <cuda_runtime.h>
#include <assert.h>
template <class T>
struct gpu_managed_alloc {
  typedef T value_type;
  gpu_managed_alloc(){};
  template <class U> gpu_managed_alloc(const gpu_managed_alloc<U>& other){}
  T* allocate(std::size_t n){
    T *ret;
    cudaMallocManaged(&ret, sizeof(T) * n);
    assert(ret);
    return ret;
  };
  void deallocate(T* p, std::size_t n){
    cudaFree(p);
  };
};
template <class T, class U>
bool operator==(const gpu_managed_alloc<T>&, const gpu_managed_alloc<U>&){
  return true; }
template <class T, class U>
bool operator!=(const gpu_managed_alloc<T>&, const gpu_managed_alloc<U>&){
  return false; }

#endif
