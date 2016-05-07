#ifndef __GPU_MEM_H
#define __GPU_MEM_H
 
#include <cstddef>
#include <cuda_runtime.h>
#include <assert.h>
#include <vector>

template <class T>
struct gpu_managed_alloc {
  typedef T value_type;
  gpu_managed_alloc(){};
  template <class U> gpu_managed_alloc(const gpu_managed_alloc<U>& other){}
  T* allocate(std::size_t n){
    T *ret;
    cudaMallocManaged(&ret, sizeof(T) * n);
    /* TODO: proper exceptions */
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

/* template <class T> */
/* class gvec : std::vector<T, gpu_managed_alloc<vec> >{ */
/*     __device__ __host__ */
/*     int &operator[](int a){ */
/*         /\* TODO: heh heh *\/ */
/*         return ((int *) _M_impl._M_start)[a]; */
/*     }; */
/* } */

#endif
