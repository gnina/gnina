#ifndef GPU_UTIL_H
#define GPU_UTIL_H

#include <stdio.h>
#include <cstdlib>
#include <iostream>

__host__ __device__ static inline void abort_on_gpu_err(void) {
  cudaError err = cudaGetLastError();
  if (cudaSuccess != err) {
    printf("cudaCheckError() failed at %s:%i : %s\n",
    __FILE__, __LINE__, cudaGetErrorString(err));
    // exit(-1);
  }
}

__host__ __device__ static inline void sync_and_errcheck(void) {
  cudaError err = cudaDeviceSynchronize();
  if (cudaSuccess != err) {
    printf("cuda async error at %s:%i : %s\n",
    __FILE__, __LINE__, cudaGetErrorString(err));
    // exit(-1);
  }
}

#ifndef __CUDA_ARCH__
// CUDA: various checks for different function calls.
#define CUDA_CHECK_GNINA(condition) \
  /* Code block avoids redefinition of cudaError_t error */ \
  do { \
    cudaError_t error = condition; \
    if(error != cudaSuccess) {                                          \
        std::cerr << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(error); abort(); } \
  } while (0)
#else
// TODO: probably don't want to make API calls on the device.
#define CUDA_CHECK_GNINA(condition) condition
#endif

cudaError definitelyPinnedMemcpy(void* dst, const void *src, size_t n,
    cudaMemcpyKind k);

#define GNINA_CUDA_NUM_THREADS (512)
#define WARPSIZE (32)
#define CUDA_THREADS_PER_BLOCK (512)

// CUDA: number of blocks for N threads with nthreads per block
__host__ __device__
inline int CUDA_GET_BLOCKS(const int N, const int nthreads) {
  return (N + nthreads - 1) / nthreads;
}

//round N up to a multiple of 32
__host__ __device__ inline int ROUND_TO_WARP(int N) {
  if (N % 32) {
    return ((N / 32) + 1) * 32;
  }
  return N;
}

#endif
