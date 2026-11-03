/*
 * gpu_util.h
 *
 * GPU utility macros and functions — backend-agnostic.
 *
 * Backends:
 *   USE_METAL=ON  → no CUDA headers; error checking via abort(); memory
 *                   transfers are plain memcpy (unified memory on Apple Silicon).
 *   USE_METAL=OFF → original CUDA implementation unchanged.
 */

#ifndef GPU_UTIL_H
#define GPU_UTIL_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

// ─────────────────────────────────────────────────────────────────────────────
#ifdef USE_METAL
// ─────────────────────────────────────────────────────────────────────────────

#include "cuda_metal_compat.h"   // cudaMemcpyKind, cudaError_t, cudaSuccess, …

// No GPU error to check on the host side: Metal surfaces errors through
// MTLCommandBuffer.status. Stubs keep existing call-sites compiling.
static inline void abort_on_gpu_err() {}

// CUDA_CHECK_GNINA: on Metal just execute the expression (it returns
// cudaError_t from device_buffer stubs which are always cudaSuccess).
#define CUDA_CHECK_GNINA(condition) do { (void)(condition); } while (0)

// definitelyPinnedMemcpy: on Apple Silicon, Metal shared buffers and regular
// host pointers all live in the same unified physical memory.
// A plain memcpy is the correct and complete implementation.
inline cudaError_t definitelyPinnedMemcpy(void* dst, const void* src,
                                           size_t n, cudaMemcpyKind /*kind*/) {
    memcpy(dst, src, n);
    return cudaSuccess;
}

// ─────────────────────────────────────────────────────────────────────────────
#else   // CUDA backend (original)
// ─────────────────────────────────────────────────────────────────────────────

#include <cuda_runtime.h>

__host__ __device__ static inline void abort_on_gpu_err(void) {
    cudaError err = cudaGetLastError();
    if (cudaSuccess != err) {
        printf("cudaCheckError() failed at %s:%i : %s\n",
               __FILE__, __LINE__, cudaGetErrorString(err));
    }
}

#ifndef __CUDA_ARCH__
#define CUDA_CHECK_GNINA(condition) \
    do { \
        cudaError_t _err = (condition); \
        if (_err != cudaSuccess) { \
            std::cerr << __FILE__ << ":" << __LINE__ << ": " \
                      << cudaGetErrorString(_err); \
            std::abort(); \
        } \
    } while (0)
#else
#define CUDA_CHECK_GNINA(condition) condition
#endif

// Declared here; defined in gpucode.cu (uses pinned staging buffer + async copies).
cudaError definitelyPinnedMemcpy(void* dst, const void* src, size_t n,
                                  cudaMemcpyKind k);

#endif // USE_METAL

// ─────────────────────────────────────────────────────────────────────────────
// Constants and utility functions — identical for both backends.
// Annotations are stripped to nothing by cuda_metal_compat.h on Metal.
// ─────────────────────────────────────────────────────────────────────────────

#define GNINA_CUDA_NUM_THREADS  512
#define WARPSIZE                32
#define CUDA_THREADS_PER_BLOCK  512

// Ceiling division: number of blocks for N elements with nthreads per block.
__host__ __device__
inline int CUDA_GET_BLOCKS(const int N, const int nthreads) {
    return (N + nthreads - 1) / nthreads;
}

// Round N up to the nearest multiple of 32 (warp size / SIMD group size).
__host__ __device__ inline int ROUND_TO_WARP(int N) {
    return (N % 32) ? ((N / 32) + 1) * 32 : N;
}

#endif // GPU_UTIL_H
