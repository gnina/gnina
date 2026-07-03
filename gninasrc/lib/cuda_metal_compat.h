/*
 * cuda_metal_compat.h
 *
 * Compatibility header that stubs out CUDA types, keywords, and API calls
 * so that .cu files can be compiled as plain C++ (with LANGUAGE CXX in CMake)
 * when building the Metal backend (USE_METAL=ON) on Apple Silicon.
 *
 * Included automatically instead of <cuda_runtime.h> when USE_METAL is defined.
 * Real GPU work is done by Metal shaders (.metal files) dispatched from the host.
 */

#ifndef CUDA_METAL_COMPAT_H
#define CUDA_METAL_COMPAT_H

#ifdef USE_METAL

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>

// ─────────────────────────────────────────────────────────────────────────────
// CUDA keyword stubs
// All device annotations become no-ops: device functions are compiled as
// regular C++ and will be re-implemented as Metal shaders in later phases.
// ─────────────────────────────────────────────────────────────────────────────
#define __global__
#define __device__
#define __host__
#define __shared__
#define __forceinline__  inline
#define __align__(n)     alignas(n)

// Synchronisation: no-op in single-threaded host compilation
#define __syncthreads()   do {} while (0)
#define __threadfence()   do {} while (0)
#define __threadfence_block() do {} while (0)

// Thread/block index stubs — host code that queries these outside a kernel
// should not exist; these catch accidental usage.
struct _CudaDim3Stub { unsigned int x = 0, y = 0, z = 0; };
static inline _CudaDim3Stub _gnina_stub_idx() { return {}; }
#define threadIdx  (_gnina_stub_idx())
#define blockIdx   (_gnina_stub_idx())
#define blockDim   (_CudaDim3Stub{512,1,1})
#define gridDim    (_CudaDim3Stub{1,1,1})

// ─────────────────────────────────────────────────────────────────────────────
// CUDA vector types
// libmolgrid/common.h (fnachon Metal fork) also defines these.
// Use a shared sentinel macro so only whichever header is included first wins.
// ─────────────────────────────────────────────────────────────────────────────
#ifndef GNINA_CUDA_VECTOR_TYPES_DEFINED
#define GNINA_CUDA_VECTOR_TYPES_DEFINED
struct float3 {
    float x, y, z;
};
struct float4 {
    float x, y, z, w;
};
struct int2   { int x, y; };
struct int3   { int x, y, z; };
struct uint2  { unsigned int x, y; };
struct uint3  { unsigned int x, y, z; };
struct ulong3 { unsigned long x, y, z; };

inline float3 make_float3(float x, float y, float z)          { return {x, y, z}; }
inline float4 make_float4(float x, float y, float z, float w) { return {x, y, z, w}; }
inline uint3  make_uint3(unsigned x, unsigned y, unsigned z)   { return {x, y, z}; }
inline int3   make_int3(int x, int y, int z)                   { return {x, y, z}; }
#endif // GNINA_CUDA_VECTOR_TYPES_DEFINED

// dim3: CUDA grid/block dimension type (always needed by gnina code)
#ifndef CUDA_METAL_DIM3_DEFINED
#define CUDA_METAL_DIM3_DEFINED
struct dim3 {
    unsigned int x, y, z;
    constexpr dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1)
        : x(x), y(y), z(z) {}
};
#endif

// ─────────────────────────────────────────────────────────────────────────────
// CUDA error types
// ─────────────────────────────────────────────────────────────────────────────
typedef int cudaError_t;
typedef int cudaError;
static constexpr cudaError_t cudaSuccess = 0;
static constexpr cudaError_t cudaErrorMemoryAllocation = 2;

enum cudaMemcpyKind {
    cudaMemcpyHostToHost     = 0,
    cudaMemcpyHostToDevice   = 1,
    cudaMemcpyDeviceToHost   = 2,
    cudaMemcpyDeviceToDevice = 3,
    cudaMemcpyDefault        = 4
};

typedef void* cudaStream_t;
static constexpr cudaStream_t cudaStreamPerThread = nullptr;

inline const char* cudaGetErrorString(cudaError_t) { return "USE_METAL build"; }
inline cudaError_t cudaGetLastError()              { return cudaSuccess; }
inline cudaError_t cudaDeviceSynchronize()         { return cudaSuccess; }
inline cudaError_t cudaStreamSynchronize(cudaStream_t) { return cudaSuccess; }

// ─────────────────────────────────────────────────────────────────────────────
// CUDA device types and attributes
// ─────────────────────────────────────────────────────────────────────────────
enum cudaDeviceAttr {
    cudaDevAttrComputeMode = 0
};
enum cudaComputeMode {
    cudaComputeModeDefault    = 0,
    cudaComputeModeProhibited = 2
};
struct cudaDeviceProp {
    char   name[256];
    size_t totalGlobalMem;
    int    major, minor;
    int    multiProcessorCount;
    int    maxThreadsPerBlock;
    int    warpSize;
};
inline cudaError_t cudaGetDeviceProperties(cudaDeviceProp* prop, int /*device*/) {
    if (prop) {
        *prop = cudaDeviceProp{};
        prop->totalGlobalMem = 16ULL * 1024 * 1024 * 1024;
        prop->warpSize = 32;
        prop->maxThreadsPerBlock = 1024;
        prop->major = 9; prop->minor = 0;
    }
    return cudaSuccess;
}
inline cudaError_t cudaSetDevice(int) { return cudaSuccess; }
inline cudaError_t cudaGetDevice(int* d) { if (d) *d = 0; return cudaSuccess; }
inline cudaError_t cudaDeviceGetAttribute(int* val, cudaDeviceAttr, int) {
    if (val) *val = cudaComputeModeDefault;
    return cudaSuccess;
}
inline int cudaGetDeviceCount(int* n) { if (n) *n = 1; return cudaSuccess; }

inline cudaError_t cudaMemGetInfo(size_t* free, size_t* total) {
    // Will be overridden by device_buffer.cpp (Metal recommendedMaxWorkingSetSize)
    *free  = 8ULL * 1024 * 1024 * 1024; // 8 GB placeholder
    *total = 8ULL * 1024 * 1024 * 1024;
    return cudaSuccess;
}

// cudaMalloc / cudaFree are intentionally NOT stubbed here:
// they are replaced by device_buffer (Metal MTLBuffer arena) throughout gnina.
// Any remaining direct cudaMalloc calls must be converted to device_malloc().

// cudaMemcpy stubs: on Apple Silicon, MTLStorageModeShared means GPU and CPU
// share the same physical memory — a plain memcpy is equivalent.
inline cudaError_t cudaMemcpy(void* dst, const void* src, size_t n, cudaMemcpyKind) {
    memcpy(dst, src, n);
    return cudaSuccess;
}
inline cudaError_t cudaMemcpyAsync(void* dst, const void* src, size_t n,
                                    cudaMemcpyKind, cudaStream_t = nullptr) {
    memcpy(dst, src, n);
    return cudaSuccess;
}
inline cudaError_t cudaMemset(void* ptr, int val, size_t n) {
    memset(ptr, val, n);
    return cudaSuccess;
}
inline cudaError_t cudaMemsetAsync(void* ptr, int val, size_t n,
                                    cudaStream_t = nullptr) {
    memset(ptr, val, n);
    return cudaSuccess;
}

// ─────────────────────────────────────────────────────────────────────────────
// Pinned memory stubs — unnecessary on Apple Silicon (unified memory)
// ─────────────────────────────────────────────────────────────────────────────
enum { cudaHostAllocDefault = 0 };
inline cudaError_t cudaHostAlloc(void** ptr, size_t n, unsigned /*flags*/) {
    *ptr = malloc(n);
    return *ptr ? cudaSuccess : cudaErrorMemoryAllocation;
}
inline cudaError_t cudaFreeHost(void* ptr) {
    free(ptr);
    return cudaSuccess;
}

// ─────────────────────────────────────────────────────────────────────────────
// Atomic stubs (single-threaded host code only)
// Real atomics live in Metal shaders (.metal) compiled in Phases 2–7.
// ─────────────────────────────────────────────────────────────────────────────
inline float    atomicAdd(float*    addr, float    val) { float    old = *addr; *addr += val; return old; }
inline int      atomicAdd(int*      addr, int      val) { int      old = *addr; *addr += val; return old; }
inline unsigned atomicAdd(unsigned* addr, unsigned val) { unsigned old = *addr; *addr += val; return old; }
inline double   atomicAdd(double*   addr, double   val) { double   old = *addr; *addr += val; return old; }

inline unsigned long long atomicCAS(unsigned long long* addr,
                                     unsigned long long cmp,
                                     unsigned long long val) {
    unsigned long long old = *addr;
    if (old == cmp) *addr = val;
    return old;
}

// ─────────────────────────────────────────────────────────────────────────────
// Warp / SIMD shuffle stubs
// Used inside device functions that become regular C++ in host compilation.
// Real implementations are in Metal shaders using simd_shuffle_down().
// ─────────────────────────────────────────────────────────────────────────────
inline float __shfl_down_sync(unsigned /*mask*/, float val, unsigned /*delta*/,
                               int /*width*/ = 32) { return val; }
inline int   __shfl_down_sync(unsigned /*mask*/, int   val, unsigned /*delta*/,
                               int /*width*/ = 32) { return val; }

// ─────────────────────────────────────────────────────────────────────────────
// CUDA_ARCH guard macro: on Metal there is no device compilation pass.
// Code guarded by #ifndef __CUDA_ARCH__ runs as-is (host only).
// ─────────────────────────────────────────────────────────────────────────────
// __CUDA_ARCH__ is never defined in the Metal build, so existing
//   #ifndef __CUDA_ARCH__ ... #else ... #endif  guards already work correctly.

#endif // USE_METAL
#endif // CUDA_METAL_COMPAT_H
