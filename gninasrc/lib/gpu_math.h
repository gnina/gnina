/*
 * gpu_math.h
 *
 * GPU math types and utilities — backend-agnostic.
 *
 * Backends:
 *   USE_METAL=ON  → float3/float4 from cuda_metal_compat.h; device annotations
 *                   stripped to nothing; array3d_gpu uses memcpy (unified mem).
 *   USE_METAL=OFF → original CUDA implementation (cuda_runtime.h, __device__, …).
 */

#ifndef GPU_MATH_H
#define GPU_MATH_H

#include <float.h>
#include "array3d.h"
#include "common.h"
#include "gpu_util.h"
#include "device_buffer.h"

#ifdef USE_METAL
  #include "cuda_metal_compat.h"
  // float3, float4, make_float3, make_float4 come from cuda_metal_compat.h
#else
  #include <cuda_runtime.h>
  #include <random>
#endif

// ─────────────────────────────────────────────────────────────────────────────
// gfloat3 — wrapper over float3 with operator[] and arithmetic overloads.
// On CUDA: compiles as __host__ __device__ struct.
// On Metal: same plain C++ struct (annotations are empty macros).
// ─────────────────────────────────────────────────────────────────────────────
struct gfloat3 : float3 {
    gfloat3() = default;

    __host__ __device__ __inline__ gfloat3(float3 f) : float3(f) {}

    __host__ __device__ inline gfloat3(float x, float y, float z)
        : float3(make_float3(x, y, z)) {}

    __host__ __device__ explicit gfloat3(vec v)
        : float3(make_float3(v[0], v[1], v[2])) {}

    __host__ __device__
    float& operator[](int b) { return b == 0 ? x : b == 1 ? y : z; }

    __host__ __device__
    const float& operator[](int b) const { return b == 0 ? x : b == 1 ? y : z; }

    gfloat3& operator=(const gfloat3&) = default;

    __host__ __device__ inline float3& operator=(const vec& b) {
        x = b[0]; y = b[1]; z = b[2];
        return *this;
    }

    __host__ __device__ inline bool operator==(const float3& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    __host__ __device__ inline bool operator!=(const float3& rhs) const {
        return !(*this == rhs);
    }
};

inline std::ostream& operator<<(std::ostream& os, const gfloat3& f) {
    os << f.x << "," << f.y << "," << f.z;
    return os;
}

// ─────────────────────────────────────────────────────────────────────────────
// Warp / SIMD shuffle and atomic helpers
//
// CUDA path:  real __shfl_down_sync / atomicAdd intrinsics.
// Metal path: stubs for host compilation; real implementations go in .metal
//             shaders (simd_shuffle_down, atomic_fetch_add_explicit) — Phase 2+.
// ─────────────────────────────────────────────────────────────────────────────
#ifndef USE_METAL
#ifdef __CUDACC__

__device__ __inline__ float shuffle_down(float val, int offset) {
#if __CUDACC_VER_MAJOR__ >= 9
    return __shfl_down_sync(0xffffffff, val, offset);
#else
    return __shfl_down(val, offset);
#endif
}

__device__ __inline__ static gfloat3 shuffle_down(const gfloat3& a, int delta) {
#if __CUDACC_VER_MAJOR__ >= 9
    return gfloat3(__shfl_down_sync(0xffffffff, a.x, delta),
                   __shfl_down_sync(0xffffffff, a.y, delta),
                   __shfl_down_sync(0xffffffff, a.z, delta));
#else
    return gfloat3(__shfl_down(a.x, delta),
                   __shfl_down(a.y, delta),
                   __shfl_down(a.z, delta));
#endif
}

template<class T>
__device__ inline static T pseudoAtomicAdd(T* address, T value) {
    return T(atomicAdd(&((*address)[0]), value[0]),
             atomicAdd(&((*address)[1]), value[1]),
             atomicAdd(&((*address)[2]), value[2]),
             atomicAdd(&((*address)[3]), value[3]));
}

#endif // __CUDACC__

// Double-precision atomicAdd emulation for pre-Pascal CUDA GPUs.
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ inline double atomicAdd(double* address, double val) {
    unsigned long long int* addr_ull = (unsigned long long int*)address;
    unsigned long long int old = *addr_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(addr_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

#else // USE_METAL — host-only stubs (real impls live in .metal shaders, Phase 2+)

inline float shuffle_down(float val, int /*offset*/) { return val; }
inline gfloat3 shuffle_down(const gfloat3& a, int /*delta*/) { return a; }

template<class T>
inline T pseudoAtomicAdd(T* address, T value) {
    // Single-threaded stub: just add (safe on host, used for compilation only).
    (*address)[0] += value[0];
    (*address)[1] += value[1];
    (*address)[2] += value[2];
    (*address)[3] += value[3];
    return *address;
}

#endif // USE_METAL

// ─────────────────────────────────────────────────────────────────────────────
// Scalar helpers
// ─────────────────────────────────────────────────────────────────────────────
inline bool almostEqual(float a, float b) {
    float diff = std::fabs(a - b);
    if (a == b) return true;
    if (a == 0 || b == 0 || diff < FLT_MIN)
        return diff < (FLT_EPSILON * FLT_MIN);
    return diff / std::min(std::fabs(a) + std::fabs(b), FLT_MAX) < FLT_EPSILON;
}

__host__ __device__ inline static float dot(gfloat3 a, gfloat3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__host__ __device__ inline static float& get(gfloat3& a, int b) {
    return b == 0 ? a.x : b == 1 ? a.y : a.z;
}

__host__ __device__ inline static const float& get(const gfloat3& a, int b) {
    return b == 0 ? a.x : b == 1 ? a.y : a.z;
}

__host__ __device__ inline static gfloat3 operator-(const gfloat3& a) {
    return gfloat3(-a.x, -a.y, -a.z);
}

__host__ __device__ inline static gfloat3 operator+(const gfloat3& a,
                                                      const gfloat3& b) {
    return gfloat3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ __device__ inline static gfloat3 operator-(const gfloat3& a,
                                                      const gfloat3& b) {
    return gfloat3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ __device__ inline static gfloat3 operator+=(gfloat3& a,
                                                       const gfloat3& b) {
    return a = a + b;
}

__host__ __device__ inline static gfloat3 operator*(const gfloat3& a,
                                                      const gfloat3& b) {
    return gfloat3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

template<typename T>
__host__ __device__ inline static gfloat3 operator*(gfloat3 a, T b) {
    return gfloat3(a.x * b, a.y * b, a.z * b);
}

template<typename T>
__host__ __device__ inline static float3 operator*(T b, float3 a) {
    return gfloat3(a.x * b, a.y * b, a.z * b);
}

// ─────────────────────────────────────────────────────────────────────────────
// array3d_gpu — 3-D array stored in the device_buffer slab.
//
// On CUDA: the data pointer lives in device memory; __device__ accessors are
//          callable from GPU kernels.
// On Metal: the data pointer lives in an MTLStorageModeShared buffer whose
//           contents() is valid for both CPU and GPU.  The same memcpy path
//           works for upload.  __device__ annotations are stripped to nothing,
//           so the accessor methods are plain C++ calleable from host code.
//           Metal shader access to the same buffer is set up in Phase 4.
// ─────────────────────────────────────────────────────────────────────────────
template<typename T, typename U>
class array3d_gpu {
    sz i, j, k;
    T* data{};

public:
    array3d_gpu(const array3d<U>& carr)
        : i(carr.m_i), j(carr.m_j), k(carr.m_k) {
        CUDA_CHECK_GNINA(thread_buffer.alloc(&data, i * j * k * sizeof(T)));
        // definitelyPinnedMemcpy: memcpy on Metal, async DMA on CUDA.
        definitelyPinnedMemcpy(data, &carr.m_data[0],
                               sizeof(T) * carr.m_data.size(),
                               cudaMemcpyHostToDevice);
    }

    __device__ sz dim0() const { return i; }
    __device__ sz dim1() const { return j; }
    __device__ sz dim2() const { return k; }

    __device__ sz dim(sz idx) const {
        switch (idx) {
            case 0: return i;
            case 1: return j;
            case 2: return k;
            default: assert(false); return 0;
        }
    }

    __device__ T& operator()(sz ii, sz jj, sz kk) {
        return data[ii + i * (jj + j * kk)];
    }

    __device__ const T& operator()(sz ii, sz jj, sz kk) const {
        return data[ii + i * (jj + j * kk)];
    }
};

#endif // GPU_MATH_H
