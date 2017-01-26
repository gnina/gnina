#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <float.h>

#include <cuda_runtime.h>

#include "common.h"

/* This exists solely to provide constructor and [] operator
   funcs. Produced binaries are identical to those using vanilla
   float3. */
struct gfloat3 : float3{
    __host__ __device__ __inline__
    gfloat3() = default;
	__host__ __device__ __inline__
    gfloat3( float3 f): float3(f) {}
    __host__ __device__ __inline__ 
    gfloat3(float x, float y, float z) : float3(make_float3(x,y,z)){};
    
    __host__ __device__
    float& operator[](int b){
        return b == 0 ? x :
               b == 1 ? y :
               z;
    };
    
    __host__ __device__
    const float& operator[](int b) const{
        return b == 0 ? x :
               b == 1 ? y :
               z;
    };

    __host__ __device__
    gfloat3 &operator=(const gfloat3 &b) = default;

    
    __host__ __device__ __inline__
    float3 &operator=(const vec &b) {
        x = b[0];
        y = b[1];
        z = b[2];
        return *this;
    }

};

#define float3 gfloat3

//Both the shuffle and atomicAdd provided below are not strictly what they say
//they are. They are convenience functions that allow, for example, templated
//code to work correctly even if the types are unsupported by CUDA, but they
//work by applying hardware operations separately to individual elements of the
//respective types.

#ifdef __CUDACC__

__device__ __inline__ static
float3 __shfl_down(const float3 &a, int delta) {
    return float3(__shfl_down(a.x, delta),
                  __shfl_down(a.y, delta),
                  __shfl_down(a.z, delta));
}

template<class T>
__device__ __inline__ static
T pseudoAtomicAdd(T* address, T value) {
    return T(atomicAdd(&((*address)[0]), value[0]),
            atomicAdd(&((*address)[1]), value[1]),
            atomicAdd(&((*address)[2]), value[2]),
            atomicAdd(&((*address)[3]), value[3]));
}

#endif

inline bool almostEqual(float a, float b) {
    float absA = std::fabs(a);
    float absB = std::fabs(b);
    float diff = std::fabs(a-b);

    if (a == b) 
        return true;
    else if (a == 0 || b == 0 || diff < FLT_MIN) 
        return diff < (FLT_EPSILON * FLT_MIN);
    else 
        return diff / std::min((absA + absB), FLT_MAX) < FLT_EPSILON;
}


__host__ __device__ __inline__ static
float dot(float3 a, float3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

__host__ __device__ __inline__ static
float &get(float3 &a, int b){
    return
           b == 0 ? a.x :
           b == 1 ? a.y :
           a.z;
}

__host__ __device__ __inline__ static
const float &get(const float3 &a, int b){
    return
           b == 0 ? a.x :
           b == 1 ? a.y :
           a.z;
}

__host__ __device__ __inline__ static
float3 operator-(const float3 &a) {
	return float3(-a.x, -a.y, -a.z);
}

__host__ __device__ __inline__ static
float3 operator+(const float3 &a, const float3 &b) {
	return float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ __device__ __inline__ static
float3 operator-(const float3 &a, const float3 &b) {
	return float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ __device__ __inline__ static
float3 operator+=(float3 &a, const float3 &b) {
	return a = a + b;
}

template<typename T>
__host__ __device__  __inline__ static
float3 operator*(float3 a, T b) {
	return float3(a.x * b, a.y * b, a.z * b);
}

template<typename T>
__host__ __device__ __inline__ static
float3 operator*(T b, float3 a) {
	return float3(a.x * b, a.y * b, a.z * b);
}


#endif
