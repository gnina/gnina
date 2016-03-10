#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <cuda_runtime.h>

__host__ __device__ __inline__ float3 operator+(float3 a, float3 b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ __device__ __inline__ float3 operator+=(float3 &a, float3 b) {
	return a = a + b;
}

template<typename T>
__host__ __device__  __inline__ float3 operator*(float3 a, T b) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

template<typename T>
__host__ __device__ __inline__ float3 operator*(T b, float3 a) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

//static inline float operator|(float3 a, int b) {
//	return *((float*)&a + b);
//}

template float3 operator*<float>(float3 a, float b);
template float3 operator*<float>(float b, float3 a);

#endif
