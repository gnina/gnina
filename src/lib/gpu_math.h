#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <cuda_runtime.h>
#include "gpu_util.h"

host device __inline__
float3 operator+(const float3 &a, const float3 &b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

host device __inline__
float3 operator+=(float3 &a, const float3 &b) {
	return a = a + b;
}

template<typename T>
host device  __inline__
float3 operator*(float3 a, T b) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

template<typename T>
host device __inline__
float3 operator*(T b, float3 a) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

//static inline float operator|(float3 a, int b) {
//	return *((float*)&a + b);
//}


#endif
