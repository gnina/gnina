#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <cuda_runtime.h>

float3 operator+(float3 a, float3 b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

float3 operator+=(float3 &a, float3 b) {
	return *a = *a + b;
}

#endif
