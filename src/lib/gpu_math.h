#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <cuda_runtime.h>
#include "gpu_util.h"

#define float3(x, y, z) make_float3(x, y, z)

host device __inline__ static
float3 __shfl_down(const float3 &a, int delta) {
    return float3(__shfl_down(a.x, delta),
                  __shfl_down(a.y, delta),
                  __shfl_down(a.z, delta));
}


host device __inline__ static
float &get(float3 &a, int b){
    return
           b == 0 ? a.x :
           b == 1 ? a.y :
           a.z;
}

host device __inline__ static
float3 operator+(const float3 &a, const float3 &b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

host device __inline__ static
float3 operator-(const float3 &a, const float3 &b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

host device __inline__ static
float3 operator+=(float3 &a, const float3 &b) {
	return a = a + b;
}

template<typename T>
host device  __inline__ static
float3 operator*(float3 a, T b) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

template<typename T>
host device __inline__ static
float3 operator*(T b, float3 a) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

//static inline float operator|(float3 a, int b) {
//	return *((float*)&a + b);
//}


#endif
