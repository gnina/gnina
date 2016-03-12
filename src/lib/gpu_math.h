#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <cuda_runtime.h>

struct gfloat3 : float3{
    __host__ __device__
    gfloat3(void){};
    __host__ __device__ 
    gfloat3(float x, float y, float z) : float3(make_float3(x,y,z)){};
    
    __host__ __device__
    float& operator[](int b) const{
        return ((float *) this)[b];
    };
};

#define float3 gfloat3

#ifdef __CUDACC__

__device__ __inline__ static
float3 __shfl_down(const float3 &a, int delta) {
    return float3(__shfl_down(a.x, delta),
                  __shfl_down(a.y, delta),
                  __shfl_down(a.z, delta));
}

#endif

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

//static inline float operator|(float3 a, int b) {
//	return *((float*)&a + b);
//}


#endif
