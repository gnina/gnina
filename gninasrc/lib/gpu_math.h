#ifndef GPU_MATH_H
#define GPU_MATH_H
#include <float.h>
#include "array3d.h"

#include <cuda_runtime.h>
#include "device_buffer.h"
#include <random>

#include "common.h"
#include "array3d.h"
#include "gpu_util.h"

/* This exists solely to provide constructor and [] operator
   funcs. Produced binaries are identical to those using vanilla
   float3. */
struct gfloat3 : float3{
    gfloat3() = default;
    __host__ __device__ __inline__
    gfloat3( float3 f): float3(f) {}
    __host__ __device__ inline 
    gfloat3(float x, float y, float z) : float3(make_float3(x,y,z)){};
    __host__ __device__ inline 
    gfloat3(vec v) : float3(make_float3(v[0],v[1],v[2])){};
    
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

    gfloat3 &operator=(const gfloat3 &b) = default;
    
    __host__ __device__ inline
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

__device__ __inline__ float shuffle_down(float val, int offset) {
  //wrapper for sync for normal flaot type
#if __CUDACC_VER_MAJOR__ >= 9
  return __shfl_down_sync(0xffffffff, val, offset);
#else
  return __shfl_down(val,offset);
#endif
}

__device__ __inline__ static
float3 shuffle_down(const float3 &a, int delta) {
#if __CUDACC_VER_MAJOR__ >= 9
    return float3(__shfl_down_sync(0xffffffff,a.x, delta),
                  __shfl_down_sync(0xffffffff,a.y, delta),
                  __shfl_down_sync(0xffffffff,a.z, delta));
#else
    return float3(__shfl_down(a.x, delta),
		  __shfl_down(a.y, delta),
		  __shfl_down(a.z, delta));
#endif
}

template<class T>
__device__ inline static
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


__host__ __device__ inline static
float dot(float3 a, float3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

__host__ __device__ inline static
float &get(float3 &a, int b){
    return
           b == 0 ? a.x :
           b == 1 ? a.y :
           a.z;
}

__host__ __device__ inline static
const float &get(const float3 &a, int b){
    return
           b == 0 ? a.x :
           b == 1 ? a.y :
           a.z;
}

__host__ __device__ inline static
float3 operator-(const float3 &a) {
	return float3(-a.x, -a.y, -a.z);
}

__host__ __device__ inline static
float3 operator+(const float3 &a, const float3 &b) {
	return float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__host__ __device__ inline static
float3 operator-(const float3 &a, const float3 &b) {
	return float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__host__ __device__ inline static
float3 operator+=(float3 &a, const float3 &b) {
	return a = a + b;
}

__host__ __device__ inline static
float3 operator*(const float3 &a, const float3 &b) {
    return make_float3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

template<typename T>
__host__ __device__  inline static
float3 operator*(float3 a, T b) {
	return float3(a.x * b, a.y * b, a.z * b);
}

template<typename T>
__host__ __device__ inline static
float3 operator*(T b, float3 a) {
	return float3(a.x * b, a.y * b, a.z * b);
}

template<typename T, typename U>
class array3d_gpu {
    cudaArray *data;
    cudaTextureObject_t tex;
    sz m_i, m_j, m_k;
public:
    array3d_gpu(const array3d<U>& carr) : m_i(carr.m_i), m_j(carr.m_j), m_k(carr.m_k) {
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
        cudaExtent volumeSize = make_cudaExtent(m_i, m_j, m_k);
        CUDA_CHECK_GNINA(cudaMalloc3DArray(&data, &channelDesc, volumeSize));

        // copy data to 3D array
        cudaMemcpy3DParms copyParams = {0};
        copyParams.srcPtr   = make_cudaPitchedPtr((void *)&carr.m_data[0], volumeSize.width*sizeof(T), volumeSize.width, volumeSize.height);
        copyParams.dstArray = data;
        copyParams.extent   = volumeSize;
        copyParams.kind     = cudaMemcpyHostToDevice;
        CUDA_CHECK_GNINA(cudaMemcpy3DAsync(&copyParams));

        // create texture object
        cudaResourceDesc resDesc;
        memset(&resDesc, 0, sizeof(resDesc));
        resDesc.resType = cudaResourceTypeLinear;
        resDesc.res.linear.devPtr = data;
        resDesc.res.linear.desc = channelDesc;
        resDesc.res.linear.sizeInBytes = m_i * m_j * m_k * sizeof(float);

        cudaTextureDesc texDesc;
        memset(&texDesc, 0, sizeof(texDesc));
        texDesc.readMode = cudaReadModeNormalizedFloat;
        texDesc.normalizedCoords = 1;

        cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL);
    }

	__device__ sz dim0() const { return m_i; }
	__device__ sz dim1() const { return m_j; }
	__device__ sz dim2() const { return m_k; }
	__device__ sz dim(sz i) const {
		switch(i) {
			case 0: return m_i;
			case 1: return m_j;
			case 2: return m_k;
			default: assert(false); return 0; // to get rid of the warning
		}
	}

	__device__ T& operator()(sz i, sz j, sz k)       { return m_data[i + m_i*(j + m_j*k)]; }
	__device__ const T& operator()(sz i, sz j, sz k) const { return m_data[i + m_i*(j + m_j*k)]; }
};

#endif
