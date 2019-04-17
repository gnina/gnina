#include "../lib/gpu_math.h"
#define CUDA_NUM_THREADS 512
#define warpSize 32

// device functions for warp-based reduction using shufl operations
// TODO: should probably just be factored out into gpu_math or gpu_util
template<class T>
__device__   __forceinline__ T warp_sum(T mySum) {
  for (int offset = warpSize >> 1; offset > 0; offset >>= 1)
    mySum += shuffle_down(mySum, offset);
  return mySum;
}

__device__ __forceinline__
bool isNotDiv32(unsigned int val) {
  return val & 31;
}

/* requires blockDim.x <= 1024, blockDim.y == 1 */
template<class T>
__device__   __forceinline__ T block_sum(T mySum) {
  const unsigned int lane = threadIdx.x & 31;
  const unsigned int wid = threadIdx.x >> 5;

  __shared__ T scratch[32];

  mySum = warp_sum(mySum);
  if (lane == 0) scratch[wid] = mySum;
  __syncthreads();

  if (wid == 0) {
    mySum = (threadIdx.x < blockDim.x >> 5) ? scratch[lane] : 0;
    mySum = warp_sum(mySum);
    if (threadIdx.x == 0 && isNotDiv32(blockDim.x))
      mySum += scratch[blockDim.x >> 5];
  }
  return mySum;
}

__global__
void gpu_l2(const float* optgrid, const float* screengrid, float* scoregrid, 
    size_t M, size_t N, size_t gsize) {
  unsigned tidx = threadIdx.x;
  // optimized grids
  for (size_t i=0; i<M; ++i) {
    // conformers to screen against
    for (size_t j=0; j<N; ++j) {
      float sum = 0.;
      for (size_t k=0; k<gsize; k+=CUDA_NUM_THREADS) {
        float diff = optgrid[i * gsize + k] - screengrid[j * gsize + k];
        float sqdiff = diff * diff;
        sum += sqdiff;
      }
    float total = block_sum<float>(sum);
    if (tidx == 0)
      scoregrid[i * N + j] = sqrtf(total);
    }
  }
}
