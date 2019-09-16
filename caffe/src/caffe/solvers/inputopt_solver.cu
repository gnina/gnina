#include "caffe/sgd_solvers.hpp"
#define CUDA_NUM_THREADS 256
#define CUDA_NUM_BLOCKS 48

namespace caffe {

  template <typename Dtype>
  __global__ void threshold_blob(Dtype* offset_tblob, size_t blobsize, 
      Dtype threshold_value) {
    unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned nthreads = blockDim.x * gridDim.x;
    if (threshold_value) {
      Dtype clip_min = threshold_value > 0 ? -threshold_value : threshold_value;
      Dtype clip_max = threshold_value > 0 ? threshold_value : -threshold_value;
      for (size_t k=tidx; k<blobsize; k+=nthreads) {
        *(offset_tblob+k) = fmaxf(clip_min, fminf(*(offset_tblob+k), clip_max));
      }
    }
    else {
      for (size_t k=tidx; k<blobsize; k+=nthreads) {
        if (*(offset_tblob+k) < 0) {
          *(offset_tblob+k) = 0;
        }
      }
    }
  }

  template <typename Dtype>
  void InputOptSolver<Dtype>::DoThresholdGPU(Dtype* offset_tblob, size_t blobsize, 
      Dtype threshold_value) {
    unsigned block_multiple = blobsize / CUDA_NUM_THREADS;
    unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
    threshold_blob<<<nblocks, CUDA_NUM_THREADS>>>(offset_tblob, blobsize,
        threshold_value);
  }

  template void InputOptSolver<float>::DoThresholdGPU(float*, size_t, float);
  template void InputOptSolver<double>::DoThresholdGPU(double*, size_t, double);
}
