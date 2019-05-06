#include "caffe/sgd_solvers.hpp"
#define CUDA_NUM_THREADS 256
#define CUDA_NUM_BLOCKS 48

namespace caffe {

  template <typename Dtype>
  __global__ void threshold_blob(Dtype* offset_tblob, size_t blobsize) {
    unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned nthreads = blockDim.x * gridDim.x;
    for (size_t k=tidx; k<blobsize; k+=nthreads) {
      if (*(offset_tblob+k) < 0) {
        *(offset_tblob+k) = 0;
      }
    }
  }

  template <typename Dtype>
  void InputOptSGDSolver<Dtype>::DoThresholdGPU(Dtype* offset_tblob, size_t blobsize) {
    unsigned block_multiple = blobsize / CUDA_NUM_THREADS;
    unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
    threshold_blob<<<nblocks, CUDA_NUM_THREADS>>>(offset_tblob, blobsize);
  }

  template void InputOptSGDSolver<float>::DoThresholdGPU(float*, size_t);
  template void InputOptSGDSolver<double>::DoThresholdGPU(double*, size_t);
}
