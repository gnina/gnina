#include <vector>

#include "caffe/filler.hpp"
#include "caffe/layers/matmul_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void MatMulLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  const Dtype* A_data = bottom[0]->gpu_data();
  const Dtype* B_data = bottom[1]->gpu_data();
  Dtype* C_data = top[0]->mutable_gpu_data();

  caffe_gpu_gemm_batch<Dtype>(
    transpose_A ? CblasTrans : CblasNoTrans,
    transpose_B ? CblasTrans : CblasNoTrans,
    M, N, K,
    (Dtype) 1., A_data, B_data,
    (Dtype) 0., C_data,
    batch_size);
}

template <typename Dtype>
void MatMulLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
  const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  const Dtype* A_data = bottom[0]->gpu_data();
  const Dtype* B_data = bottom[1]->gpu_data();
  const Dtype* C_diff = top[0]->gpu_diff();
  Dtype* A_diff = bottom[0]->mutable_gpu_diff();
  Dtype* B_diff = bottom[1]->mutable_gpu_diff();

  if (propagate_down[0])
  {
    caffe_gpu_gemm_batch<Dtype>(
      CblasNoTrans,
      transpose_B ? CblasNoTrans : CblasTrans,
      M, K, N,
      (Dtype) 1., C_diff, B_data,
      (Dtype) 0., A_diff,
      batch_size);
  }
  if (propagate_down[1])
  {
    caffe_gpu_gemm_batch<Dtype>(
      transpose_A ? CblasNoTrans : CblasTrans,
      CblasNoTrans,
      K, N, M,
      (Dtype) 1., A_data, C_diff,
      (Dtype) 0., B_diff,
      batch_size);
  }
}

INSTANTIATE_LAYER_GPU_FUNCS(MatMulLayer);

}  // namespace caffe
