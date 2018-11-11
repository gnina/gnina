#include <vector>

#include "caffe/filler.hpp"
#include "caffe/layers/matmul_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void MatMulLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  transpose_A = this->layer_param_.matmul_param().transpose_a();
  transpose_B = this->layer_param_.matmul_param().transpose_b();
  Reshape(bottom, top);
}

template <typename Dtype>
void MatMulLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  CHECK_EQ(bottom.size(), 2);
  CHECK_EQ(top.size(), 1);

  CHECK(bottom[0]->num_axes() > 2);
  CHECK(bottom[1]->num_axes() > 2);

  CHECK_EQ(bottom[0]->shape(0), bottom[1]->shape(0));
  batch_size = bottom[0]->shape(0);

  const int A_rows = bottom[0]->shape(1);
  const int A_cols = bottom[0]->count(2);

  const int B_rows = bottom[1]->shape(1);
  const int B_cols = bottom[1]->count(2);

  if (transpose_A)
  {
    M = A_cols;
    K = A_rows;
  }
  else
  {
    M = A_rows;
    K = A_cols;
  }

  if (transpose_B)
  {
    CHECK_EQ(K, B_cols);
    N = B_rows;
  }
  else
  {
    CHECK_EQ(K, B_rows);
    N = B_cols;
  }

  vector<int> top_shape(3);
  top_shape[0] = batch_size;
  top_shape[1] = M;
  top_shape[2] = N;
  top[0]->Reshape(top_shape);
}

template <typename Dtype>
void MatMulLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  const Dtype* A_data = bottom[0]->cpu_data();
  const Dtype* B_data = bottom[1]->cpu_data();
  Dtype* C_data = top[0]->mutable_cpu_data();

  caffe_cpu_gemm_batch<Dtype>(
    transpose_A ? CblasTrans : CblasNoTrans,
    transpose_B ? CblasTrans : CblasNoTrans,
    M, N, K,
    (Dtype) 1., A_data, B_data,
    (Dtype) 0., C_data,
    batch_size);
}

template <typename Dtype>
void MatMulLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
  const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  const Dtype* A_data = bottom[0]->cpu_data();
  const Dtype* B_data = bottom[1]->cpu_data();
  const Dtype* C_diff = top[0]->cpu_diff();
  Dtype* A_diff = bottom[0]->mutable_cpu_diff();
  Dtype* B_diff = bottom[1]->mutable_cpu_diff();

  if (propagate_down[0])
  {
    caffe_cpu_gemm_batch<Dtype>(
      CblasNoTrans,
      transpose_B ? CblasNoTrans : CblasTrans,
      M, K, N,
      (Dtype) 1., C_diff, B_data,
      (Dtype) 0., A_diff,
      batch_size);
  }
  if (propagate_down[1])
  {
    caffe_cpu_gemm_batch<Dtype>(
      transpose_A ? CblasNoTrans : CblasTrans,
      CblasNoTrans,
      K, N, M,
      (Dtype) 1., A_data, C_diff,
      (Dtype) 0., B_diff,
      batch_size);
  }
}

template <typename Dtype>
void MatMulLayer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top,
  const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  //TODO
}

#ifdef CPU_ONLY
STUB_GPU(MatMulLayer);
#endif

INSTANTIATE_CLASS(MatMulLayer);
REGISTER_LAYER_CLASS(MatMul);

}  // namespace caffe
