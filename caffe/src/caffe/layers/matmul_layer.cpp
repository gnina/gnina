#include <vector>

#include "caffe/filler.hpp"
#include "caffe/layers/matmul_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void MatMulLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  transpose_ = this->layer_param_.matmul_param().transpose();
  Reshape(bottom, top);
}

template <typename Dtype>
void MatMulLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  CHECK_EQ(bottom.size(), 2);
  CHECK_EQ(top.size(), 1);

  const int row_axis = bottom[0]->CanonicalAxisIndex(
    this->layer_param_.matmul_param().axis());
  CHECK(row_axis >= 0);

  const int col_axis = row_axis + 1;
  CHECK(bottom[0]->num_axes() > col_axis);
  CHECK(bottom[1]->num_axes() > col_axis);

  CHECK_EQ(bottom[0]->count(0, row_axis), bottom[1]->count(0, row_axis));
  CHECK_EQ(bottom[0]->count(col_axis), bottom[1]->shape(row_axis));

  n_pairs = bottom[0]->count(0, row_axis);
  M_ = bottom[0]->shape(row_axis);
  K_ = bottom[0]->count(col_axis);
  N_ = bottom[1]->shape(col_axis);

  vector<int> top_shape = bottom[1]->shape();
  top_shape[row_axis] = M_;
  top[0]->Reshape(top_shape);
}

template <typename Dtype>
void MatMulLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top)
{
  const Dtype* A_data = bottom[0]->cpu_data();
  const Dtype* B_data = bottom[1]->cpu_data();
  Dtype* C_data = top[0]->mutable_cpu_data();

  for (int i = 0; i < n_pairs; ++i) // C_data = A_data * B_data
  {
    caffe_cpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, M_, N_, K_,
      (Dtype) 1., A_data + i*M_*K_, B_data + i*K_*N_,
      (Dtype) 0., C_data + i*M_*N_);
  }
}

template <typename Dtype>
void MatMulLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
  const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  const Dtype* C_diff = top[0]->cpu_diff();
  const Dtype* A_data = bottom[0]->cpu_data();
  const Dtype* B_data = bottom[1]->cpu_data();
  Dtype* A_diff = bottom[0]->mutable_cpu_diff();
  Dtype* B_diff = bottom[1]->mutable_cpu_diff();

  caffe_set<Dtype>(bottom[0]->count(), (Dtype) 0., A_diff);
  caffe_set<Dtype>(bottom[1]->count(), (Dtype) 0., B_diff);

  if (propagate_down[0])
  {
    for (int i = 0; i < n_pairs; ++i) // A_diff = C_diff * B_data.T
    {
      caffe_cpu_gemm<Dtype>(CblasNoTrans, CblasTrans, M_, K_, N_,
        (Dtype) 1., C_diff + i*M_*N_, B_data + i*K_*N_,
        (Dtype) 0., A_diff + i*M_*K_);
    }
  }
  if (propagate_down[1])
  {
    for (int i = 0; i < n_pairs; ++i) // B_diff = A_data.T * C_diff
    {
      caffe_cpu_gemm<Dtype>(CblasTrans, CblasNoTrans, K_, N_, M_,
        (Dtype) 1., A_data + i*M_*K_, C_diff + i*M_*N_,
        (Dtype) 0., B_diff + i*K_*N_);
    }
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
