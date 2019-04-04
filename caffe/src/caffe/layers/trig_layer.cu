// Cos neuron activation function layer.
// Adapted from TanH layer code

#include <vector>

#include "caffe/layers/trig_layer.hpp"

namespace caffe {

template <typename Dtype>
__global__ void CosForward(const int n, const Dtype* in, Dtype* out) {
  CUDA_KERNEL_LOOP(index, n) {
    out[index] = cos(in[index]);
  }
}

template <typename Dtype>
__global__ void SinForward(const int n, const Dtype* in, Dtype* out) {
  CUDA_KERNEL_LOOP(index, n) {
    out[index] = sin(in[index]);
  }
}

template <typename Dtype>
void TrigLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  const Dtype* bottom_data = bottom[0]->gpu_data();
  Dtype* top_data = top[0]->mutable_gpu_data();
  const int count = bottom[0]->count();
  // NOLINT_NEXT_LINE(whitespace/operators)
  if(func == TrigParameter_Function_COS)
    CosForward<Dtype><<<CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS>>>(
      count, bottom_data, top_data);
  else if(func == TrigParameter_Function_SIN)
    SinForward<Dtype><<<CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS>>>(
          count, bottom_data, top_data);
  else
    CHECK(0) << "Unsupported trig function";
  CUDA_POST_KERNEL_CHECK;
}

template <typename Dtype>
__global__ void CosBackward(const int n, const Dtype* in_diff,
    const Dtype* in_data, Dtype* out_diff) {
  CUDA_KERNEL_LOOP(index, n) {
    Dtype x = in_data[index];
    out_diff[index] = in_diff[index] * (-sin(x));
  }
}

template <typename Dtype>
__global__ void SinBackward(const int n, const Dtype* in_diff,
    const Dtype* in_data, Dtype* out_diff) {
  CUDA_KERNEL_LOOP(index, n) {
    Dtype x = in_data[index];
    out_diff[index] = in_diff[index] * (cos(x));
  }
}


template <typename Dtype>
void TrigLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
  if (propagate_down[0]) {
    const Dtype* bottom_data = bottom[0]->gpu_data();
    const Dtype* top_diff = top[0]->gpu_diff();
    Dtype* bottom_diff = bottom[0]->mutable_gpu_diff();
    const int count = bottom[0]->count();
    // NOLINT_NEXT_LINE(whitespace/operators)
    if(func == TrigParameter_Function_COS)
      CosBackward<Dtype><<<CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS>>>(
        count, top_diff, bottom_data, bottom_diff);
    else if(func == TrigParameter_Function_SIN)
      SinBackward<Dtype><<<CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS>>>(
        count, top_diff, bottom_data, bottom_diff);
    else
        CHECK(0) << "Unsupported trig function";
    CUDA_POST_KERNEL_CHECK;
  }
}

INSTANTIATE_LAYER_GPU_FUNCS(TrigLayer);


}  // namespace caffe
