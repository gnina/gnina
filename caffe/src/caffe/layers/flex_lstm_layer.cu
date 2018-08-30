#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"

namespace caffe {
// using Flex_LSTMLayer<typename Dtype>::AccessPattern;
// 
// template <typename Dtype>
// __device__ Dtype sigmoid(const Dtype x) {
  // return Dtype(1) / (Dtype(1) + exp(-x));
// }
// 
// template <typename Dtype>
// __device__ Dtype tanh(const Dtype x) {
  // return Dtype(2) * sigmoid(Dtype(2) * x) - Dtype(1);
// }
// 
// template <typename Dtype>
// __global__ void LSTMActsForward(const int nthreads, const int dim,
                                // const Dtype* X, Dtype* X_acts) {
  // CUDA_KERNEL_LOOP(index, nthreads) {
    // const int x_dim = 4 * dim;
    // const int d = index % x_dim;
    // if (d < 3 * dim) {
      // X_acts[index] = sigmoid(X[index]);
    // } else {
      // X_acts[index] = tanh(X[index]);
    // }
  // }
// }
// 
// template <typename Dtype>
// __global__ void LSTMUnitForward(const int nthreads, const int dim,
    // const Dtype* C_prev, const Dtype* X, const Dtype* cont,
    // Dtype* C, Dtype* H) {
  // CUDA_KERNEL_LOOP(index, nthreads) {
    // const int n = index / dim;
    // const int d = index % dim;
    // const Dtype* X_offset = X + 4 * dim * n;
    // const Dtype i = X_offset[d];
    // const Dtype f = X_offset[1 * dim + d];
    // const Dtype o = X_offset[2 * dim + d];
    // const Dtype g = X_offset[3 * dim + d];
    // const Dtype c_prev = C_prev[index];
    // const Dtype c = cont[n] * f * c_prev + i * g;
    // C[index] = c;
    // const Dtype tanh_c = tanh(c);
    // H[index] = o * tanh_c;
  // }
// }
// 
// template<typename Dtype>
// template<AccessPattern apat>
// void Flex_LSTMLayer::Forward_gpu {
  // rather than explicitly unrolling the net and doing separate kernel
  // launches for each dummy layer, we'll handle everything in a single kernel
  // launch and share data via warp shuffles and shared memory whenever
  // possible
  // const int count = top[1]->count();
  // const Dtype* C_prev = bottom[0]->gpu_data();
  // const Dtype* X = bottom[1]->gpu_data();
  // const Dtype* cont = bottom[2]->gpu_data();
  // Dtype* X_acts = X_acts_.mutable_gpu_data();
  // Dtype* C = top[0]->mutable_gpu_data();
  // Dtype* H = top[1]->mutable_gpu_data();
  // const int X_count = bottom[1]->count();
  // NOLINT_NEXT_LINE(whitespace/operators)
  // LSTMActsForward<Dtype><<<CAFFE_GET_BLOCKS(X_count), CAFFE_CUDA_NUM_THREADS>>>(
      // X_count, hidden_dim_, X, X_acts);
  // CUDA_POST_KERNEL_CHECK;
  // NOLINT_NEXT_LINE(whitespace/operators)
  // LSTMUnitForward<Dtype><<<CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS>>>(
      // count, hidden_dim_, C_prev, X_acts, cont, C, H);
  // CUDA_POST_KERNEL_CHECK;
// }
}  // namespace caffe

