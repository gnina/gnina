#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"

namespace caffe {

template <typename Dtype>
__global__ void LSTMFlexForward(const int nthreads, const Dtype* src, Dtype* dest,
    AccessPattern pattern, unsigned batch_size, unsigned ntypes, unsigned
    subgrid_dim, unsigned dim, unsigned current_timestep, unsigned cube_stride,
    unsigned example_size) {
  if (pattern == AccessPattern::strided_cube) {
    strided_cube_data_handler<Dtype> handler;
    handler.GetData(src, dest, batch_size, ntypes, 
        subgrid_dim, dim, current_timestep, cube_stride, example_size);
  }
}

template <typename Dtype>
__global__ void LSTMFlexBackward(const int nthreads, const Dtype* src, Dtype* dest,
    Dtype* total_diff, const Dtype* partial_diff, AccessPattern pattern, unsigned batch_size, 
    unsigned ntypes, unsigned subgrid_dim, unsigned dim, unsigned current_timestep,
    unsigned cube_stride, unsigned example_size) {
  if (pattern == AccessPattern::strided_cube) {
    strided_cube_data_handler<Dtype> handler;
    handler.GetData(src, dest, batch_size, ntypes, 
        subgrid_dim, dim, current_timestep, cube_stride, example_size);
    //also accumulate gradients for the current timestep in the right location
    handler.AccumulateDiff(partial_diff, total_diff, batch_size, 
        ntypes, subgrid_dim, dim, current_timestep, cube_stride, example_size);
  }
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  const int count = top[0]->count();
  const Dtype* src = bottom[0]->gpu_data();
  Dtype* dest = top[0]->mutable_gpu_data();
  unsigned subgrid_size = batch_size * ntypes * subgrid_dim * subgrid_dim * subgrid_dim;
  LSTMFlexForward<Dtype><<<CAFFE_GET_BLOCKS(subgrid_size),
    CAFFE_CUDA_NUM_THREADS>>>(count, src, dest, pattern, batch_size, ntypes,
        subgrid_dim, dim, current_timestep, cube_stride, example_size);
  CUDA_POST_KERNEL_CHECK;
  if (current_timestep != num_timesteps - 1)
    ++current_timestep;
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
  Dtype* total_diff = bottom[0]->mutable_gpu_diff();
  const Dtype* partial_diff = top[0]->gpu_diff();
  if (current_timestep == num_timesteps-1) {
    //TODO: this is a synchronous call, right?
    CUDA_CHECK(cudaMemset(total_diff, 0, bottom[0]->count()));
  }
  const int count = top[0]->count();
  const Dtype* src = bottom[0]->gpu_data();
  Dtype* dest = top[0]->mutable_gpu_data();
  unsigned subgrid_size = batch_size * ntypes * subgrid_dim * subgrid_dim * subgrid_dim;
  LSTMFlexBackward<Dtype><<<CAFFE_GET_BLOCKS(subgrid_size),
    CAFFE_CUDA_NUM_THREADS>>>(count, src, dest, total_diff, partial_diff, pattern, 
        batch_size, ntypes, subgrid_dim, dim, current_timestep, cube_stride,
        example_size);
  CUDA_POST_KERNEL_CHECK;
  if (current_timestep != 0)
    --current_timestep;
}

INSTANTIATE_LAYER_GPU_FUNCS(LSTMDataGetterLayer);

}  // namespace caffe
