#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"

namespace caffe {

__device__ void flat_idx_to_multi(unsigned tidx, unsigned subgrid_dim, 
    unsigned batch_size, unsigned ntypes, unsigned& i, unsigned& j, unsigned& k,
    unsigned& grid, unsigned& batch_idx) {
  //we have a flat index (our thread index) and map this to a valid 3d
  //index that falls within the subgrid
  unsigned idx = tidx;
  k = idx % subgrid_dim;
  idx /= subgrid_dim;
  j = idx % subgrid_dim;
  idx /= subgrid_dim;
  i = idx % subgrid_dim;
  idx /= subgrid_dim;
  grid = idx % ntypes;
  idx /= ntypes;
  batch_idx = idx % batch_size;
}

template <typename Dtype>
__global__ void LSTMFlexForward(const int nthreads, const Dtype* src, Dtype* dest,
    AccessPattern pattern, unsigned batch_size, unsigned ntypes, unsigned
    subgrid_dim, unsigned dim, unsigned current_timestep, unsigned cube_stride,
    unsigned example_size) {
  switch(pattern) {
    case AccessPatterns::strided_cube:
      {
        //use the current_timestep to find the location of the first value in
        //the subcube we're going to use at this timestep; this is our starting
        //offset
        unsigned overall_size = dim * dim * dim;
        unsigned factor = (((dim - subgrid_dim) / cube_stride) + 1);
        unsigned x_offset = ((current_timestep / (factor * factor)) % factor) * cube_stride;
        unsigned y_offset = ((current_timestep / factor) % factor) * cube_stride;
        unsigned z_offset = (current_timestep % factor) * cube_stride;

        unsigned subgrid_count = batch_size * ntypes * subgrid_dim * subgrid_dim * subgrid_dim;
        CUDA_KERNEL_LOOP(tidx, subgrid_count) {
          //where in the grid is this index?
          unsigned i;
          unsigned j;
          unsigned k;
          unsigned grid;
          unsigned batch_idx;
          flat_idx_to_multi(tidx, subgrid_dim, batch_size, ntypes, i, j, k,
              grid, batch_idx);
          //what overall index does that correspond to?
          unsigned subgrid_idx = (((batch_idx * ntypes + grid) * subgrid_dim + i) * 
              subgrid_dim + j) * subgrid_dim + k;
          unsigned overall_idx = batch_idx * example_size + grid * overall_size +
              x_offset * dim * dim + y_offset * dim + z_offset + 
              ((i * dim) + j) * dim + k;
          dest[subgrid_idx] = src[overall_idx];
        }
        break;
      }
    default:
      {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i == 0)
          assert(pattern < AccessPatterns::num_patterns);
      }
  }
}

template <typename Dtype>
void LSTMKernelWrapper(const int cube_size, const Dtype* src, Dtype* dest,
    AccessPattern pattern, unsigned batch_size, unsigned ntypes, unsigned
    subgrid_dim, unsigned dim, unsigned current_timestep, unsigned cube_stride,
    unsigned example_size) {
  //Used for tests of LSTMFlexForward
  LSTMFlexForward<Dtype><<<CAFFE_GET_BLOCKS(cube_size),
    CAFFE_CUDA_NUM_THREADS>>>(cube_size, src, dest,
        pattern, batch_size, ntypes, subgrid_dim,
        dim, current_timestep, cube_stride, example_size);
}

template <typename Dtype>
__global__ void LSTMFlexBackward(const int nthreads, const Dtype* src, Dtype* dest,
    Dtype* total_diff, const Dtype* partial_diff, AccessPattern pattern, unsigned batch_size, 
    unsigned ntypes, unsigned subgrid_dim, unsigned dim, unsigned current_timestep,
    unsigned cube_stride, unsigned example_size) {
  switch(pattern) {
    case AccessPatterns::strided_cube:
      {
        unsigned overall_size = dim * dim * dim;
        //to be used for accumulating diff for current subcube blob
        unsigned factor = (((dim - subgrid_dim) / cube_stride) + 1);
        unsigned x_offset = ((current_timestep / (factor * factor)) % factor) * cube_stride;
        unsigned y_offset = ((current_timestep / factor) % factor) * cube_stride;
        unsigned z_offset = (current_timestep % factor) * cube_stride;

        //to be used to compute indices for current_x blob update (to be ready
        //for the previous timestep), if we aren't at the first timestep
        unsigned x_offset_prev; 
        unsigned y_offset_prev; 
        unsigned z_offset_prev; 
        if (current_timestep > 0) {
          x_offset_prev = x_offset - cube_stride;
          y_offset_prev = y_offset - cube_stride;
          z_offset_prev = z_offset - cube_stride;
        }

        unsigned subgrid_count = batch_size * ntypes * subgrid_dim * subgrid_dim * subgrid_dim;
        //we have a flat index (our thread index) and map this to a valid 3d
        //index that falls within the subgrid. 
        CUDA_KERNEL_LOOP(tidx, subgrid_count) {
          //where in the grid is this index?
          unsigned i;
          unsigned j;
          unsigned k;
          unsigned grid;
          unsigned batch_idx;
          flat_idx_to_multi(tidx, subgrid_dim, batch_size, ntypes, i, j, k,
              grid, batch_idx);
          //what overall index does that correspond to?
          unsigned subgrid_idx = (((batch_idx * ntypes + grid) * subgrid_dim + i) * 
              subgrid_dim + j) * subgrid_dim + k;
          unsigned overall_idx = batch_idx * example_size + grid * overall_size +
              x_offset * dim * dim + y_offset * dim + z_offset + 
              ((i * dim) + j) * dim + k;
          //accumulate diff
          atomicAdd(&total_diff[overall_idx], partial_diff[subgrid_idx]);
          if (current_timestep > 0) {
            //also update current data blob to be accurate for previous timestep
            unsigned overall_idx_prev = batch_idx * example_size + grid * overall_size +
                x_offset_prev * dim * dim + y_offset_prev * dim + z_offset_prev + 
                ((i * dim) + j) * dim + k;
            dest[subgrid_idx] = src[overall_idx_prev];
          }
        }
        break;
      }
    default:
      {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i == 0)
          assert(pattern < AccessPatterns::num_patterns);
      }
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
}

INSTANTIATE_LAYER_GPU_FUNCS(LSTMDataGetterLayer);

template void LSTMKernelWrapper<float>(const int nthreads, const float* src, float* dest,
    AccessPattern pattern, unsigned batch_size, unsigned ntypes, unsigned
    subgrid_dim, unsigned dim, unsigned current_timestep, unsigned cube_stride,
    unsigned example_size);

}  // namespace caffe
