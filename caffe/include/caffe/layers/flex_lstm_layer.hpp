#ifndef CAFFE_FLEX_LSTM_LAYER_HPP_
#define CAFFE_FLEX_LSTM_LAYER_HPP_

#include "caffe/layers/lstm_layer.hpp"
#include "gninasrc/lib/gpu_math.h"
/**
 * @brief Alternative LSTM layer implementation that supports custom data access
 *        patterns.
 */

namespace caffe {
enum AccessPattern {
  /* 0 */ strided_cube, /* Access data in cubic grid with fractional
                           filter-width stride */
  num_patterns
};

template <typename Dtype>
  class FlexLSTMLayer : public LSTMLayer<Dtype> {
  public:
    explicit FlexLSTMLayer(const LayerParameter& param) : LSTMLayer<Dtype>(param) {}
    virtual inline const char* type() const { return "FlexLSTM"; }

  protected:
    virtual void FillUnrolledNet(NetParameter* net_param) const;
    virtual void RecurrentInputShapes(vector<BlobShape>* shapes) const;
};

/**
 * @brief A helper for FlexLSTMLayer; it updates the contents of the current
 *        data blob during forward and backward.
 */
template <typename Dtype>
class LSTMDataGetterLayer : public Layer<Dtype> {
  public:
    explicit LSTMDataGetterLayer(const LayerParameter& param) : Layer<Dtype>(param) {}
    virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);

    virtual inline const char* type() const { return "LSTMDataGetter"; }
    virtual inline int ExactNumBottomBlobs() const { return 2; } //data, x
    virtual inline int ExactNumTopBlobs() const { return 1; } //x

  protected:
    virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
    virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

    AccessPattern pattern;
    unsigned num_timesteps;
    unsigned batch_size;
    unsigned ntypes;
    unsigned dim;
    unsigned subgrid_dim;
    unsigned cube_stride;
    unsigned example_size;
    unsigned hidden_dim;
    static unsigned current_timestep;
};

template <typename Dtype> unsigned LSTMDataGetterLayer<Dtype>::current_timestep = 0;

//alternatively could use enum as idx into array of function pointers, I
//*think* this will be cleaner though esp with device ptrs. 
//polymorphism doesn't really help here because we can't construct these
//ahead of time and then copy to the device
template <typename Dtype>
struct data_handler {
   __host__ __device__ virtual void GetData(const Dtype* src, Dtype* dest, 
       unsigned batch_size, unsigned ntypes, unsigned subgrid_dim, unsigned dim, 
       unsigned current_timestep, unsigned cube_stride, unsigned example_size) = 0;
   __host__ __device__ virtual void AccumulateDiff(const Dtype* src, Dtype* dest, 
       unsigned batch_size, unsigned ntypes, unsigned subgrid_dim, unsigned dim, 
       unsigned current_timestep, unsigned cube_stride, unsigned example_size) = 0;
   __host__ __device__ virtual ~data_handler() {}
};

template <typename Dtype>
struct strided_cube_data_handler : public data_handler<Dtype> {
  __host__ __device__ virtual void GetData(const Dtype* src, Dtype* dest,
       unsigned batch_size, unsigned ntypes, unsigned subgrid_dim, unsigned dim, 
       unsigned current_timestep, unsigned cube_stride, unsigned example_size) {
    unsigned overall_size = dim * dim * dim;
    unsigned factor = (((dim - subgrid_dim) / cube_stride) + 1);
    unsigned x_offset = ((current_timestep / (factor * factor)) % factor) * cube_stride;
    unsigned y_offset = ((current_timestep / factor) % factor) * cube_stride;
    unsigned z_offset = (current_timestep % factor) * cube_stride;
    //extract a single "timestep" corresponding to the correct stride
#ifndef __CUDA_ARCH__
    for (unsigned batch_idx=0; batch_idx < batch_size; ++batch_idx) {
      for (unsigned grid=0; grid < ntypes; ++grid) {
        for (unsigned i=0; i<subgrid_dim; ++i) {
          for (unsigned j=0; j<subgrid_dim; ++j) {
            for (unsigned k=0; k<subgrid_dim; ++k) {
              dest[(((batch_idx * ntypes + grid) * subgrid_dim + i) * subgrid_dim + j) * 
                subgrid_dim + k] = src[batch_idx * example_size + grid * overall_size + 
                x_offset * dim * dim + y_offset * dim + z_offset + 
                ((i * dim) + j) * dim + k];
            }
          }
        }
      }
    }
#else
    //we stop when we're out of the subgrid
    unsigned subgrid_size = ntypes * subgrid_dim * subgrid_dim * subgrid_dim;
    unsigned subgrid_count = subgrid_size * batch_size;
    CUDA_KERNEL_LOOP(tidx, subgrid_count) {
      //where in the grid is this index?
      unsigned idx = tidx;
      unsigned k = idx % subgrid_dim;
      idx /= subgrid_dim;
      unsigned j = idx % subgrid_dim;
      idx /= subgrid_dim;
      unsigned i = idx % subgrid_dim;
      idx /= subgrid_dim;
      unsigned grid = idx % ntypes;
      idx /= ntypes;
      unsigned batch_idx = idx % batch_size;
      //what overall index does that correspond to?
      dest[(((batch_idx * ntypes + grid) * subgrid_dim + i) * subgrid_dim + j) * 
        subgrid_dim + k] = src[batch_idx * example_size + grid * overall_size + 
        x_offset * dim * dim + y_offset * dim + z_offset + 
        ((i * dim) + j) * dim + k];
    }
#endif
}

 __host__ __device__ virtual void AccumulateDiff(const Dtype* src, Dtype* dest, 
       unsigned batch_size, unsigned ntypes, unsigned subgrid_dim, unsigned dim, 
       unsigned current_timestep, unsigned cube_stride, unsigned example_size) {
   unsigned overall_size = dim * dim * dim;
   unsigned factor = (((dim - subgrid_dim) / cube_stride) + 1);
   unsigned x_offset = ((current_timestep / (factor * factor)) % factor) * cube_stride;
   unsigned y_offset = ((current_timestep / factor) % factor) * cube_stride;
   unsigned z_offset = (current_timestep % factor) * cube_stride;
#ifndef __CUDA_ARCH__
   for (unsigned batch_idx=0; batch_idx < batch_size; ++batch_idx) {
     for (unsigned grid=0; grid < ntypes; ++grid) {
       for (unsigned i=0; i<subgrid_dim; ++i) {
         for (unsigned j=0; j<subgrid_dim; ++j) {
           for (unsigned k=0; k<subgrid_dim; ++k) {
             dest[batch_idx * example_size + grid * overall_size + 
                x_offset * dim * dim + y_offset * dim + z_offset + 
                ((i * dim) + j) * dim + k] += 
             src[(((batch_idx * ntypes + grid) * subgrid_dim + i) * subgrid_dim + j) * 
               subgrid_dim + k];
           }
         }
       }
     }
   }
#else
   unsigned subgrid_size = ntypes * subgrid_dim * subgrid_dim * subgrid_dim;
   unsigned subgrid_count = subgrid_size * batch_size;
   CUDA_KERNEL_LOOP(tidx, subgrid_count) {
     unsigned idx = tidx;
     unsigned k = idx % subgrid_dim;
     idx /= subgrid_dim;
     unsigned j = idx % subgrid_dim;
     idx /= subgrid_dim;
     unsigned i = idx % subgrid_dim;
     idx /= subgrid_dim;
     unsigned grid = idx % ntypes;
     idx /= ntypes;
     unsigned batch_idx = idx % batch_size;
     atomicAdd(&dest[batch_idx * example_size + grid * overall_size + 
        x_offset * dim * dim + y_offset * dim + z_offset + 
        ((i * dim) + j) * dim + k], 
     		src[(((batch_idx * ntypes + grid) * subgrid_dim + i) * subgrid_dim + j) * 
     		 subgrid_dim + k]);
   }
#endif
 }
  __host__ __device__ virtual ~strided_cube_data_handler() {}
};

template <typename Dtype>
void LSTMKernelWrapper(const int nthreads, const Dtype* src, Dtype* dest,
    AccessPattern pattern, unsigned batch_size, unsigned ntypes, unsigned
    subgrid_dim, unsigned dim, unsigned current_timestep, unsigned cube_stride,
    unsigned example_size);

} // namespace caffe

#endif // CAFFE_FLEX_LSTM_LAYER_HPP_
