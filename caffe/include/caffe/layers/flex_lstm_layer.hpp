#ifndef CAFFE_FLEX_LSTM_LAYER_HPP_
#define CAFFE_FLEX_LSTM_LAYER_HPP_

#include "caffe/layers/lstm_layer.hpp"
#include "gninasrc/lib/gpu_math.h"

void test_strided_cube_datagetter();

namespace caffe {
  namespace AccessPatterns {
    enum pattern {
      /* 0 */ strided_cube, /* Access data in cubic grid with fractional
                               filter-width stride */
      num_patterns
    };
  }
  typedef AccessPatterns::pattern AccessPattern;

/**
 * @brief Alternative LSTM layer implementation that supports custom data access
 *        patterns. Like the standard Caffe LSTMLayer implementation, it is
 *        implemented by unrolling the LSTM graph through time, but it also
 *        provides support for flexible data access patterns that are
 *        impossible or resource-intensive to implement using the standard LSTM
 *        layer. An example use case is when timesteps may overlap each other,
 *        as in the "strided cube" access pattern. Rather than perform the data
 *        duplication necessary to use the standard LSTM layer, using this Flex
 *        LSTM layer has a lightweight memory footprint due to reuse of a
 *        shared buffer that is updated at each timestep by the helper
 *        LSTMDataGetter layer. Another example use case is when the data
 *        dimension changes at each timestep; again the LSTMDataGetter is
 *        responsible for managing the blob dimensions and performing
 *        interpolation as necessary. The Flex LSTM layer inserts
 *        LSTMDataGetter layers as needed when it unrolls the computation. 
 */
template <typename Dtype>
  class FlexLSTMLayer : public LSTMLayer<Dtype> {
  public:
    explicit FlexLSTMLayer(const LayerParameter& param) : LSTMLayer<Dtype>(param) {}
    virtual inline const char* type() const { return "FlexLSTM"; }
    virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top);
    virtual void Reshape(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top);

  protected:
    virtual void FillUnrolledNet(NetParameter* net_param) const;
    virtual void RecurrentInputShapes(vector<BlobShape>* shapes) const;
};

/**
 * @brief A helper for FlexLSTMLayer; it updates the contents of the current
 *        data blob during forward and backward and accumulates the diff (if
 *        required) in the correct location during backward. 
 */
template <typename Dtype>
class LSTMDataGetterLayer : public Layer<Dtype> {
  public:
    explicit LSTMDataGetterLayer(const LayerParameter& param) : Layer<Dtype>(param), 
    num_timesteps(0), dim(0), subgrid_dim(0), cube_stride(0), example_size(0), 
    current_timestep(0) {}
    virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);

    virtual inline const char* type() const { return "LSTMDataGetter"; }
    virtual inline int ExactNumBottomBlobs() const { return 1; } //x
    virtual inline int ExactNumTopBlobs() const { return 1; } //current_x

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
    unsigned current_timestep;

    friend void ::test_strided_cube_datagetter();
};

} // namespace caffe

#endif // CAFFE_FLEX_LSTM_LAYER_HPP_
