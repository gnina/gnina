#ifndef CAFFE_FLEX_LSTM_LAYER_HPP_
#define CAFFE_FLEX_LSTM_LAYER_HPP_

#include "caffe/layers/lstm_layer.hpp"
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
    virtual inline int ExactNumBottomBlobs() const { return 3; } //data, seqcont, h
    virtual inline int ExactNumTopBlobs() const { return 2; } //x, h_conted

    virtual inline bool AllowForceBackward(const int bottom_index) const {
      // Can't propagate to sequence continuation indicators.
      return bottom_index != 1;
    }

  protected:
    template <typename AccessPattern> void GetData(const Dtype* src, Dtype* dest);
    template <typename AccessPattern> void AccumulateDiff(const Dtype* src, Dtype* dest);

    virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
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
    static unsigned current_timestep;
};

} // namespace caffe

#endif // CAFFE_FLEX_LSTM_LAYER_HPP_
