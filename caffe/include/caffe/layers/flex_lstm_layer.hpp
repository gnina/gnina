#ifndef CAFFE_FLEX_LSTM_LAYER_HPP_
#define CAFFE_FLEX_LSTM_LAYER_HPP_

#include "caffe/layers/lstm_layer.hpp"
/**
 * @brief Extends the LSTM layer implementation to support custom data access
 *        patterns.
 */

namespace caffe {

template <typename Dtype>
  class Flex_LSTMLayer : public LSTMLayer <Dtype> {
  public:
    explicit Flex_LSTMLayer(const LayerParameter& param) : Layer<Dtype>(param) {}
    virtual inline const char* type() const { return "Flex_LSTMUnit"; }
  protected:
    virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);

    enum AccessPattern {
      /* 0 */ strided_cube, /* Access data in cubic grid with fractional
                               filter-width stride */
      num_patterns
    };

    AccessPattern pattern;
    template<AccessPattern apat> void Forward_cpu(const vector<Blob<Dtype>*>& bottom, 
        const vector<Blob<Dtype>*>& top, bool gpu);

    template<AccessPattern apat> void Forward_gpu(const vector<Blob<Dtype>*>& bottom, 
        const vector<Blob<Dtype>*>& top, bool gpu);

    template<AccessPattern apat> void Backward_cpu(const vector<Blob<Dtype>*>& bottom, 
        const vector<Blob<Dtype>*>& top, bool gpu);

    template<AccessPattern apat> void Backward_gpu(const vector<Blob<Dtype>*>& bottom, 
        const vector<Blob<Dtype>*>& top, bool gpu);
};

} // namespace caffe

#endif // CAFFE_FLEX_LSTM_LAYER_HPP_
