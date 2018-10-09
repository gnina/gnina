#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"

namespace caffe {

template <typename Dtype>
void FlexLSTMLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
}

template <typename Dtype>
void FlexLSTMLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
}

INSTANTIATE_LAYER_GPU_FUNCS(FlexLSTMLayer);

}  // namespace caffe
