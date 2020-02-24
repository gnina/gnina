// Cos neuron activation function layer.
// Adapted from TanH layer code

#include <vector>

#include "caffe/layers/trig_layer.hpp"

namespace caffe {

template <typename Dtype>
void TrigLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  func = this->layer_param_.trig_param().function();
}

template <typename Dtype>
void TrigLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  const Dtype* bottom_data = bottom[0]->cpu_data();
  Dtype* top_data = top[0]->mutable_cpu_data();
  const int count = bottom[0]->count();
  for (int i = 0; i < count; ++i) {
    if(func == TrigParameter_Function_COS)
      top_data[i] = cos(bottom_data[i]);
    else if(func == TrigParameter_Function_SIN)
      top_data[i] = sin(bottom_data[i]);
    else
      CHECK(0) << "Unsupported trig function";
  }
}

template <typename Dtype>
void TrigLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
  if (propagate_down[0]) {
    const Dtype* bottom_data = bottom[0]->cpu_data();
    const Dtype* top_diff = top[0]->cpu_diff();
    Dtype* bottom_diff = bottom[0]->mutable_cpu_diff();
    const int count = bottom[0]->count();
    Dtype x, deriv;
    for (int i = 0; i < count; ++i) {
      x = bottom_data[i];
      if(func == TrigParameter_Function_COS)
        deriv = -sin(x);
      else if(func == TrigParameter_Function_SIN)
        deriv = cos(x);
      else
        CHECK(0) << "Unsupported trig function";
      bottom_diff[i] = top_diff[i] * deriv;
    }
  }
}

#ifdef CPU_ONLY
STUB_GPU(TrigLayer);
#endif

INSTANTIATE_CLASS(TrigLayer);
REGISTER_LAYER_CLASS(Trig);

}  // namespace caffe
