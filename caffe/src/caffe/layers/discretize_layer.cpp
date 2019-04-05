#include <vector>

#include "caffe/filler.hpp"
#include "caffe/layers/discretize_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void DiscretizeLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  N_ = this->layer_param_.discretize_param().bins();
  min_value = this->layer_param_.discretize_param().min_value();
  max_value = this->layer_param_.discretize_param().max_value();
  CHECK_GT(N_, 0) << "DiscretizeLayer bins must be positive.";
  CHECK_GT(max_value,min_value) << "DiscretizeLayer max_value must be greater than min_value.";

  this->param_propagate_down_.resize(this->blobs_.size(), false);
}

template <typename Dtype>
void DiscretizeLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  // Figure out the dimensions
  vector<int> top_shape = bottom[0]->shape();
  top[0]->Reshape(top_shape);

}

template <typename Dtype>
void DiscretizeLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  const Dtype* bottom_data = bottom[0]->cpu_data();
  Dtype* top_data = top[0]->mutable_cpu_data();
  memset(top_data, 0, sizeof(Dtype)*top[0]->count());
  CHECK_EQ(bottom[0]->shape(0), top[0]->shape(0)) << "First dimensions of bottom and top in DiscretizeLayer do not match";
  unsigned M_ = bottom[0]->count();
  double span = max_value-min_value;
  for (int n = 0; n < M_; ++n) {
    double value = bottom_data[n];
    value -= min_value;
    if(value < 0) value = 0;
    value /= span;
    unsigned bin = value*N_;
    if(bin >= N_) bin = N_-1;
    top_data[n] = bin;
  }

}

#ifdef CPU_ONLY
STUB_GPU(DiscretizeLayer);
#endif

INSTANTIATE_CLASS(DiscretizeLayer);
REGISTER_LAYER_CLASS(Discretize);

}  // namespace caffe
