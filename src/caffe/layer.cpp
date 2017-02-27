#include "caffe/layer.hpp"

namespace caffe {


template <typename Dtype>
void Layer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top, 
                const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom){}



INSTANTIATE_CLASS(Layer);

}  // namespace caffe
