#include <vector>

#include "caffe/layers/affinity_loss_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void AffinityLossLayer<Dtype>::Reshape(
  const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  LossLayer<Dtype>::Reshape(bottom, top);
  CHECK_EQ(bottom[0]->count(1), bottom[1]->count(1))
      << "Inputs must have the same dimension.";
  diff_.ReshapeLike(*bottom[0]);
}

template <typename Dtype>
void AffinityLossLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top) {
  int count = bottom[0]->count();
  Dtype sum = 0.0;
  Dtype cutoff = this->layer_param_.affinity_loss_param().cutoff();
  const Dtype *labels = bottom[1]->cpu_data();
  const Dtype *preds = bottom[0]->cpu_data();
  Dtype *d = diff_.mutable_cpu_data();

  for(unsigned i = 0; i < count; i++) {
	 Dtype label = labels[i];
	 Dtype pred = preds[i];
	 if(label > 0) { //normal euclidean
		 Dtype diff = pred-label;
		 d[i] = diff;
		 sum += diff*diff;
	 } else if(label == 0 && pred > cutoff) { //hinge like
		 Dtype diff = pred-cutoff;
		 d[i] = diff;
		 sum += diff*diff;
	 } else { //if negative or less than cutoff, ignore
		 d[i] = 0;
	 }
  }

  Dtype loss = sum / bottom[0]->num() / Dtype(2);
  top[0]->mutable_cpu_data()[0] = loss;
}

template <typename Dtype>
void AffinityLossLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
  for (int i = 0; i < 2; ++i) {
    if (propagate_down[i]) {
      const Dtype sign = (i == 0) ? 1 : -1;
      const Dtype alpha = sign * top[0]->cpu_diff()[0] / bottom[i]->num();
      caffe_cpu_axpby(
          bottom[i]->count(),              // count
          alpha,                              // alpha
          diff_.cpu_data(),                   // a
          Dtype(0),                           // beta
          bottom[i]->mutable_cpu_diff());  // b
    }
  }
}


INSTANTIATE_CLASS(AffinityLossLayer);
REGISTER_LAYER_CLASS(AffinityLoss);

}  // namespace caffe
