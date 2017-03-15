#include <algorithm>
#include <cmath>
#include <vector>

#include "caffe/layers/rank_loss_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void RankLossLayer<Dtype>::Reshape(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  LossLayer<Dtype>::Reshape(bottom, top);

  //check single labels
  const vector<int>& labels = bottom[1]->shape();
  for(unsigned i = 1, n = labels.size(); i < n; i++) {
    //skip batch size (i = 0)
    CHECK_EQ(labels[i], 1);
  }

  //must have only two classes
  const vector<int>& preds = bottom[0]->shape();
  int numclasses = 1;
  for(unsigned i = 1, n = preds.size(); i < n; i++) {
    //skip batch size (i = 0)
    numclasses *= preds[i];
  }
  CHECK_LE(numclasses,2) << "RankLossLayer requires one or two classes";
}

//compute pairwise rank loss
template <typename Dtype>
void RankLossLayer<Dtype>::Forward_cpu(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  const Dtype* bottom_data = bottom[0]->cpu_data();
  const Dtype* bottom_label = bottom[1]->cpu_data();
  int num = bottom[0]->num();
  int dim = bottom[0]->count() / bottom[0]->num();
  Dtype loss = 0;
  for (int i = 0; i < num; ++i) {
    int label = static_cast<int>(bottom_label[i]);
    //get probability of label = 1
    //dim should be 1 or 2
    Dtype prob = bottom_data[i * dim +1];
    //TODO
    abort();
  }
  top[0]->mutable_cpu_data()[0] = loss / num;
}

template <typename Dtype>
void RankLossLayer<Dtype>::Backward_cpu(
    const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
  if (propagate_down[1]) {
    LOG(FATAL) << this->type()
               << " Layer cannot backpropagate to label inputs.";
  }
  if (propagate_down[0]) {
    const Dtype* bottom_data = bottom[0]->cpu_data();
    const Dtype* bottom_label = bottom[1]->cpu_data();
    Dtype* bottom_diff = bottom[0]->mutable_cpu_diff();
    int num = bottom[0]->num();
    int dim = bottom[0]->count() / bottom[0]->num();
    caffe_set(bottom[0]->count(), Dtype(0), bottom_diff);
    const Dtype scale = top[0]->cpu_diff()[0] / num;

    for (int i = 0; i < num; ++i) {
      int label = static_cast<int>(bottom_label[i]);
      //TODO
      abort();
    } //for
  } //if
}

INSTANTIATE_CLASS(RankLossLayer);
REGISTER_LAYER_CLASS(RankLossLayer);

}  // namespace caffe
