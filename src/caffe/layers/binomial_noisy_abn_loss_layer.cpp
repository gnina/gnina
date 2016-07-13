#include <algorithm>
#include <cmath>
#include <vector>

#include "caffe/layers/binomial_noisy_abn_loss_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void BinomialNoisyABNLossLayer<Dtype>::Reshape(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  LossLayer<Dtype>::Reshape(bottom, top);
  theta0 = this->layer_param_.binomial_noisy_abn_loss_param().theta0();
  theta1 = this->layer_param_.binomial_noisy_abn_loss_param().theta1();

  CHECK_LT(theta0, 1.0) << "theta0 must between 0 and 1";
  CHECK_LT(theta1, 1.0) << "theta1 must between 0 and 1";
  CHECK_GE(theta0, 0.0) << "theta0 must between 0 and 1";
  CHECK_GE(theta1, 0.0) << "theta1 must between 0 and 1";

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
  CHECK_EQ(numclasses,2) << "BinomialNoisyABNLossLayer requires exactly two classes";
}

//compute noise modified cross entropy
template <typename Dtype>
void BinomialNoisyABNLossLayer<Dtype>::Forward_cpu(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  const Dtype* bottom_data = bottom[0]->cpu_data();
  const Dtype* bottom_label = bottom[1]->cpu_data();
  int num = bottom[0]->num();
  int dim = bottom[0]->count() / bottom[0]->num();
  Dtype loss = 0;
  for (int i = 0; i < num; ++i) {
    int label = static_cast<int>(bottom_label[i]);
    //get probability of label = 1
    Dtype prob = bottom_data[i * dim +1];
    CHECK_LE(prob, 1.0) << "BinomialNoisyABNLossLayer not receiving probabilities.";
    CHECK_GE(prob, 0.0) << "BinomialNoisyABNLossLayer not receiving probabilities.";
    prob = theta0*(1-prob)+(1-theta1)*prob;
    if(label == 0) prob = 1-prob; //invert
    prob = std::max(prob, Dtype(kLOG_THRESHOLD));
    loss -= log(prob);
  }
  top[0]->mutable_cpu_data()[0] = loss / num;
}

template <typename Dtype>
void BinomialNoisyABNLossLayer<Dtype>::Backward_cpu(
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
      //get probability of being 1
      Dtype prob = bottom_data[i * dim + 1];
      Dtype diff = 0;
      if(label == 0) {
        if(theta1 == 0 && prob == 1)
         diff = 0;
        else
          diff = theta1*prob/( (1-theta0)*(1-prob) + theta1*prob );
      }
      else if(label == 1) {
        if(theta0 == 0 && prob == 0)
        diff = 1;
        else
        diff = (1-theta1)*prob/ ( (1-theta1)*prob+theta0*(1-prob));
      }
      else {
        LOG(FATAL) << this->type() << " Require binary labels.";
      }
      diff = diff - prob;
      bottom_diff[i * dim + 1 ] = -scale * diff;
      bottom_diff[i * dim ] = scale*diff;
    } //for
  } //if
}

INSTANTIATE_CLASS(BinomialNoisyABNLossLayer);
REGISTER_LAYER_CLASS(BinomialNoisyABNLoss);

}  // namespace caffe
