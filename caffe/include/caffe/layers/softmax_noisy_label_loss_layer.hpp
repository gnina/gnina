#ifndef CAFFE_SOFTMAX_NOISY_LABEL_LOSS_LAYER_HPP_
#define CAFFE_SOFTMAX_NOISY_LABEL_LOSS_LAYER_HPP_

#include <vector>
#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "caffe/layers/loss_layer.hpp"
#include "caffe/layers/softmax_layer.hpp"


namespace caffe {
	
template <typename Dtype>
class SoftmaxWithNoisyLabelLossLayer : public LossLayer<Dtype> {
 public:
  enum NoiseType {
    NoiseFree = 0,
    RandomNoise,
    ConfusingNoise,
    NumNoiseType
  };

 public:
  explicit SoftmaxWithNoisyLabelLossLayer(const LayerParameter& param)
      : LossLayer<Dtype>(param),
        softmax_layer_y_(new SoftmaxLayer<Dtype>(param)),
        softmax_layer_z_(new SoftmaxLayer<Dtype>(param)) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "SoftmaxWithNoisyLabelLoss"; }
  virtual inline int ExactNumBottomBlobs() const { return 3; }
  virtual inline int ExactNumTopBlobs() const { return -1; }
  // Can have multiple top blobs.
  // 1. The softmax loss
  // 2. p(y | x), i.e., p(clean_label | image)
  // 3. p(z | x), i.e., p(noise_type | image)
  // 4. p(y, z | y_tilde, x), i.e., p(clean_label, noise_type | noisy_label, image)
  virtual inline int MinNumTopBlobs() const { return 1; }
  virtual inline int MaxNumTopBlobs() const { return 4; }

 protected:
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  virtual inline bool AllowForceBackward(const int bottom_index) const {
    return bottom_index != 2;
  }

  virtual void BackProp(const Blob<Dtype>& prob,
      const Blob<Dtype>& true_prob, Dtype lr, Blob<Dtype>* diff);

  shared_ptr<SoftmaxLayer<Dtype> > softmax_layer_y_;
  shared_ptr<SoftmaxLayer<Dtype> > softmax_layer_z_;

  vector<Blob<Dtype>*> softmax_bottom_vec_y_;
  vector<Blob<Dtype>*> softmax_bottom_vec_z_;
  vector<Blob<Dtype>*> softmax_top_vec_y_;
  vector<Blob<Dtype>*> softmax_top_vec_z_;

  Blob<Dtype> prob_y_;
  Blob<Dtype> prob_z_;
  Blob<Dtype> posterior_;

  int M_; // batch size
  int N_; // number of classes

  Dtype lr_z_;
};

}
#endif
