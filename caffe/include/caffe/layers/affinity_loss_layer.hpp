/*
 * affinity_loss_layer.hpp
 *
 *  Created on: Feb 5, 2017
 *      Author: dkoes
 */

#ifndef CAFFE_AFFINITY_LOSS_LAYER_HPP_
#define CAFFE_AFFINITY_LOSS_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "caffe/layers/loss_layer.hpp"

namespace caffe {

/**
 * @brief Custom layer. Inputs are predicted affinity and true affinity.
 * Euclidean loss is computed if there is a positive true affinity.
 * If the true affinity is negative, a hinge like L2 loss is computed
 * with respect to the absolute value (i.e., the predicted value is penalized
 * for being greater than the "true" value).
 * If the true affinity is zero, no penalty is assessed.
 *
 * This is *not* symmetric. The *second* input is assumed to be the true affinity.
 *
 */
template <typename Dtype>
class AffinityLossLayer : public LossLayer<Dtype> {
 public:
  explicit AffinityLossLayer(const LayerParameter& param)
      : LossLayer<Dtype>(param), diff_(), nranklosspairs(0) {}

  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline bool AllowForceBackward(const int bottom_index) const {
    return true;
  }

  virtual inline int ExactNumBottomBlobs() const {
    return 2 + this->layer_param_.affinity_loss_param().weighted();
  }

  virtual inline const char* type() const { return "AffinityLoss"; }

 protected:
  /// @copydoc AffinityLossLayer
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  /**
   * @brief Computes the error gradient w.r.t. the inputs.
   *
   */
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  virtual void Backward_relevance(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  Dtype compute_pair_loss(const vector<Blob<Dtype>*>& bottom, unsigned i, unsigned j, Dtype mult);
  void compute_pair_gradient(const vector<Blob<Dtype>*>& top,
          const vector<Blob<Dtype>*>& bottom, unsigned i, unsigned j, Dtype multi);

  Blob<Dtype> diff_;
  unsigned nranklosspairs;
};

}  // namespace caffe

#endif /* CAFFE_AFFINITY_LOSS_LAYER_HPP_ */
