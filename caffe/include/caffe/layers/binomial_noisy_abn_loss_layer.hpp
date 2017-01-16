#ifndef CAFFE_BINOMIAL_NOISY_ABN_LOSS_LAYER_HPP_
#define CAFFE_BINOMIAL_NOISY_ABN_LOSS_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "caffe/layers/loss_layer.hpp"

namespace caffe {

/**
 * @brief Implements the asymmetric Bernoulli noise model (ABN) as described by
 * Mnih and Hinton (http://machinelearning.wustl.edu/mlpapers/papers/ICML2012Mnih_318).
 *
 * This computes a cross entropy between two classes under the assumption of
 * an asymmetric Bernoulli noise model specified with the parameters
 * theta0 (probability of a false example being observed as true) and
 * theta1 (probability of true example being observed as false).
 *
 * Assumes the input is a two class probability distribution.
 *
 * TODO: Create a layer integrated with Softmax.
 *

 * @param bottom input Blob vector (length 2)
 *   -# @f$ (N \times 2) @f$
 *      the predictions @f$ \hat{p} @f$, a Blob with values in
 *      @f$ [0, 1] @f$ indicating the predicted probability of the
 *      2 classes.  Each prediction vector @f$ \hat{p}_n @f$
 *      should sum to 1 as in a probability distribution: @f$
 *      \forall n \sum\limits_{k=1}^K \hat{p}_{nk} = 1 @f$.
 *   -# @f$ (N \times 1 \times 1 \times 1) @f$
 *      the labels @f$ l @f$, an integer-valued Blob with values
 *      @f$ l_n \in [0, 1] @f$
 *      indicating the correct class label.
 * @param top output Blob vector (length 1)
 *   -# @f$ (1 \times 1 \times 1 \times 1) @f$
 *      the computed noise-adjusted cross entropy loss.
 *      @f$
 */
template <typename Dtype>
class BinomialNoisyABNLossLayer : public LossLayer<Dtype> {
 public:
  explicit BinomialNoisyABNLossLayer(const LayerParameter& param)
      : LossLayer<Dtype>(param), theta0(0), theta1(0) {
    theta0 = this->layer_param_.binomial_noisy_abn_loss_param().theta0();
    theta1 = this->layer_param_.binomial_noisy_abn_loss_param().theta1();
  }
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "BinomialNoisyABNLoss"; }

 protected:
  /// @copydoc BinomialNoisyABNLossLayer
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  /**
   * @brief Computes the ABN cross entropy loss error gradient w.r.t. the
   *        predictions.
   *
   * Gradients cannot be computed with respect to the label inputs (bottom[1]),
   * so this method ignores bottom[1] and requires !propagate_down[1], crashing
   * if propagate_down[1] is set.
   *
   * @param top output Blob vector (length 1), providing the error gradient with
   *      respect to the outputs
   * @param propagate_down see Layer::Backward.
   *      propagate_down[1] must be false as we can't compute gradients with
   *      respect to the labels.
   * @param bottom input Blob vector (length 2)
   *   -# @f$ (N \times C \times H \times W) @f$
   *      the predictions @f$ \hat{p} @f$; Backward computes diff
   *      @f$ \frac{\partial E}{\partial \hat{p}} @f$
   *   -# @f$ (N \times 1 \times 1 \times 1) @f$
   *      the labels -- ignored as we can't compute their error gradients
   */
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

 private:
  double theta0;  // probability of a false becoming true
  double theta1; //probability of a true becoming false
};

}  // namespace caffe

#endif  // CAFFE_BINOMIAL_NOISY_ABN_LOSS_LAYER_HPP_
