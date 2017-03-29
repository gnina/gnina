#ifndef CAFFE_RANK_LOSS_LAYER_HPP_
#define CAFFE_RANK_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "caffe/layers/loss_layer.hpp"

namespace caffe {

/**
 * @brief Implements the RankNet pairwise loss described by Burges et al
 * (http://dl.acm.org/citation.cfm?id=1102363).
 *
 * This essentially computes the logistic loss of the _difference_ between
 * pairs of examples.  Pairs are either immediately adjacent or the n^2 pairs are considered.
 *
 *

 * @param bottom input Blob vector (length 2)
 *   -# @f$ (N \times (1 or2) ) @f$
 *      the predictions @f$ \hat{p} @f$, a Blob with values in
 *      @f$ [0, 1] @f$ indicating the predicted probability of the
 *      2 classes (for binary classification) or a single score.
 *   -# @f$ (N \times 1 \times 1 \times 1) @f$
 *      the labels @f$ l @f$, for binary classification an integer-valued Blob with values
 *      @f$ l_n \in [0, 1] @f$
 *      indicating the correct class label.  For single score, floating point values.
 *      In either case, these are used to determine the correct ordering between pairs.
 * @param top output Blob vector (length 1)
 *   -# @f$ (1 \times 1 \times 1 \times 1) @f$
 *      the computed  loss.
 *      @f$
 */
template <typename Dtype>
class RankLossLayer : public LossLayer<Dtype> {
 public:
  explicit RankLossLayer(const LayerParameter& param)
      : LossLayer<Dtype>(param), allpairs(false) {
    allpairs = this->layer_param_.rank_loss_param().allpairs();
  }
  virtual ~RankLossLayer() {}

  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "RankLossLayer"; }

 protected:
  /// @copydoc BinomialNoisyABNLossLayer
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  /**
   * @brief Computes the rank loss error gradient w.r.t. the
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
  bool allpairs;  // evaluate all possible pairs in batch rather than only adjacent

  Dtype compute_pair_loss(const vector<Blob<Dtype>*>& bottom, unsigned i, unsigned j);
  void compute_pair_gradient(const vector<Blob<Dtype>*>& top,
		  const vector<Blob<Dtype>*>& bottom, unsigned i, unsigned j);
};

}  // namespace caffe

#endif  // CAFFE_RANK_LOSS_LAYER_HPP_
