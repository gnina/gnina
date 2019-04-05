#ifndef CAFFE_DISCRETIZE_LAYER_HPP_
#define CAFFE_DISCRETIZE_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

namespace caffe {

/**
 * @brief A layer that discretizes the numerical input into categorical labels.
 * Bins are created by specifying the number of bins and min and max values.
 * Over or under flowing values are put in the highest or lowest bin.
 * The output is an integer vector suitable for use with the softmaxwithloss layer.
 *
 */
template <typename Dtype>
class DiscretizeLayer : public Layer<Dtype> {
 public:
  explicit DiscretizeLayer(const LayerParameter& param)
      : Layer<Dtype>(param), N_(0), min_value(0), max_value(0) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "DiscretizeLayer"; }
  virtual inline int ExactNumBottomBlobs() const { return 1; }
  virtual inline int ExactNumTopBlobs() const { return 1; }

 protected:
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  //don't actually do anything, although I suppose a gradient calculation is possible
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {}

  int N_;
  double min_value;
  double max_value;
};

}  // namespace caffe

#endif  // CAFFE_DISCRETIZE_LAYER_HPP_
