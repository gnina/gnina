#ifndef CAFFE_NDIM_DATA_LAYER_HPP_
#define CAFFE_NDIM_DATA_LAYER_HPP_

#include <string>
#include <utility>
#include <vector>

#include "caffe/blob.hpp"
#include "caffe/data_transformer.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/base_data_layer.hpp"
#include "caffe/proto/caffe.pb.h"

namespace caffe {

/**
 * @brief Provides data to the Net from n-dimension  files of raw floating point data.
 *
 * TODO(dox): thorough documentation for Forward and proto params.
 */
template <typename Dtype>
class NDimDataLayer : public BasePrefetchingDataLayer<Dtype> {
 public:
  explicit NDimDataLayer(const LayerParameter& param)
      : BasePrefetchingDataLayer<Dtype>(param), actives_pos_(0),
        decoys_pos_(0), all_pos_(0), example_size(0), num_rotations(0), current_rotation(0) {}
  virtual ~NDimDataLayer();
  virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "NDimData"; }
  virtual inline int ExactNumBottomBlobs() const { return 0; }
  virtual inline int ExactNumTopBlobs() const { return 2; }

 protected:
  shared_ptr<Caffe::RNG> prefetch_rng_;
  virtual void Shuffle();
  virtual void load_batch(Batch<Dtype>* batch);

  virtual void load_data_from_files(Dtype*, const std::string& root, const vector<std::string>& files);
  const vector<int> blob2vec(const BlobShape& b) const;

  virtual void rotate_data(Dtype *data, unsigned rot);
  vector<vector<std::string> > actives_;
  vector<vector<std::string> > decoys_;
  vector<std::pair< vector<std::string>, int> > all_;
  int actives_pos_, decoys_pos_, all_pos_;
  unsigned example_size;
  unsigned num_rotations;
  unsigned current_rotation;
  vector<int> top_shape;

};


}  // namespace caffe

#endif  // CAFFE_NDIM_DATA_LAYER_HPP_
