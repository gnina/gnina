/*
 * cnn_torch_scorer.h
 *
 *  Created on: Dec 20, 2023
 *      Author: dkoes
 */

#ifndef SRC_LIB_CNN_TORCH_SCORER_H_
#define SRC_LIB_CNN_TORCH_SCORER_H_

#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <memory>
#include <torch/script.h>
#include <vector>

#include "tee.h"
#include "dl_scorer.h"
#include "model.h"
#include "torch_model.h"

/* This class evaluates protein-ligand poses according to a provided
 * torch convolutional neural net (CNN) model.
 */
template<bool isCUDA>
class CNNTorchScorer : public DLScorer {
public:
  typedef float Dtype;

private:
  // scratch vectors to avoid memory reallocation
  std::vector<gfloat3> gradient, gradient_lig, gradient_rec;
  std::vector<std::shared_ptr<TorchModel<isCUDA> > > models;

  void getGradient(std::shared_ptr<TorchModel<isCUDA> >, std::vector<gfloat3>& gradient);
public:
  CNNTorchScorer() {}
  virtual ~CNNTorchScorer() {}

  CNNTorchScorer(const cnn_options &opts, tee *log = nullptr);

  bool initialized() const { return models.size() > 0; }

  bool has_affinity() const; // return true if can predict affinity

  float score(model &m, float &variance); // score only - no gradient
  float score(model &m, bool compute_gradient, float &affinity, float &loss, float &variance);

  fl get_grid_dim() const;
  fl get_grid_res() const;

  void set_bounding_box(grid_dims &box) const;

  std::shared_ptr<DLScorer> fresh_copy() const { return std::make_shared<CNNTorchScorer>(cnnopts); }

protected:
};

#endif /* SRC_LIB_CNN_SCORER_H_ */
