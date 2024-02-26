/*
 * dl_scorer.h
 *
 *  Base clasee for neural network scoring functions.
 *
 *  Created on: Dec 21, 2023
 *      Author: dkoes
 */

#ifndef SRC_LIB_DL_SCORER_H_
#define SRC_LIB_DL_SCORER_H_

#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <memory>
#include <vector>

#include "model.h"

// this must be called in every thread that uses CNNScorer
extern int initializeCUDA(int device);

class DLScorer {

protected:
  std::shared_ptr<boost::recursive_mutex> mtx; //todo, enable parallel scoring

  // Receptor and ligand information
  std::vector<float3> ligand_coords, receptor_coords;
  std::vector<smt> ligand_smtypes, receptor_smtypes;
  std::vector<int> ligand_map, receptor_map; // from index in above arrays to index in movable atoms

  std::size_t num_atoms = 0; // Number of movable atoms + inflexible atoms
  vec current_center;        // center last time set_center was called, if min frame is moving, the mgrid center will be
                             // changing
  cnn_options cnnopts;

  // Set ligand and receptor atoms and coordinates from model
  virtual void setLigand(const model &m);
  virtual void setReceptor(const model &m);

public:
  DLScorer(): mtx(new boost::recursive_mutex), current_center(NAN,NAN,NAN) {}
  DLScorer(const cnn_options& opts) : mtx(new boost::recursive_mutex), current_center(NAN,NAN,NAN), cnnopts(opts)  {}

  virtual ~DLScorer() {}

  virtual bool initialized() const = 0;

  virtual bool has_affinity() const = 0; // return true if can predict affinity

  virtual const cnn_options &options() const { return cnnopts; }

  virtual float score(model &m, float &variance) = 0; // score only - no gradient
  virtual float score(model &m, bool compute_gradient, float &affinity, float &loss, float &variance) = 0;

  // readjust center
  virtual void set_center_from_model(model &m);

  virtual vec get_center() const { return current_center; }

  // disable receptor movement (e.g. for score only)
  void freeze_receptor() {
    if (isnan(cnnopts.cnn_center[0]))
      cnnopts.move_minimize_frame = true; // only move if center not specified
    cnnopts.fix_receptor = true;
  }

  // if model has bounding box, set it (leave bbox unchanged otherwise)
  virtual void set_bounding_box(grid_dims &box) const = 0;

  virtual std::shared_ptr<DLScorer> fresh_copy() const = 0;
};

#endif /* SRC_LIB_DL_SCORER_H_ */
