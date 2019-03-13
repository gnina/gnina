#include <string>
#include <vector>

#include "caffe/sgd_solvers.hpp"
#include "caffe/util/hdf5.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/upgrade_proto.hpp"

namespace caffe {
  // Update is very simple, no reason (I think) to normalize/regularize updates
  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::ComputeUpdateValue(int param_id, Dtype rate) {}

  // When snapshotting, include net data blob
  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::SnapshotSolverState(const string& model_filename) {}

  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::SnapshotSolverStateToBinaryProto(const string& model_filename) {}

  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::SnapshotSolverStateToHDF5(const string& model_filename) {}

  // If we have a data blob from a checkpoint, update the net
  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::RestoreSolverStateFromHDF5() {}

  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::RestoreSolverStateFromBinaryProto(const string& state_file) {
    SGDSolver<Dtype>::RestoreSolverStateFromBinaryProto(state_file);
  }
} // namespace caffe
