#include <string>
#include <vector>

#include "caffe/sgd_solvers.hpp"
#include "caffe/util/hdf5.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/upgrade_proto.hpp"

namespace caffe {
  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::InputOptSGDPreSolve() {
    // do this better, also check std::map docs this is probably wrong as written
    if (this->net_->blob_names_index_.find("data") == this->net_->blob_names_index_.end())
      LOG(FATAL) << "Net doesn't have a data blob";
    int data_idx = this->net_->blob_names_index_["data"];
    input_blob = this->net_->blobs_()[data_idx];
    SGDSolver<Dtype>::PreSolve();
  }

  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::ComputeUpdateValue(Dtype rate) {
    const vector<float>& net_params_lr = this->net_->params_lr();
    Dtype momentum = this->param_.momentum();
    switch (Caffe::mode()) {
    case Caffe::CPU: {
      caffe_cpu_axpby(input_blob->count(), local_rate,
                input_blob->cpu_diff(), momentum,
                history_[param_id]->mutable_cpu_data());
      caffe_copy(input_blob->count(),
          history_[param_id]->cpu_data(),
          input_blob->mutable_cpu_diff());
      break;
    }
    case Caffe::GPU: {
  #ifndef CPU_ONLY
      sgd_update_gpu(input_blob->count(),
          input_blob->mutable_gpu_diff(),
          history_[param_id]->mutable_gpu_data(),
          momentum, local_rate);
  #else
      NO_GPU;
  #endif
      break;
    }
    default:
      LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
    }
  }

  // Update is very simple, no reason (I think) to normalize/regularize updates
  template <typename Dtype>
  void InputOptSGDSolver<Dtype>::ApplyUpdate() {
    Dtype rate = GetLearningRate();
    if (this->param_.display() && this->iter_ % this->param_.display() == 0) {
      LOG_IF(INFO, Caffe::root_solver()) << "Iteration " << this->iter_
          << ", lr = " << rate;
    }
    // do we want this? I guess it doesn't hurt to have it enabled, even if
    // param.clip_gradients() is always 0
    ClipGradients();
    ComputeUpdateValue(rate);
    this->net_->Update();
  }

  template <typename Dtype>
  void InputOptSGDSolver<Dtype>::ClipGradients() {
    const Dtype clip_gradients = this->param_.clip_gradients();
    if (clip_gradients < 0) { return; }
    Dtype sumsq_diff = input_blob->sumsq_diff();
    const Dtype l2norm_diff = std::sqrt(sumsq_diff);
    if (l2norm_diff > clip_gradients) {
      Dtype scale_factor = clip_gradients / l2norm_diff;
      LOG(INFO) << "Gradient clipping: scaling down gradients (L2 norm "
          << l2norm_diff << " > " << clip_gradients << ") "
          << "by scale factor " << scale_factor;
      input_blob->scale_diff(scale_factor);
    }
  }

  // When snapshotting, include net data blob...can probably refactor base
  // class to simplify or even remove this
  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::SnapshotSolverStateToBinaryProto(const string& model_filename) {
    SolverState state;
    state.set_iter(this->iter_);
    state.set_learned_net(model_filename);
    state.set_current_step(this->current_step_);
    state.clear_history();
    for (int i = 0; i < history_.size(); ++i) {
      // Add history
      BlobProto* history_blob = state.add_history();
      history_[i]->ToProto(history_blob);
    }
    BlobProto* datablob = state.add_datablob();
    input_blob->ToProto(datablob);
    string snapshot_filename = Solver<Dtype>::SnapshotFilename(".solverstate");
    LOG(INFO)
      << "Snapshotting solver state to binary proto file " << snapshot_filename;
    WriteProtoToBinaryFile(state, snapshot_filename.c_str());
  }

  template <typename Dtype> 
  void InputOptSGDSolver<Dtype>::SnapshotSolverStateToHDF5(const string& model_filename) {
    string snapshot_filename =
        Solver<Dtype>::SnapshotFilename(".solverstate.h5");
    LOG(INFO) << "Snapshotting solver state to HDF5 file " << snapshot_filename;
    hid_t file_hid = H5Fcreate(snapshot_filename.c_str(), H5F_ACC_TRUNC,
        H5P_DEFAULT, H5P_DEFAULT);
    CHECK_GE(file_hid, 0)
        << "Couldn't open " << snapshot_filename << " to save solver state.";
    hdf5_save_int(file_hid, "iter", this->iter_);
    hdf5_save_string(file_hid, "learned_net", model_filename);
    hdf5_save_int(file_hid, "current_step", this->current_step_);
    hid_t history_hid = H5Gcreate2(file_hid, "history", H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT);
    CHECK_GE(history_hid, 0)
        << "Error saving solver state to " << snapshot_filename << ".";
    for (int i = 0; i < history_.size(); ++i) {
      ostringstream oss;
      oss << i;
      hdf5_save_nd_dataset<Dtype>(history_hid, oss.str(), *history_[i]);
    }
    hid_t input_hid = H5Gcreate2(file_hid, "inputblob", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_GD(input_hid, 0)
        << "Error saving solver state to " << snapshot_filename << ".";
    ostringstream oss;
    hdf5_save_nd_dataset<Dtype>(input_hid, oss.str(), *input_blob);

    H5Gclose(history_hid);
    H5Fclose(file_hid);
  }
} // namespace caffe
