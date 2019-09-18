#include <string>
#include <vector>

#include "caffe/sgd_solvers.hpp"
#include "caffe/util/hdf5.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/upgrade_proto.hpp"

namespace caffe {
  template <typename Dtype>
  void InputOptSolver<Dtype>::LInfInit(const SolverParameter& param, 
      const shared_ptr<Net<Dtype> >& net, float eps) {
    this->param_ = param;
    CHECK_GE(this->param_.average_loss(), 1) << "average_loss should be non-negative.";
    this->iter_ = 0;
    this->current_step_ = 0;
    this->net_ = net;
    threshold_value_ = eps;
    // net should have already been set up
    CHECK_NOTNULL(net.get());
    InputOptPreSolve();
  }

  template <typename Dtype> 
  void InputOptSolver<Dtype>::InputOptPreSolve() {
    const auto& blob_names = this->net_->blob_names();
    int data_idx = -1;
    for (size_t i=0; i<blob_names.size(); ++i) {
      if (blob_names[i] == "data") {
        data_idx = i;
        break;
      }
    }
    if (data_idx < 0)
      LOG(FATAL) << "Net doesn't have a data blob";
    data_idx_ = data_idx;
    input_blob_ = this->net_->blobs()[data_idx];
    this->history_.clear();
    this->update_.clear();
    this->temp_.clear();
    const vector<int>& shape = input_blob_->shape();
    this->history_.push_back(shared_ptr<Blob<Dtype> >(new Blob<Dtype>(shape)));
    this->update_.push_back(shared_ptr<Blob<Dtype> >(new Blob<Dtype>(shape)));
    this->temp_.push_back(shared_ptr<Blob<Dtype> >(new Blob<Dtype>(shape)));
  }

  template <typename Dtype>
  PoolingLayer<Dtype>* InputOptSolver<Dtype>::ToggleMaxToAve() {
    const vector<shared_ptr<Layer<Dtype> > >& layers = this->net_->layers();
    PoolingLayer<Dtype> *pool = NULL;
    for (unsigned i = 1, nl = layers.size(); i < nl; i++) {
      pool = dynamic_cast<PoolingLayer<Dtype>*>(layers[i].get());
      if (pool)
        break; 
      else
        if (layers[i]->type() == string("Convolution"))
          break; 
        else
          if (layers[i]->type() == string("InnerProduct")) break;
    }
  
    if (pool) {
      if (pool->pool() == PoolingParameter_PoolMethod_MAX) {
        pool->set_pool(PoolingParameter_PoolMethod_AVE);
      } else {
        pool = NULL; //no need to reset to max
      }
    }
    return pool;
  }

  template <typename Dtype>
  void InputOptSolver<Dtype>::ThresholdBlob(shared_ptr<Blob<Dtype> >& tblob) {
    // if threshold_value_ is nonzero, clip to [-val, +val]
    size_t blobsize = tblob->count() - nrec_types * npoints_ - nlig_types * npoints_;
    switch (Caffe::mode()) {
    case Caffe::CPU: {
      Dtype* tptr = tblob->mutable_cpu_data();
      Dtype* offset_tptr = tptr + nrec_types * npoints_;
        if (this->threshold_value_ == 0) {
    #pragma omp parallel for
        for (size_t i=0; i<blobsize; ++i) {
          if (*(offset_tptr + i) < 0)
            *(offset_tptr + i) = 0;
        }
      }
        else {
          Dtype clip_min = threshold_value_ > 0 ? -threshold_value_ : threshold_value_;
          Dtype clip_max = threshold_value_ > 0 ? threshold_value_ : -threshold_value_;
    #pragma omp parallel for
          for (size_t i=0; i<blobsize; ++i) {
            *(offset_tptr+i) = std::max(clip_min, std::min(*(offset_tptr+i), clip_max));
          }
        }
      break;
    }
    case Caffe::GPU: {
  #ifndef CPU_ONLY
      Dtype* tptr = tblob->mutable_gpu_data();
      Dtype* offset_tptr = tptr + nrec_types * npoints_;
      DoThresholdGPU(offset_tptr, blobsize, this->threshold_value_);
  #else
      NO_GPU;
  #endif
      break;
    }
    default:
      LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
    }
  }

  template <typename Dtype> 
  void InputOptSolver<Dtype>::Step(int iters) {
    const int start_iter = this->iter_;
    const int stop_iter = this->iter_ + iters;
    int average_loss = this->param_.average_loss();
    this->losses_.clear();
    this->smoothed_loss_ = 0;
    this->iteration_timer_.Start();

    while (this->iter_ < stop_iter) {
      // zero-init the params
      this->net_->ClearParamDiffs();
      if (this->param_.test_interval() && this->iter_ % this->param_.test_interval() == 0
          && (this->iter_ > 0 || this->param_.test_initialization())) {
        if (Caffe::root_solver()) {
          this->TestAll();
        }
        if (this->requested_early_exit_) {
          // Break out of the while loop because stop was requested while testing.
          break;
        }
      }

      for (int i = 0; i < this->callbacks_.size(); ++i) {
        this->callbacks_[i]->on_start();
      }
      const bool display = this->param_.display() && this->iter_ % this->param_.display() == 0;
      this->net_->set_debug_info(display && this->param_.debug_info());
      // accumulate the loss and gradient
      Dtype loss = 0;
      for (int i = 0; i < this->param_.iter_size(); ++i) {
        loss += this->net_->ForwardFrom(data_idx_+1);
        PoolingLayer<Dtype>* pool = ToggleMaxToAve();
        this->net_->Backward();
        if (pool) 
          pool->set_pool(PoolingParameter_PoolMethod_MAX);
      }
      loss /= this->param_.iter_size();
      // average the loss across iterations for smoothed reporting
      this->UpdateSmoothedLoss(loss, start_iter, average_loss);
      if (display) {
        float lapse = this->iteration_timer_.Seconds();
        float per_s = (this->iter_ - this->iterations_last_) / (lapse ? lapse : 1);
        LOG_IF(INFO, Caffe::root_solver()) << "Iteration " << this->iter_
            << " (" << per_s << " iter/s, " << lapse << "s/"
            << this->param_.display() << " iters), loss = " << this->smoothed_loss_;
        this->iteration_timer_.Start();
        this->iterations_last_ = this->iter_;
        const vector<Blob<Dtype>*>& result = this->net_->output_blobs();
        int score_index = 0;
        for (int j = 0; j < result.size(); ++j) {
          const Dtype* result_vec = result[j]->cpu_data();
          const string& output_name =
              this->net_->blob_names()[this->net_->output_blob_indices()[j]];
          const Dtype loss_weight =
              this->net_->blob_loss_weights()[this->net_->output_blob_indices()[j]];
          for (int k = 0; k < result[j]->count(); ++k) {
            ostringstream loss_msg_stream;
            if (loss_weight) {
              loss_msg_stream << " (* " << loss_weight
                              << " = " << loss_weight * result_vec[k] << " loss)";
            }
            LOG_IF(INFO, Caffe::root_solver()) << "    Train net output #"
                << score_index++ << ": " << output_name << " = "
                << result_vec[k] << loss_msg_stream.str();
          }
        }
      }
      for (int i = 0; i < this->callbacks_.size(); ++i) {
        this->callbacks_[i]->on_gradients_ready();
      }
      ApplyUpdate();

      // Increment the internal iter_ counter -- its value should always indicate
      // the number of times the weights have been updated.
      ++this->iter_;

      SolverAction::Enum request = this->GetRequestedAction();

      // Save a snapshot if needed.
      if ((this->param_.snapshot()
           && this->iter_ % this->param_.snapshot() == 0
           && Caffe::root_solver()) ||
           (request == SolverAction::SNAPSHOT)) {
        this->Snapshot();
      }
      if (SolverAction::STOP == request) {
        this->requested_early_exit_ = true;
        // Break out of training loop.
        break;
      }
    }
  }

#ifndef CPU_ONLY
template <typename Dtype>
void sgd_update_gpu(int N, Dtype* g, Dtype* h, Dtype momentum,
    Dtype local_rate);
#endif

  template <typename Dtype> 
  void InputOptSolver<Dtype>::ComputeUpdateValue(Dtype rate) {
    Dtype momentum = this->param_.momentum();
    size_t blobsize = input_blob_->count() - nrec_types * npoints_ - 
      nlig_types * npoints_;
    size_t ptr_offset = nrec_types * npoints_;
    // if we're stepping in LInf, update is rate * sgn(g)
    // so replace g with sgn(g) (in-place); also momentum should be 0 (set in
    // SolverParameter by MinMaxSolver)
    switch (Caffe::mode()) {
    case Caffe::CPU: {
      Dtype* diff_ptr = input_blob_->mutable_cpu_diff();
      if (nrec_types)
        std::fill(diff_ptr, diff_ptr+ptr_offset, 0);
      else if (nlig_types)
        std::fill(diff_ptr+blobsize, diff_ptr+blobsize+nlig_types*npoints_, 0);
      if (l_inf_) 
        caffe_cpu_sign(blobsize, input_blob_->cpu_diff(), 
            input_blob_->mutable_cpu_diff());
      caffe_cpu_axpby(blobsize, rate, input_blob_->cpu_diff() + ptr_offset,
          momentum, this->history_[0]->mutable_cpu_data() + ptr_offset);
      caffe_copy(blobsize, this->history_[0]->cpu_data() + ptr_offset,
          diff_ptr + ptr_offset);
      break;
    }
    case Caffe::GPU: {
  #ifndef CPU_ONLY
      Dtype* diff_ptr = input_blob_->mutable_gpu_diff();
      if (nrec_types) 
          CUDA_CHECK(cudaMemset(diff_ptr, 0, nrec_types * npoints_ * sizeof(Dtype)));
      if (nlig_types)
          CUDA_CHECK(cudaMemset(diff_ptr+blobsize, 0, nlig_types * npoints_ * sizeof(Dtype)));
      if (l_inf_) 
        caffe_gpu_sign(blobsize, input_blob_->gpu_diff(), 
            input_blob_->mutable_gpu_diff());
      sgd_update_gpu(blobsize, input_blob_->mutable_gpu_diff() + ptr_offset,
          this->history_[0]->mutable_gpu_data() + ptr_offset, momentum, rate);
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
  void InputOptSolver<Dtype>::ApplyUpdate() {
    Dtype rate = this->GetLearningRate(); // if L_inf stepping, this is a_
    if (this->param_.display() && this->iter_ % this->param_.display() == 0) {
      LOG_IF(INFO, Caffe::root_solver()) << "Iteration " << this->iter_
          << ", inputopt_lr = " << rate;
    }
    // do we want this? I guess it doesn't hurt to have it enabled, even if
    // param.clip_gradients() is always 0
    ClipGradients();
    ComputeUpdateValue(rate);
    if (l_inf_) 
      input_blob_->AscentUpdate();
    else
      input_blob_->Update();
    if (threshold_update_)
      ThresholdBlob(input_blob_);
  }

  template <typename Dtype>
  void InputOptSolver<Dtype>::ClipGradients() {
    const Dtype clip_gradients = this->param_.clip_gradients();
    if (clip_gradients < 0) { return; }
    Dtype sumsq_diff = input_blob_->sumsq_diff();
    const Dtype l2norm_diff = std::sqrt(sumsq_diff);
    if (l2norm_diff > clip_gradients) {
      Dtype scale_factor = clip_gradients / l2norm_diff;
      LOG(INFO) << "Gradient clipping: scaling down gradients (L2 norm "
          << l2norm_diff << " > " << clip_gradients << ") "
          << "by scale factor " << scale_factor;
      input_blob_->scale_diff(scale_factor);
    }
  }

  // When snapshotting, include net data blob...can probably refactor base
  // class to simplify or even remove this
  template <typename Dtype> 
  void InputOptSolver<Dtype>::SnapshotSolverStateToBinaryProto(const string& model_filename) {
    SolverState state;
    state.set_iter(this->iter_);
    state.set_learned_net(model_filename);
    state.set_current_step(this->current_step_);
    state.clear_history();
    for (int i = 0; i < this->history_.size(); ++i) {
      // Add history
      BlobProto* history_blob = state.add_history();
      this->history_[i]->ToProto(history_blob);
    }
    input_blob_->ToProto(state.mutable_datablob());
    string snapshot_filename = Solver<Dtype>::SnapshotFilename(".solverstate");
    LOG(INFO)
      << "Snapshotting solver state to binary proto file " << snapshot_filename;
    WriteProtoToBinaryFile(state, snapshot_filename.c_str());
  }

  template <typename Dtype> 
  void InputOptSolver<Dtype>::SnapshotSolverStateToHDF5(const string& model_filename) {
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
    for (int i = 0; i < this->history_.size(); ++i) {
      ostringstream oss;
      oss << i;
      hdf5_save_nd_dataset<Dtype>(history_hid, oss.str(), *this->history_[i]);
    }
    hid_t input_hid = H5Gcreate2(file_hid, "inputblob", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_GE(input_hid, 0)
        << "Error saving solver state to " << snapshot_filename << ".";
    ostringstream oss;
    hdf5_save_nd_dataset<Dtype>(input_hid, oss.str(), *input_blob_);

    H5Gclose(history_hid);
    H5Fclose(file_hid);
  }

INSTANTIATE_CLASS(InputOptSolver);
REGISTER_SOLVER_CLASS(InputOpt);

} // namespace caffe
