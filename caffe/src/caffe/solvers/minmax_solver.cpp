#include "caffe/sgd_solvers.hpp"

namespace caffe {
template <typename Dtype>
void MinMaxSolver<Dtype>::MinMaxPreSolve() {
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
}

template <typename Dtype>
void MinMaxSolver<Dtype>::Step(int iters) {
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
    // go forward from data layer to set input blob 
    this->net_->ForwardFromTo(data_idx_, data_idx_);
    // do PGD via the InputOpt Solver adversary
    adversary_.ResetIter();
    adversary_.Step(this->param_.k());
    
    // accumulate the loss and gradient
    Dtype loss = 0;
    for (int i = 0; i < this->param_.iter_size(); ++i) {
      // go forward from layer after data layer (don't touch data, which was
      // perturbed by adversary)
      loss += this->net_->ForwardFrom(data_idx_+1); 
      this->net_->Backward();
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
    this->ApplyUpdate();

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

INSTANTIATE_CLASS(MinMaxSolver);
REGISTER_SOLVER_CLASS(MinMax);

} // namespace caffe
