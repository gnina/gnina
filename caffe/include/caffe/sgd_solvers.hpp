#ifndef CAFFE_SGD_SOLVERS_HPP_
#define CAFFE_SGD_SOLVERS_HPP_

#include <string>
#include <vector>

#include "caffe/solver.hpp"
#include "caffe/layers/pooling_layer.hpp"

namespace caffe {

/**
 * @brief Optimizes the parameters of a Net using
 *        stochastic gradient descent (SGD) with momentum.
 */
template <typename Dtype>
class SGDSolver : public Solver<Dtype> {
 public:
  explicit SGDSolver(const SolverParameter& param)
      : Solver<Dtype>(param) { PreSolve(); }
  explicit SGDSolver(const string& param_file)
      : Solver<Dtype>(param_file) { PreSolve(); }
  explicit SGDSolver() {}
  virtual inline const char* type() const { return "SGD"; }

  const vector<shared_ptr<Blob<Dtype> > >& history() { return history_; }

 protected:
  void PreSolve();
  Dtype GetLearningRate();
  virtual void ApplyUpdate();
  virtual void Normalize(int param_id);
  virtual void Regularize(int param_id);
  virtual void ComputeUpdateValue(int param_id, Dtype rate);
  virtual void ClipGradients();
  virtual void SnapshotSolverState(const string& model_filename);
  virtual void SnapshotSolverStateToBinaryProto(const string& model_filename);
  virtual void SnapshotSolverStateToHDF5(const string& model_filename);
  virtual void RestoreSolverStateFromHDF5(const string& state_file);
  virtual void RestoreSolverStateFromBinaryProto(const string& state_file);
  // history maintains the historical momentum data.
  // update maintains update related data and is not needed in snapshots.
  // temp maintains other information that might be needed in computation
  //   of gradients/updates and is not needed in snapshots
  vector<shared_ptr<Blob<Dtype> > > history_, update_, temp_;

  DISABLE_COPY_AND_ASSIGN(SGDSolver);
};

template <typename Dtype>
class NesterovSolver : public SGDSolver<Dtype> {
 public:
  explicit NesterovSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param) {}
  explicit NesterovSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file) {}
  virtual inline const char* type() const { return "Nesterov"; }

 protected:
  virtual void ComputeUpdateValue(int param_id, Dtype rate);

  DISABLE_COPY_AND_ASSIGN(NesterovSolver);
};

template <typename Dtype>
class AdaGradSolver : public SGDSolver<Dtype> {
 public:
  explicit AdaGradSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param) { constructor_sanity_check(); }
  explicit AdaGradSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file) { constructor_sanity_check(); }
  virtual inline const char* type() const { return "AdaGrad"; }

 protected:
  virtual void ComputeUpdateValue(int param_id, Dtype rate);
  void constructor_sanity_check() {
    CHECK_EQ(0, this->param_.momentum())
        << "Momentum cannot be used with AdaGrad.";
  }

  DISABLE_COPY_AND_ASSIGN(AdaGradSolver);
};


template <typename Dtype>
class RMSPropSolver : public SGDSolver<Dtype> {
 public:
  explicit RMSPropSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param) { constructor_sanity_check(); }
  explicit RMSPropSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file) { constructor_sanity_check(); }
  virtual inline const char* type() const { return "RMSProp"; }

 protected:
  virtual void ComputeUpdateValue(int param_id, Dtype rate);
  void constructor_sanity_check() {
    CHECK_EQ(0, this->param_.momentum())
        << "Momentum cannot be used with RMSProp.";
    CHECK_GE(this->param_.rms_decay(), 0)
        << "rms_decay should lie between 0 and 1.";
    CHECK_LT(this->param_.rms_decay(), 1)
        << "rms_decay should lie between 0 and 1.";
  }

  DISABLE_COPY_AND_ASSIGN(RMSPropSolver);
};

template <typename Dtype>
class AdaDeltaSolver : public SGDSolver<Dtype> {
 public:
  explicit AdaDeltaSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param) { AdaDeltaPreSolve(); }
  explicit AdaDeltaSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file) { AdaDeltaPreSolve(); }
  virtual inline const char* type() const { return "AdaDelta"; }

 protected:
  void AdaDeltaPreSolve();
  virtual void ComputeUpdateValue(int param_id, Dtype rate);

  DISABLE_COPY_AND_ASSIGN(AdaDeltaSolver);
};

/**
 * @brief AdamSolver, an algorithm for first-order gradient-based optimization
 *        of stochastic objective functions, based on adaptive estimates of
 *        lower-order moments. Described in [1].
 *
 * [1] D. P. Kingma and J. L. Ba, "ADAM: A Method for Stochastic Optimization."
 *     arXiv preprint arXiv:1412.6980v8 (2014).
 */
template <typename Dtype>
class AdamSolver : public SGDSolver<Dtype> {
 public:
  explicit AdamSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param) { AdamPreSolve();}
  explicit AdamSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file) { AdamPreSolve(); }
  virtual inline const char* type() const { return "Adam"; }

 protected:
  void AdamPreSolve();
  virtual void ComputeUpdateValue(int param_id, Dtype rate);

  DISABLE_COPY_AND_ASSIGN(AdamSolver);
};

/**
 * @brief Solver that optimizes input grids using weights from a pre-trained
 * network, i.e. "dreaming." Compared with the base SGD solver, the InputOpt
 * Solver stores and loads the input data blob when checkpointing, and it only
 * updates the input data blob when updating. Can also be used as part of
 * adversarially robust training.
 *
 */
template <typename Dtype>
class InputOptSolver : public SGDSolver<Dtype> {
 public:
  explicit InputOptSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param), data_idx_(-1), input_blob_(nullptr), l_inf_(false), 
        threshold_update_(true),  exclude_ligand_(false), exclude_receptor_(false), 
        threshold_value_(Dtype(0)), nrec_types(0), nlig_types(0), 
        npoints_(0) { InputOptPreSolve(); } 
  explicit InputOptSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file), data_idx_(-1), input_blob_(nullptr), l_inf_(false), 
        threshold_update_(true), exclude_ligand_(false), exclude_receptor_(false), 
        threshold_value_(Dtype(0)),
        nrec_types(0), nlig_types(0), npoints_(0) { InputOptPreSolve(); }
  explicit InputOptSolver(const SolverParameter& param, shared_ptr<Net<Dtype> > net_, 
        vector<shared_ptr<Net<Dtype> > > test_nets_, Dtype eps_) : 
        data_idx_(-1), input_blob_(nullptr), l_inf_(true), threshold_update_(true), 
        exclude_ligand_(false), exclude_receptor_(true), threshold_value_(eps_), 
        nrec_types(0), nlig_types(0), npoints_(0) {
       InputOptPreSolve();
  }
  virtual inline const char* type() const { return "InputOpt"; }
  void ResetIter() { this->iter_ = 0; }
  shared_ptr<Blob<Dtype> > input_blob() const { return input_blob_; }
  void DoThreshold(bool threshold_update) { this->threshold_update_ = threshold_update; }
  void SetThresholdValue(Dtype threshold_value) { threshold_value_ = threshold_value; }
  void SetNrecTypes(unsigned ntypes) { nrec_types = ntypes; }
  void SetNligTypes(unsigned ntypes) { nlig_types = ntypes; }
  void SetNpoints(unsigned npoints) { npoints_ = npoints; }
  void SetLInf(bool l_inf) { l_inf_ = l_inf; }

 protected:
  void InputOptPreSolve();
  virtual void Step(int iters);
  virtual void ComputeUpdateValue(Dtype rate);
  virtual void ComputeUpdateValue(int param_id, Dtype rate) { NOT_IMPLEMENTED; }
  virtual void ApplyUpdate();
  virtual void ClipGradients();
  virtual void SnapshotSolverStateToBinaryProto(const string& model_filename);
  virtual void SnapshotSolverStateToHDF5(const string& model_filename);
  PoolingLayer<Dtype>* ToggleMaxToAve();
  void ThresholdBlob(shared_ptr<Blob<Dtype> >& tblob);
  void DoThresholdGPU(Dtype* offset_tblob, size_t blobsize, Dtype threshold_value);
  int data_idx_;
  shared_ptr<Blob<Dtype> > input_blob_;
  bool l_inf_; // take max steps within l_inf ball
  bool threshold_update_;
  bool exclude_ligand_;
  bool exclude_receptor_;
  Dtype threshold_value_; 
  unsigned nrec_types;
  unsigned nlig_types;
  unsigned npoints_;

  DISABLE_COPY_AND_ASSIGN(InputOptSolver);
};

/**
 * @brief Attempts to learn a more robust model that is resistant to
 * adversarial attacks using Projected Gradient Descent, training on 
 * the worst-case attacks that can be generated using the provided 
 * training examples. Described in [1]. 
 *
 * [1] Madry et al, "Towards Deep Learning Models Resistant to Adversarial
 * Attacks." arXiv preprint arXiv:1706.06083 (2019). 
 *
 * It does this by containing an internal inputopt solver for properly perturbing 
 * each training example, then using those "attacks" for training.
 *
 */
template <typename Dtype>
class MinMaxSolver : public SGDSolver<Dtype> {
 public:
  explicit MinMaxSolver(const SolverParameter& param)
      : SGDSolver<Dtype>(param), a_(0.1), eps_(0.3), k_(40) { 
        MinMaxPreSolve(); 
        SolverParameter tmparam = this->param();
        tmparam.set_lr_policy("fixed");
        tmparam.set_base_lr(a_);
        tmparam.set_momentum(Dtype(0));
        adversary_ = InputOptSolver<Dtype>(tmparam, this->net(), 
        this->test_nets(), eps_);
      }
  explicit MinMaxSolver(const string& param_file)
      : SGDSolver<Dtype>(param_file), adversary_(param_file), 
      a_(0.1), eps_(0.3), k_(40) { 
        MinMaxPreSolve(); 
        SolverParameter tmparam = this->param();
        tmparam.set_lr_policy("fixed");
        tmparam.set_base_lr(a_);
        tmparam.set_momentum(Dtype(0));
        adversary_ = InputOptSolver<Dtype>(tmparam, this->net(), 
        this->test_nets(), eps_);
      }
  virtual inline const char* type() const { return "MinMax"; }

 protected:
  void MinMaxPreSolve();
  virtual void Step(int iters); // looks like generic SGDSolver::Step() except it calls out to the InputOptSolver before starting, and does ForwardBackward starting and ending after the data layer 
  int data_idx_ = -1;
  InputOptSolver<Dtype> adversary_;
  float a_; // maximization problem step size
  float eps_; // maximization problem distance from initial point
  int k_; // maximization problem number of steps

  DISABLE_COPY_AND_ASSIGN(MinMaxSolver);
};

}  // namespace caffe

#endif  // CAFFE_SGD_SOLVERS_HPP_
