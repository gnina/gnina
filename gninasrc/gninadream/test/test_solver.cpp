#include "test_solver.h"
#include <vector>
#include "common.h"
#include "caffe/blob.hpp"
#include "caffe/net.hpp"

void test_iopt_update(caffe::Solver<float>* solver) {
  caffe::InputOptSGDSolver<float>* iopt_solver = dynamic_cast<caffe::InputOptSGDSolver<float>* >(solver);
  if (!iopt_solver) {
    throw usage_error("InputOptSolver test passed non-InputOptSolver ");
  }
  // check that updates touch the input blob and only the input blob
  std::vector<std::shared_ptr<caffe::Blob<float> > > init_params;
  caffe::Blob<float> init_grid;
  caffe::shared_ptr<caffe::Net<float> > net = iopt_solver->net();
  const auto& params = net->learnable_params();
  size_t nparams = params.size();
  for (size_t i=0; i<nparams; ++i) {
    init_params.push_back(std::make_shared<caffe::Blob<float> >());
    init_params[i]->CopyFrom(*params[i], false, true);
  }
  const auto& input_blob = iopt_solver->input_blob();
  init_grid.CopyFrom(*input_blob, false, true);
  iopt_solver->Solve();
  // no params got updated
  for (size_t i=0; i<nparams; ++i) {
    const float* current = params[i]->cpu_data();
    const float* init = init_params[i]->cpu_data();
    int count = params[i]->count();
    bool failed = false;
#pragma omp parallel for shared(failed)
    for (size_t j=0; j<count; ++j) {
      if ((*(current+j) - *(init+j)) != 0) 
        failed = true;
    }
    BOOST_CHECK_EQUAL(failed, false);
  }
  // top blob did get updated
  int count = input_blob->count();
  const float* current = input_blob->cpu_data();
  const float* init = init_grid.cpu_data();
  bool failed = true;
#pragma omp parallel for shared(failed)
  for (size_t j=0; j<count; ++j) {
    if ((*(current+j) - *(init+j)) != 0)
      failed = false;
  }
  BOOST_CHECK_EQUAL(failed, false);
}

void test_iopt_improvement(caffe::Solver<float>* solver) {
  caffe::InputOptSGDSolver<float>* iopt_solver = dynamic_cast<caffe::InputOptSGDSolver<float>* >(solver);
  if (!iopt_solver) {
    throw usage_error("InputOptSolver test passed non-InputOptSolver ");
  }
  // check that after several iterations the loss has improved
  iopt_solver->Solve();
}
