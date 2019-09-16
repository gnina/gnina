#include "test_solver.h"
#include <vector>
#include "common.h"
#include "caffe/blob.hpp"
#include "caffe/net.hpp"
#include "caffe/layers/molgrid_data_layer.hpp"

typedef caffe::BaseMolGridDataLayer<float, GridMaker> mgridT;

void test_iopt_update(caffe::Solver<float>* solver) {
  caffe::InputOptSolver<float>* iopt_solver = dynamic_cast<caffe::InputOptSolver<float>* >(solver);
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
  caffe::InputOptSolver<float>* iopt_solver = dynamic_cast<caffe::InputOptSolver<float>* >(solver);
  if (!iopt_solver) {
    throw usage_error("InputOptSolver test passed non-InputOptSolver ");
  }
  caffe::shared_ptr<caffe::Net<float> > net = iopt_solver->net();
  // do a forward pass to get initial loss
  float start_loss;
  net->Forward(&start_loss);
  std::cout << "start loss is " << start_loss << "\n";
  // do a full solve
  iopt_solver->Solve();
  // do final forward pass to get loss after several training updates
  float end_loss;
  net->Forward(&end_loss);
  std::cout << "end loss is " << end_loss << "\n";
  BOOST_CHECK_LT(end_loss, start_loss);
}

void test_iopt_exclude_rec(caffe::Solver<float>* solver) {
  caffe::InputOptSolver<float>* iopt_solver = dynamic_cast<caffe::InputOptSolver<float>* >(solver);
  if (!iopt_solver) {
    throw usage_error("InputOptSolver test passed non-InputOptSolver ");
  }
  // check that when exclude_receptor is set updates touch only the ligand channels
  caffe::Blob<float> init_grid;
  caffe::shared_ptr<caffe::Net<float> > net = iopt_solver->net();

  const auto& input_blob = iopt_solver->input_blob();
  init_grid.CopyFrom(*input_blob, false, true);
  const vector<caffe::shared_ptr<caffe::Layer<float> > >& layers = net->layers();
  mgridT* mgrid = dynamic_cast<caffe::BaseMolGridDataLayer<float, GridMaker>*>(layers[0].get());
  if (mgrid == NULL) {
    throw usage_error("First layer of model must be MolGridDataLayer.");
  }
  iopt_solver->SetNrecTypes(mgrid->getNumRecTypes());
  iopt_solver->SetNpoints(mgrid->getNumGridPoints());
  iopt_solver->Solve();
  
  unsigned nlt = mgrid->getNumLigTypes();
  unsigned nrt = mgrid->getNumRecTypes();
  unsigned npts = mgrid->getNumGridPoints();

  // rec grid points were not updated
  const float* current = input_blob->cpu_data();
  const float* init = init_grid.cpu_data();
  for (size_t i=0; i<nrt; ++i) {
    bool failed = false;
#pragma omp parallel for shared(failed)
    for (size_t j=0; j<npts; ++j) {
      if ((current[i*npts + j] - init[i*npts+j]) >= 0.001) {
        std::cout << current[i*npts + j] - init[i*npts+j] << "\n";
        failed = true;
      }
    }
    std::cout << "failure index " << i << "\n";
    BOOST_CHECK_EQUAL(failed, false);
  }
  // lig grid points were updated
  unsigned nt = nlt + nrt;
  bool failed = true;
  for (size_t i=nrt; i<nt; ++i) {
#pragma omp parallel for shared(failed)
    for (size_t j=0; j<npts; ++j) {
      if ((current[i*npts + j] - init[i*npts+j]) != 0) 
        failed = false;
    }
  }
  BOOST_CHECK_EQUAL(failed, false);
}
