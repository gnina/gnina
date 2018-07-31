#include <assert.h>
#include "test_utils.h"
#include "test_cnn.h"
#include "atom_constants.h"
#include "gridmaker.h"
#include "cnn_scorer.h"
#include <cuda_runtime.h>
#include "quaternion.h"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/rng.hpp"
#include <boost/multi_array/multi_array_ref.hpp>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

extern parsed_args p_args;

void test_set_atom_gradients() {
  // randomly generate gridpoint gradients, accumulate for atoms, and then compare
  p_args.log << "CNN Set Atom Gradients Test \n";
  p_args.log << "Using random seed: " << p_args.seed << '\n';
  p_args.log << "Iteration " << p_args.iter_count << '\n';
  std::mt19937 engine(p_args.seed);
  std::normal_distribution<> diff_dist(0, 5);
  auto gen = std::bind(diff_dist, engine);

  std::vector<atom_params> mol_atoms;
  std::vector<smt> mol_types;
  make_mol(mol_atoms, mol_types, engine, 0, 1, 5000, 11.5, 11.5, 11.5);
  cnn_options cnnopts;
  cnnopts.cnn_scoring = true;
  model m;
  CNNScorer cnn_scorer(cnnopts, m);
  typedef CNNScorer::Dtype Dtype;
  unsigned ntypes = smt::NumTypes;
  caffe::GenericMolGridDataLayer<Dtype>* mgrid = 
    dynamic_cast<caffe::GenericMolGridDataLayer<Dtype>*>(cnn_scorer.mgrid);
  assert(mgrid);
  mgrid->batch_transform.resize(1);
  mgrid->batch_transform[0] =
      caffe::BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform();
  caffe::BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform& transform =
      mgrid->batch_transform[0];
  vec center(0, 0, 0);
  for (size_t i = 0; i < mol_atoms.size(); ++i) {
    atom_params& ainfo = mol_atoms[i];
    transform.mol.atoms.push_back(
        make_float4(ainfo.coords.x, ainfo.coords.y, ainfo.coords.z,
            xs_radius(mol_types[i])));
    transform.mol.whichGrid.push_back(mol_types[i]);
    transform.mol.gradient.push_back(make_float3(0, 0, 0));
    center += vec(ainfo.coords.x, ainfo.coords.y, ainfo.coords.z);
  }
  center /= mol_atoms.size();
  transform.center = center;
  transform.Q = quaternion(1, 0, 0, 0);
  caffe::Caffe::set_random_seed(p_args.seed);
  transform.set_random_quaternion(caffe::caffe_rng());
  //initialize gmaker
  GridMaker gmaker;
  bool spherize = false;
  double dim = round(mgrid->dimension / mgrid->resolution) + 1;
  gmaker.initialize(mgrid->resolution, mgrid->dimension, mgrid->radiusmultiple,
      mgrid->binary, spherize);
  gmaker.setCenter(center[0], center[1], center[2]);

  //randomly intialize gridpoint gradients
  std::vector<float> diff((ntypes) * dim * dim * dim);
  generate(begin(diff), end(diff), gen);

  //set up and calculate CPU atom gradients
  caffe::BaseMolGridDataLayer<Dtype, GridMaker>::Grids grids(&diff[0], 
          boost::extents[ntypes][dim][dim][dim]);
  caffe::BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform cpu_transform =
      mgrid->batch_transform[0];
  gmaker.setAtomGradientsCPU(cpu_transform.mol.atoms, cpu_transform.mol.whichGrid,
          cpu_transform.Q.boost(), grids, cpu_transform.mol.gradient);

  //calculate GPU atom gradients
  float* gpu_diff;
  cudaMalloc(&gpu_diff, diff.size() * sizeof(float));
  cudaMemcpy(gpu_diff, &diff[0], diff.size() * sizeof(float),
      cudaMemcpyHostToDevice);
  mgrid->setAtomGradientsGPU(gmaker, gpu_diff, 1);

  //compare results
  for (size_t i = 0; i < mol_atoms.size(); ++i) {
    for (size_t j = 0; j < 3; ++j) {
      p_args.log << "CPU " << cpu_transform.mol.gradient[i][j] << " GPU "
          << transform.mol.gradient[i][j] << "\n";
      BOOST_REQUIRE_SMALL(
          cpu_transform.mol.gradient[i][j] - transform.mol.gradient[i][j],
          (float )0.01);
    }
  }
  cudaFree(gpu_diff);
}

void test_subcube_grids() {
  //randomly generate mol, randomly choose a subgrid size, generate full grid
  //and subcube grid, check that the grid points match
  p_args.log << "CNN Subcube Grids Test \n";
  p_args.log << "Using random seed: " << p_args.seed << '\n';
  p_args.log << "Iteration " << p_args.iter_count << '\n';
  std::mt19937 engine(p_args.seed);

  //first find subcube size and full grid size (via multiplicative factor
  //to ensure that full grid and subgrids have the same number of points) 
  //idk i'm just picking a max arbitrarily here
  std::uniform_int_distribution<unsigned> subcubedim_dist(1, 15); 
  std::uniform_int_distribution<unsigned> gridfactor_dist(1, 10);
  std::uniform_int_distribution<unsigned> res_dist(1, 20);
  unsigned factor = gridfactor_dist(engine);
  float resolution = res_dist(engine) * 0.05;
  float subcube_dim = subcubedim_dist(engine) * resolution;
  float dimension = (subcube_dim + resolution) * factor + subcube_dim;
  unsigned dim = std::round(dimension / resolution) + 1; //number of grid points on a side
  typedef CNNScorer::Dtype Dtype;
  unsigned ntypes = smt::NumTypes;

  //allow batch size to be > 1
  std::uniform_int_distribution<unsigned> batch_dist(1, 17);
  unsigned batch_size = batch_dist(engine);

  //set params that don't change first
  cnn_options cnnopts;
  cnnopts.cnn_scoring = true;
  model m;
  CNNScorer cnn_scorer(cnnopts, m);
  caffe::GenericMolGridDataLayer<Dtype>* mgrid = 
    dynamic_cast<caffe::GenericMolGridDataLayer<Dtype>*>(cnn_scorer.mgrid);
  assert(mgrid);
  caffe::MolGridDataParameter* param = mgrid->layer_param_.mutable_molgrid_data_param();
  param->set_dimension(dimension);
  param->set_resolution(resolution);
  param->set_batch_size(batch_size);
  param->set_subgrid_dim(subcube_dim);
  mgrid->batch_transform.resize(batch_size);

  unsigned example_size = ntypes * dim * dim * dim;
  unsigned gsize = batch_size * example_size;

  //full grid
  Dtype* data = new Dtype[gsize];
  GridMaker gmaker;
  gmaker.initialize(*param);
  //cpu subgrids
  Dtype* rnndata = new Dtype[gsize];
  RNNGridMaker rnngmaker;
  rnngmaker.initialize(*param);
  rnngmaker.ntypes = ntypes;
  //gpu subgrids (don't share gridmaker because it's batch-aware)
  Dtype* rnndata_gpu = new Dtype[gsize];
  Dtype* gpu_grids;
  CUDA_CHECK(cudaMalloc(&gpu_grids, sizeof(Dtype) * gsize));
  RNNGridMaker rnngmaker_gpu;
  rnngmaker_gpu.initialize(*param);
  rnngmaker_gpu.ntypes = ntypes;

  unsigned& grids_per_dim = rnngmaker.grids_per_dim;
  unsigned& subgrid_dim_in_points = rnngmaker.subgrid_dim_in_points;
  unsigned ngrids = grids_per_dim * grids_per_dim * grids_per_dim;

  for (unsigned batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    //make mol
    std::vector<atom_params> mol_atoms;
    std::vector<smt> mol_types;
    make_mol(mol_atoms, mol_types, engine, 0, 1, 5000, dimension/2, dimension/2, dimension/2);

    mgrid->batch_transform[batch_idx] =
        caffe::BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform();
    caffe::BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform& transform =
        mgrid->batch_transform[batch_idx];
    vec center(0, 0, 0);
    for (size_t i = 0; i < mol_atoms.size(); ++i) {
      atom_params& ainfo = mol_atoms[i];
      transform.mol.atoms.push_back(
          make_float4(ainfo.coords.x, ainfo.coords.y, ainfo.coords.z,
              xs_radius(mol_types[i])));
      transform.mol.whichGrid.push_back(mol_types[i]);
      transform.mol.gradient.push_back(make_float3(0, 0, 0));
      center += vec(ainfo.coords.x, ainfo.coords.y, ainfo.coords.z);
    }
    center /= mol_atoms.size();
    transform.center = center;
    transform.Q = quaternion(1, 0, 0, 0);
    caffe::Caffe::set_random_seed(p_args.seed);
    transform.set_random_quaternion(caffe::caffe_rng());

    //full grid
    int offset = batch_idx * example_size;
    gmaker.setCenter(center[0], center[1], center[2]);
    gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q.boost(), 
        data+offset, ntypes);

    //now subcube grid, CPU first
    rnngmaker.setCenter(center[0], center[1], center[2]);
    rnngmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q.boost(), 
        &rnndata[0], ntypes);

    //set subcube grid atoms, GPU version
    unsigned natoms = transform.mol.atoms.size();
    mgrid->allocateGPUMem(natoms);
    CUDA_CHECK(cudaMemcpy(mgrid->gpu_gridatoms, &transform.mol.atoms[0], 
          natoms*sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(mgrid->gpu_gridwhich, &transform.mol.whichGrid[0], 
          natoms*sizeof(short), cudaMemcpyHostToDevice));
    rnngmaker_gpu.setCenter(center[0], center[1], center[2]);
    rnngmaker_gpu.setAtomsGPU(natoms, mgrid->gpu_gridatoms, 
        mgrid->gpu_gridwhich, transform.Q.boost(), ntypes, gpu_grids);
  }

  boost::multi_array_ref<Dtype, 6> rnngrids(rnndata, 
      boost::extents[ngrids][batch_size][ntypes][subgrid_dim_in_points]
      [subgrid_dim_in_points][subgrid_dim_in_points]);

  CUDA_CHECK(cudaMemcpy(rnndata_gpu, gpu_grids, sizeof(Dtype) * gsize, 
        cudaMemcpyDeviceToHost));
  boost::multi_array_ref<Dtype, 6> rnngrids_gpu(rnndata_gpu, 
      boost::extents[ngrids][batch_size][ntypes][subgrid_dim_in_points]
      [subgrid_dim_in_points][subgrid_dim_in_points]);

  //compare them
  for (size_t batch_idx = 0; batch_idx < batch_size; ++batch_idx) {
    for (size_t type = 0; type < ntypes; ++type) {
      for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
          for (size_t k = 0; k < dim; ++k) {
            int offset = batch_idx * example_size;
            boost::multi_array_ref<Dtype, 4> grids(data+offset, 
                boost::extents[ntypes][dim][dim][dim]);
            unsigned subgrid_idx_x = i / subgrid_dim_in_points; 
            unsigned subgrid_idx_y = j / subgrid_dim_in_points; 
            unsigned subgrid_idx_z = k / subgrid_dim_in_points; 
            unsigned rel_x = i % subgrid_dim_in_points; 
            unsigned rel_y = j % subgrid_dim_in_points; 
            unsigned rel_z = k % subgrid_dim_in_points; 
            unsigned grid_idx = (((subgrid_idx_x * grids_per_dim) + 
                  subgrid_idx_y) * grids_per_dim + subgrid_idx_z);
            p_args.log << "batch_idx: " << batch_idx << " grid_idx: " << grid_idx << 
              " type: " << type << " rel_x: " << rel_x << " rel_y: " << rel_y << 
              " rel_z: " << rel_z << " CPU full grid " << grids[type][i][j][k] << 
              " CPU subcube grid " << rnngrids[grid_idx][batch_idx][type][rel_x][rel_y][rel_z] << 
              " GPU subcube grid " << rnngrids_gpu[grid_idx][batch_idx][type][rel_x][rel_y][rel_z] << "\n";
            BOOST_REQUIRE_SMALL(
                grids[type][i][j][k] - rnngrids[grid_idx][batch_idx][type][rel_x][rel_y][rel_z],
                (float )0.01);
            BOOST_REQUIRE_SMALL(
                grids[type][i][j][k] - rnngrids_gpu[grid_idx][batch_idx][type][rel_x][rel_y][rel_z],
                (float )0.01);
          }
        }
      }
    }
  }

  delete data;
  delete rnndata;
  delete rnndata_gpu;
}
