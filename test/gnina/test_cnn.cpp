#include <assert.h>
#include <cmath>
#include "test_utils.h"
#include "test_cnn.h"
#include "atom_constants.h"
#include "cnn_scorer.h"
#include <cuda_runtime.h>
#include "quaternion.h"
#include "caffe/proto/caffe.pb.h"
#include "caffe/layers/flex_lstm_layer.hpp"
#include "caffe/util/rng.hpp"
#include "caffe/util/device_alternate.hpp"
#include <boost/multi_array/multi_array_ref.hpp>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#define TOL 0.0001f

extern parsed_args p_args;

using namespace caffe;
using namespace std;

void test_set_atom_gradients() {
  Caffe::set_mode(Caffe::GPU);
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
  cnnopts.cnn_scoring = CNNall;
  model m;
  CNNScorer cnn_scorer(cnnopts);
  typedef CNNScorer::Dtype Dtype;
  MolGridDataLayer<Dtype>* mgrid = cnn_scorer.get_mgrid();
  assert(mgrid);
  int ntypes = mgrid->getNumChannels();
  int dim = mgrid->getGridDims().x;
  mgrid->enableLigandGradients();

  set_cnn_grids(mgrid, mol_atoms, mol_types);

  vector<gfloat3> cpugrad, gpugrad;
  //randomly intialize gridpoint gradients
  Blob<Dtype> diff({1,ntypes,dim,dim,dim});
  Dtype* diffcpu = diff.mutable_cpu_diff();
  generate(diffcpu, diffcpu+diff.count(), gen);

  vector<Blob<Dtype>*> top{&diff};
  vector<Blob<Dtype>*> bottom; //not used

  //set up and calculate CPU atom gradients
  mgrid->backward(top, bottom, false);
  //store
  mgrid->getLigandGradient(0, cpugrad);

  //calculate GPU atom gradients
  mgrid->backward(top, bottom, true);
  mgrid->getLigandGradient(0, gpugrad);

  p_args.log << "NUMATOMS " << mol_atoms.size() << "\n";
  //compare results
  double sum = 0;
  for (size_t i = 0; i < mol_atoms.size(); ++i) {
    for (size_t j = 0; j < 3; ++j) {
      p_args.log << i <<"," << j << "  CPU " << cpugrad[i][j] << " GPU " << gpugrad[i][j] << "  " << mol_types[i] << " " << mol_atoms[i].coords << "\n";
      BOOST_REQUIRE_SMALL(cpugrad[i][j] - gpugrad[i][j], (float )0.01);
      sum += cpugrad[i].x+cpugrad[i].y+cpugrad[i].z;
    }
  }

  BOOST_REQUIRE_NE(sum, 0); //make sure something was actually calculated
}

void test_vanilla_grids() {
  //randomly generate mol, check that CPU and GPU grids match
  p_args.log << "CNN Vanilla Grids Test \n";
  p_args.log << "Using random seed: " << p_args.seed << '\n';
  p_args.log << "Iteration " << p_args.iter_count << '\n';
  std::mt19937 engine(p_args.seed);
  Caffe::set_mode(Caffe::GPU);

  //honestly just trying to replicate libmolgrid test to check whether it fails
  //here too, so being less stringent than I was for the other tests in terms
  //of randomization
  std::vector<atom_params> mol_atoms;
  std::vector<smt> mol_types;
  make_mol(mol_atoms, mol_types, engine);

  cnn_options cnnopts;
  cnnopts.cnn_scoring = CNNall;
  model m;
  CNNScorer cnn_scorer(cnnopts);
  typedef CNNScorer::Dtype Dtype;

  MolGridDataLayer<Dtype>* mgrid = cnn_scorer.get_mgrid();
  assert(mgrid);
  int ntypes = mgrid->getNumChannels();
  int dim = mgrid->getGridDims().x;

  set_cnn_grids(mgrid, mol_atoms, mol_types);

  vector<Blob<Dtype> > topblobs(mgrid->ExactNumTopBlobs());
  topblobs[0].Reshape({1,ntypes,dim,dim,dim});

  vector<Blob<Dtype>*> bottom;
  vector<Blob<Dtype>*> top;

  top.push_back(&topblobs[0]);
  for(unsigned i = 1; i < mgrid->ExactNumTopBlobs(); i++) {
    topblobs[i].Reshape({1,1});
    top.push_back(&topblobs[i]);
  }

  //get CPU grid
  mgrid->forward(bottom, top, false);
  //store cpu result
  vector<Dtype> cpuout(topblobs[0].count(), 0);
  copy(topblobs[0].cpu_data(),topblobs[0].cpu_data()+topblobs[0].count(), cpuout.begin());
  topblobs[0].Clear();
  //get gpu grid
  mgrid->forward(bottom, top, true);

  vector<Dtype> gpuout(topblobs[0].count(), 0);
  copy(topblobs[0].cpu_data(),topblobs[0].cpu_data()+topblobs[0].count(), gpuout.begin());

  BOOST_REQUIRE_EQUAL(cpuout.size(),gpuout.size());
  BOOST_REQUIRE_EQUAL(gpuout.size(), dim*dim*dim*ntypes);
  for(unsigned i = 0, n = gpuout.size(); i < n; i++) {
    BOOST_REQUIRE_SMALL(cpuout[i]-gpuout[i], TOL);
    p_args.log << "CPU " << cpuout[i] << " GPU " << gpuout[i] << "\n";
  }
}

//TODO TODO TODO: reimplement this functionality
#if 0
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
  CNNScorer cnn_scorer(cnnopts);
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
  SubcubeGridMaker rnngmaker;
  rnngmaker.initialize(*param);
  rnngmaker.ntypes = ntypes;
  //gpu subgrids (don't share gridmaker because it's batch-aware)
  Dtype* rnndata_gpu = new Dtype[gsize];
  Dtype* gpu_grids;
  CUDA_CHECK(cudaMalloc(&gpu_grids, sizeof(Dtype) * gsize));
  SubcubeGridMaker rnngmaker_gpu;
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

  delete[] data;
  delete[] rnndata;
  delete[] rnndata_gpu;
}

void test_strided_cube_datagetter() {
  typedef float Dtype;
  //check whether we're getting the correct subcube of the parent grid
  p_args.log << "CNN Strided Cube Handler Test \n";
  p_args.log << "Using random seed: " << p_args.seed << '\n';
  p_args.log << "Iteration " << p_args.iter_count << '\n';
  std::mt19937 engine(p_args.seed);

  //how big is subcube
  std::uniform_int_distribution<unsigned> subcubedim_dist(3, 15); 
  std::uniform_int_distribution<unsigned> res_dist(1, 20);
  float resolution = res_dist(engine) * 0.25;
  unsigned subcube_dim_pts = subcubedim_dist(engine);

  //full grid at least as big as subcube
  std::uniform_int_distribution<unsigned> griddim_dist(subcube_dim_pts, 50);
  unsigned dim = griddim_dist(engine);
  unsigned ntypes = smt::NumTypes;

  //get subcube stride
  std::uniform_int_distribution<unsigned> stride_dist(3, subcube_dim_pts);
  unsigned cube_stride = stride_dist(engine);

  //get batch size and resolve sizes of things
  std::uniform_int_distribution<unsigned> batch_dist(1, 32);
  unsigned batch_size = batch_dist(engine);
  unsigned example_size = ntypes * dim * dim * dim;
  unsigned gsize = batch_size * example_size;
  unsigned cube_size = batch_size * ntypes * subcube_dim_pts * 
    subcube_dim_pts * subcube_dim_pts;
  unsigned slices_per_dim = ((dim - subcube_dim_pts) / cube_stride) + 1;
  unsigned n_timesteps = slices_per_dim * slices_per_dim * slices_per_dim;
  p_args.log << "Num Timesteps " << n_timesteps << '\n';

  //set up flex layer
  caffe::LayerParameter param;
  caffe::LSTMDataGetterLayer<float> datagetter(param);
  datagetter.pattern = caffe::AccessPatterns::strided_cube;
  datagetter.num_timesteps = n_timesteps;
  datagetter.batch_size = batch_size;
  datagetter.ntypes = ntypes;
  datagetter.dim = dim;
  datagetter.subgrid_dim = subcube_dim_pts;
  datagetter.cube_stride = cube_stride;
  datagetter.example_size = example_size;
  datagetter.current_timestep = 0;

  //blob shapes: subcube is BxCxSub_dimxSub_dimxSub_dim
  //             full grid is BxCxDimxDimxDim
  std::vector<int> subcube_shape = {(int)batch_size, (int)ntypes, 
    (int)subcube_dim_pts, (int)subcube_dim_pts, (int)subcube_dim_pts};
  std::vector<int> full_shape = {(int)batch_size, (int)ntypes, (int)dim, 
    (int)dim, (int)dim};

  //blobs
  caffe::Blob<float> b_blob;
  caffe::Blob<float> t_blob;
  std::vector<caffe::Blob<float>*> bottom_blobs = {&b_blob};
  std::vector<caffe::Blob<float>*> top_blobs = {&t_blob};
  bottom_blobs[0]->Reshape(full_shape);
  top_blobs[0]->Reshape(subcube_shape);

  //stuff to make boost multi_array usage a little less verbose
  typedef boost::multi_array<Dtype, 5> array_t;
  typedef boost::multi_array_types::index_range range_t;
  array_t::index_gen indices;

  //randomly populate full grid
  float* in_data = bottom_blobs[0]->mutable_cpu_data();
  std::uniform_real_distribution<Dtype> dist;
  generate(in_data, in_data + gsize, bind(dist, engine));
  boost::multi_array_ref<Dtype, 5> in_data_ref(in_data, boost::extents[batch_size][ntypes] 
      [dim][dim][dim]);

  //storage for the CPU forward pass at each timestep, since we'll overwrite the
  //top blob when we do a second pass for the GPU
  caffe::Blob<float> cpu_subcube_buffer;
  cpu_subcube_buffer.Reshape(subcube_shape);

  for (unsigned ts=0; ts<n_timesteps; ++ts) {
    datagetter.current_timestep = ts;
    //update cpu subcube, copy out the result and make a multi_array ref for indexing
    datagetter.Forward_cpu(bottom_blobs, top_blobs);
    cpu_subcube_buffer.CopyFrom(*top_blobs[0]);
    boost::multi_array_ref<Dtype, 5> cpu_subcube(cpu_subcube_buffer.mutable_cpu_data(), 
        boost::extents[batch_size][ntypes] 
        [subcube_dim_pts][subcube_dim_pts][subcube_dim_pts]);
    //update gpu subcube, make ref to result
    datagetter.Forward_gpu(bottom_blobs, top_blobs);
    boost::multi_array_ref<Dtype, 5> gpu_subcube(top_blobs[0]->mutable_cpu_data(), 
        boost::extents[batch_size][ntypes] 
        [subcube_dim_pts][subcube_dim_pts][subcube_dim_pts]);
    //use timestep and slices per dim to figure out x,y,z ranges for
    //overall array
    unsigned factor = (((dim - subcube_dim_pts) / cube_stride) + 1);
    unsigned i = (ts / (factor * factor) % slices_per_dim) * cube_stride;
    unsigned j = ((ts / factor) % slices_per_dim) * cube_stride;
    unsigned k = (ts % slices_per_dim) * cube_stride;
    bool failed = false;
#pragma omp parallel for shared(failed)
    for (unsigned batch_idx=0; batch_idx < batch_size; ++batch_idx) {
      for (unsigned type=0; type<ntypes; ++type) {
        //get view of cube using boost and use this as the ground truth
        array_t::array_view<3>::type cube_view = in_data_ref[ indices[batch_idx][type][range_t(i,i+subcube_dim_pts)][range_t(j,j+subcube_dim_pts)][range_t(k,k+subcube_dim_pts)] ];
        for (unsigned x=0; x<subcube_dim_pts; ++x) 
          for (unsigned y=0; y<subcube_dim_pts; ++y) 
            for (unsigned z=0; z<subcube_dim_pts; ++z) {
              if (std::abs(cube_view[x][y][z] - cpu_subcube[batch_idx][type][x][y][z] > 0.01))
                failed = true;
              if (std::abs(cube_view[x][y][z] - gpu_subcube[batch_idx][type][x][y][z] > 0.01))
                failed = true;
            }
      }
    }
    BOOST_CHECK_EQUAL(failed, false);
  }
}
#endif
