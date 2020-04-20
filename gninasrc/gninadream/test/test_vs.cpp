#include "test_vs.h"
#include "caffe/blob.hpp"
#include "caffe/net.hpp"
#include "caffe/layers/molgrid_data_layer.hpp"
#include "../lib/cnn_scorer.h"

using namespace caffe;

void test_do_exact_vs() {
  // set up some simple grids and make sure the output is correct
  // both including and excluding the rec channels
  // 16 rec types, 19 lig types
  // use three different examples, check result when optimized grids are
  // - set equal to example
  //   - L2 should be 0
  //   - mult should be sum of grid^2, normalized
  //   - threshold should be the number of nonzero grid points in example, normalized
  // - set equal to 0
  //   - L2 should be the rmsd of the reference density
  //   - mult should be 0
  //   - threshold should be 0
  // - set equal to 1
  //   - L2 should be reference mean deviation from 1
  //   - mult should be the sum of all the density values, normalized
  //   - thresh should be be the number of nonzero grid points in example, normalized
  // - set equal to -1
  //   - L2 should be reference mean deviation from -1
  //   - mult should be sum of all density values, negated and normalized
  //   - thresh should be number of nonzero grid points in example, negated and normalized
  //   expected order by score is
  //   1-L2: 3 < 2 < 1 < 0 (i think)
  //   mult: 3 < 1 < 0 < 2 
  //   thresh: 3 < 1 < 0,2
  // comparing across the different examples is harder - even if we control for
  // the number of atoms in a ligand, the different atom types have different
  // radii, which means atoms of different types have different amounts of density, 
  // effectively rewarding/penalizing certain types compared with others. we
  // probably _should_ normalize this effect away
  
  // set up params etc
  LayerParameter param;
  std::string vsfile;
  std::vector<std::string> ligs = {"aa2ar_ligand.mol2", "cdk2_ligand.mol2", "esr1_ligand.mol2"};
  std::vector<std::string> recs = {"aa2ar_receptor.pdbqt", "cdk2_receptor.pdbqt", "esr1_receptor.pdbqt"};
  std::vector<caffe::shared_ptr<std::ostream> > out;
  bool gpu=1;
  float positive_threshold=0.25;
  float negative_threshold=0.25;
  std::vector<std::string> dist_method = {"l2", "mult", "threshold"};
  out.push_back(boost::make_shared<std::ofstream>("tmp.vsout"));
  caffe::Caffe::SetDevice(0);
  caffe::Caffe::set_mode(caffe::Caffe::GPU);
  caffe::Caffe::set_cudnn(false);

  // set up net
  NetParameter net_param;
  cnn_options cnnopts;
  cnnopts.cnn_scoring = true;
  const char *cmodel = cnn_models[cnnopts.cnn_model_name].model;
  google::protobuf::io::ArrayInputStream modeldata(cmodel, strlen(cmodel));
  bool success = google::protobuf::TextFormat::Parse(&modeldata, &net_param);
  if (!success) throw usage_error("Error with built-in cnn model "+cnnopts.cnn_model_name);
  UpgradeNetAsNeeded("default", &net_param);
  net_param.mutable_state()->set_phase(TRAIN);

  LayerParameter* first = net_param.mutable_layer(1);
  MolGridDataParameter* mgridparam = first->mutable_molgrid_data_param();
  if (mgridparam == NULL) {
    throw usage_error("First layer of model must be MolGridData.");
  }

  mgridparam->set_random_rotation(false);
  mgridparam->set_random_translate(false);
  mgridparam->set_shuffle(false);
  mgridparam->set_balanced(false);
  mgridparam->set_stratify_receptor(false);
  mgridparam->set_inmemory(true);

  const char *recmap = cnn_models[cnnopts.cnn_model_name].recmap;
  const char *ligmap = cnn_models[cnnopts.cnn_model_name].ligmap;
  mgridparam->set_mem_recmap(recmap);
  mgridparam->set_mem_ligmap(ligmap);
  caffe::Net<float> net(net_param);

  NetParameter wparam;

  const unsigned char *weights = cnn_models[cnnopts.cnn_model_name].weights;
  unsigned int nweights = cnn_models[cnnopts.cnn_model_name].num_bytes;

  google::protobuf::io::ArrayInputStream weightdata(weights,nweights);
  google::protobuf::io::CodedInputStream strm(&weightdata);
  strm.SetTotalBytesLimit(INT_MAX, 536870912);
  success = wparam.ParseFromCodedStream(&strm);
  if (!success) throw usage_error("Error with default weights.");

  net.CopyTrainedLayersFrom(wparam);

  const vector<caffe::shared_ptr<Layer<float> > >& layers = net.layers();
  mgridT* mgrid = dynamic_cast<MolGridDataLayer<float>*>(layers[0].get());
  if (mgrid == NULL) {
    throw usage_error("First layer of model must be MolGridDataLayer.");
  }
  unsigned ligGridSize = mgrid->getNumGridPoints() * mgrid->getNumLigandTypes();
  unsigned recGridSize = mgrid->getNumGridPoints()* mgrid->getNumReceptorTypes();
  std::vector<float> fillbuffer(ligGridSize, 0);
  std::vector<float> refgrid(ligGridSize, 0);

  // find data blob
  const auto& blob_names = net.blob_names();
  int data_idx = -1;
  for (size_t i=0; i<blob_names.size(); ++i) {
    if (blob_names[i] == "data") {
      data_idx = i;
      break;
    }
  }
  if (data_idx < 0)
    LOG(FATAL) << "Net doesn't have a data blob";
  caffe::shared_ptr<Blob<float> > input_blob = net.blobs()[data_idx];
  float* dest = input_blob->mutable_gpu_data();

  // set up molgetter/model
  tee log(true);
  FlexInfo finfo(log);
  std::vector<float> expected_scores;
  for (size_t i=0; i<recs.size(); ++i) {
    MolGetter mols(recs[i], std::string(), finfo, true, true, log);
    mols.setInputFile(ligs[i]);
    model m;
    mols.readMoleculeIntoModel(m);
    // set reference density equal to ligand density
    std::vector<smt> ligand_smtypes;
    std::vector<float3> ligand_coords;
    setLigand(m, ligand_coords, ligand_smtypes);
    mgrid->setLigand(ligand_coords, ligand_smtypes);
    std::vector<smt> rec_smtypes;
    std::vector<float3> rec_coords;
    setReceptor(m, rec_coords, rec_smtypes);
    mgrid->setReceptor(rec_coords, rec_smtypes);
    mgrid->setLabels(1, 10); 
    unsigned natoms = m.num_movable_atoms();
    net.ForwardFromTo(0,0);
    CUDA_CHECK_GNINA(cudaMemcpy(&refgrid[0], dest + recGridSize, ligGridSize * sizeof(float), cudaMemcpyDeviceToHost));
    std::vector<std::string> this_lig = {ligs[i]};
    for (const auto& method : dist_method) {
      do_exact_vs(*net_param.mutable_layer(1), net, this_lig[0], this_lig, out, gpu, 
          method, positive_threshold, negative_threshold);
      if (!std::strcmp(method.c_str(), "l2"))
        expected_scores.push_back(0);
      if (!std::strcmp(method.c_str(), "mult")) {
        float sum;
        cpu_mult(&refgrid[0], &refgrid[0], &sum, ligGridSize);
        expected_scores.push_back(sum / (natoms * ligGridSize));
      }
      if (!std::strcmp(method.c_str(), "threshold")) {
        float sum;
        cpu_thresh(&refgrid[0], &refgrid[0], &sum, ligGridSize, positive_threshold, negative_threshold);
        expected_scores.push_back(sum / (natoms * ligGridSize));
      }
    }

    // set reference equal to 0
    std::fill(&fillbuffer.front(), &fillbuffer.back(), 0);
    CUDA_CHECK(cudaMemcpy(dest + recGridSize, &fillbuffer[0], ligGridSize* sizeof(float), cudaMemcpyHostToDevice));
    for (const auto& method : dist_method) {
      do_exact_vs(*net_param.mutable_layer(1), net, this_lig[0], this_lig, out, gpu, 
          method, positive_threshold, negative_threshold);
      if (!std::strcmp(method.c_str(), "l2")) {
        float sum;
        cpu_l2sq(&refgrid[0], &refgrid[0], &sum, ligGridSize);
        expected_scores.push_back(std::sqrt(sum) / (natoms * ligGridSize));
      }
      if (!std::strcmp(method.c_str(), "mult"))
        expected_scores.push_back(0);
      if (!std::strcmp(method.c_str(), "threshold")) 
        expected_scores.push_back(0);
    }

    // set reference equal to 1
    std::fill(&fillbuffer.front(), &fillbuffer.back(), 1);
    CUDA_CHECK(cudaMemcpy(dest + recGridSize, &fillbuffer[0], ligGridSize* sizeof(float), cudaMemcpyHostToDevice));
    for (const auto& method : dist_method) {
      do_exact_vs(*net_param.mutable_layer(1), net, this_lig[0], this_lig, out, gpu, 
          method, positive_threshold, negative_threshold);
      if (!std::strcmp(method.c_str(), "l2")) {
        float sum;
        cpu_l2sq(&refgrid[0], &refgrid[0], &sum, ligGridSize);
        expected_scores.push_back(std::sqrt(sum) / (natoms * ligGridSize));
      }
      if (!std::strcmp(method.c_str(), "mult")) {
        float sum = 0.;
#pragma omp parallel for reduction(+:sum)
        for (size_t k=0; k<ligGridSize; ++k)
          sum = sum + refgrid[k];
        expected_scores.push_back(sum / (natoms * ligGridSize));
      }
      if (!std::strcmp(method.c_str(), "threshold")) {
        float sum = 0.;
        float val;
#pragma omp parallel for reduction(+:sum)
        for (size_t k=0; k<ligGridSize; ++k) {
          val = refgrid[k] > 0 ? 1 : 0;
          sum = sum + val;
        }
        expected_scores.push_back(sum / (natoms * ligGridSize));
      }
    }

    // set reference equal to -1
    std::fill(&fillbuffer.front(), &fillbuffer.back(), -1);
    CUDA_CHECK(cudaMemcpy(dest + recGridSize, &fillbuffer[0], ligGridSize* sizeof(float), cudaMemcpyHostToDevice));
    for (const auto& method : dist_method) {
      do_exact_vs(*net_param.mutable_layer(1), net, this_lig[0], this_lig, out, gpu, 
          method, positive_threshold, negative_threshold);
      if (!std::strcmp(method.c_str(), "l2")) {
        float sum;
        cpu_l2sq(&refgrid[0], &refgrid[0], &sum, ligGridSize);
        expected_scores.push_back(std::sqrt(sum) / (natoms * ligGridSize));
      }
      if (!std::strcmp(method.c_str(), "mult")) {
        float sum = 0.;
#pragma omp parallel for reduction(+:sum)
        for (size_t k=0; k<ligGridSize; ++k)
          sum = sum + refgrid[k];
        expected_scores.push_back(-sum / (natoms * ligGridSize));
      }
      if (!std::strcmp(method.c_str(), "threshold")) {
        float sum = 0.;
        float val;
#pragma omp parallel for reduction(+:sum)
        for (size_t k=0; k<ligGridSize; ++k) {
          val = refgrid[k] > 0 ? 1 : 0;
          sum = sum + val;
        }
        expected_scores.push_back(-sum / (natoms * ligGridSize));
      }
    }
  }

  out.clear();
  // read the scores back in...annoying
  std::ifstream infile("tmp.vsout");
  CHECK((bool)infile) << "Could not open vsout file";
  string line;
  float score;
  std::vector<float> scores;
  while (getline(infile, line)) {
    stringstream stream(line);
    stream >> score;
    scores.push_back(score);
  }
  BOOST_CHECK_EQUAL(expected_scores.size(), scores.size());

  for (size_t i=0; i<expected_scores.size(); ++i) {
    BOOST_REQUIRE_SMALL(scores[i] - expected_scores[i], (float)0.00001);
  }
}
