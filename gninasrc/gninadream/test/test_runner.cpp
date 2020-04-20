#include <boost/program_options.hpp>
#include <iostream>
#include <cuda_runtime.h>
#include "test/gnina/parsed_args.h"
#include "cnn_scorer.h"
#include "molgetter.h"
#include "tee.h"
#include "test/gnina/test_utils.h"
#include "test_loss.h"
#include "test_solver.h"
#include "test_vs.h"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include <iostream>
#define N_ITERS 5
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

namespace ut = boost::unit_test;
namespace po = boost::program_options;
using namespace caffe;
typedef MolGridDataLayer<float> mgridT;

parsed_args p_args;

void boost_loop_test(void (*func)());

struct f_iopt_solver {
  f_iopt_solver() {
    // set up basic inputopt solver with gnina default net
    std::mt19937 engine(p_args.seed);
    NetParameter net_param;
    std::string cnn_model_name = "default2017";
    const char *model = cnn_models[cnn_model_name].model;
    google::protobuf::io::ArrayInputStream modeldata(model, strlen(model));
    bool success = google::protobuf::TextFormat::Parse(&modeldata, &net_param);
    if (!success) throw usage_error("Error with built-in cnn model "+cnn_model_name);
    UpgradeNetAsNeeded("default", &net_param);
    net_param.mutable_state()->set_phase(TRAIN);
    net_param.set_force_backward(true);

    LayerParameter* first = net_param.mutable_layer(1);
    MolGridDataParameter* mgridparam = first->mutable_molgrid_data_param();
    if (mgridparam == NULL) {
      throw usage_error("First layer of model must be MolGridData.");
    }
    const char *recmap = cnn_models[cnn_model_name].recmap;
    const char *ligmap = cnn_models[cnn_model_name].ligmap;
    mgridparam->set_mem_recmap(recmap);
    mgridparam->set_mem_ligmap(ligmap);
    mgridparam->set_inmemory(true);
    mgridparam->set_batch_size(1);
    mgridparam->set_use_rec_center(true);
    mgridparam->set_ignore_ligand(true);

    SolverParameter solver_param;
    solver_param.mutable_net_param()->CopyFrom(net_param);
    solver_param.set_base_lr(1);
    solver_param.set_max_iter(100);
    solver_param.set_lr_policy("inv");
    solver_param.set_snapshot_after_train(false);
    solver_param.set_snapshot_prefix("inputopt");
    solver_param.set_type("InputOpt");

    solver = SolverRegistry<float>::CreateSolver(solver_param);
    
    boost::shared_ptr<Net<float> > net = solver->net();
    NetParameter wparam;

    const unsigned char *weights = cnn_models[cnn_model_name].weights;
    unsigned int nweights = cnn_models[cnn_model_name].num_bytes;

    google::protobuf::io::ArrayInputStream weightdata(weights,nweights);
    google::protobuf::io::CodedInputStream strm(&weightdata);
    strm.SetTotalBytesLimit(INT_MAX, 536870912);
    success = wparam.ParseFromCodedStream(&strm);
    if (!success) throw usage_error("Error with default weights.");

    net->CopyTrainedLayersFrom(wparam);
    const vector<caffe::shared_ptr<Layer<float> > >& layers = net->layers();
    mgridT* mgrid = dynamic_cast<MolGridDataLayer<float>*>(layers[0].get());
    if (mgrid == NULL) {
      throw usage_error("First layer of model must be MolGridDataLayer.");
    }

    // populate an atomv for rec and lig, set inmem
    tee log(true);
    ::FlexInfo finfo(log);
    ::MolGetter mols(p_args.rec, std::string(), finfo, true, true, log);
    mols.setInputFile("");
    ::model m;
    mols.readMoleculeIntoModel(m);
    std::vector<smt> rec_smtypes;
    std::vector<float3> rec_coords;
    setReceptor(m, rec_coords, rec_smtypes);
    mgrid->setReceptor(rec_coords, rec_smtypes);
    // std::vector<smt> ligand_smtypes;
    // std::vector<float3> ligand_coords;
    // setLigand(m, ligand_coords, ligand_smtypes);
    // mgrid->setLigand(ligand_coords, ligand_smtypes);
    mgrid->setLabels(1);
  }
  ~f_iopt_solver() {}

  Solver<float>* solver;
  atomv rec;
  atomv lig;
};

BOOST_AUTO_TEST_SUITE(loss)

BOOST_AUTO_TEST_CASE(cpu_loss) {
  test_cpu_l2();
}

BOOST_AUTO_TEST_CASE(gpu_loss) {
  test_gpu_l2();
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(solver, f_iopt_solver)

BOOST_AUTO_TEST_CASE(solver_update) {
  test_iopt_update(this->solver);
}

BOOST_AUTO_TEST_CASE(solver_improvement) {
  test_iopt_improvement(this->solver);
}

BOOST_AUTO_TEST_CASE(solver_exclude_rec) {
  test_iopt_exclude_rec(this->solver);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(vs)

BOOST_AUTO_TEST_CASE(vs_rank) {
  test_do_exact_vs();
}

BOOST_AUTO_TEST_SUITE_END()

void initializeCUDA(int device) {
  cudaError_t error;
  cudaDeviceProp deviceProp;

  error = cudaSetDevice(device);
  if (error != cudaSuccess) {
    std::cerr << "cudaSetDevice returned error code " << error << "\n";
    exit(-1);
  }

  error = cudaGetDevice(&device);

  if (error != cudaSuccess) {
    std::cerr << "cudaGetDevice returned error code " << error << "\n";
    exit(-1);
  }

  error = cudaGetDeviceProperties(&deviceProp, device);

  if (deviceProp.computeMode == cudaComputeModeProhibited) {
    std::cerr
        << "Error: device is running in <Compute Mode Prohibited>, no threads can use ::cudaSetDevice().\n";
    exit(-1);
  }

  if (error != cudaSuccess) {
    std::cerr << "cudaGetDeviceProperties returned error code " << error
        << "\n";
    exit(-1);
  }

}

bool init_unit_test() {
  // TODO: multithread running tests
  std::string logname;
  unsigned seed;
  int gpu = 0;
  po::positional_options_description positional;
  po::options_description inputs("Input");
  inputs.add_options()("seed", po::value<unsigned>(&p_args.seed),
      "seed for random number generator")("n_iters",
      po::value<unsigned>(&p_args.n_iters),
      "number of iterations to repeat relevant tests")("log",
      po::value<std::string>(&logname), "specify logfile, default is test.log")("gpu", 
      po::value<int>(&gpu), "gpu to use")
      ("rec", po::value<std::string>(&p_args.rec), 
      "specify receptor filename for input optimization solver test, default is 10gs_pocket.pdb");
  po::options_description desc, desc_simple;
  desc.add(inputs);
  desc_simple.add(inputs);
  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(ut::framework::master_test_suite().argc,
            ut::framework::master_test_suite().argv).options(desc).style(
            po::command_line_style::default_style
                ^ po::command_line_style::allow_guessing).positional(positional).run(),
        vm);
    notify(vm);
  } catch (po::error& e) {
    std::cerr << "Command line parse error: " << e.what() << '\n'
        << "\nCorrect usage:\n" << desc_simple << '\n';
    return false;
  }
  if (!vm.count("seed")) {
    p_args.seed = std::random_device()();
    p_args.many_iters = true;
  }
  caffe::Caffe::set_cudnn(false);
  if (!vm.count("n_iters")) p_args.n_iters = N_ITERS;
  if (!vm.count("log")) logname = "test.log";
  if (vm.count("gpu")) initializeCUDA(gpu);

  p_args.log.init(logname);

  p_args.params = {p_args.seed};
  if (p_args.many_iters) for (size_t i = 0; i < p_args.n_iters; ++i)
    p_args.params.push_back(std::random_device()());

  return true;
}

int main(int argc, char* argv[]) {
  //The following exists so that passing --help prints both the UTF help and our
  //specific program options - I can't see any way of doing this without
  //parsing the args twice, because UTF chomps args it thinks it owns
  unsigned _dumvar1;
  unsigned _dumvar2;
  std::string _dumvar3, _dumvar5;
  int _dumvar4;
  bool help = false;
  po::positional_options_description positional;
  po::options_description inputs("Input");
  inputs.add_options()("seed", po::value<unsigned>(&_dumvar1),
      "seed for random number generator")("n_iters",
      po::value<unsigned>(&_dumvar2),
      "number of iterations to repeat relevant tests")("log",
      po::value<std::string>(&_dumvar3),
      "specify logfile, default is test.log")("gpu", 
        po::value<int>(&_dumvar4), "gpu to use")(
      "rec", po::value<std::string>(&_dumvar5), 
      "specify receptor filename for input optimization solver test, default is 10gs_pocket.pdb");
  po::options_description info("Information");
  info.add_options()("help", po::bool_switch(&help), "print usage information");
  po::options_description desc, desc_simple;
  desc.add(inputs).add(info);
  desc_simple.add(inputs).add(info);
  po::variables_map vm;
  try {
    po::store(
        po::command_line_parser(argc, argv).options(desc).style(
            po::command_line_style::default_style
                ^ po::command_line_style::allow_guessing).positional(positional).run(),
        vm);
    notify(vm);
  } catch (po::error& e) {
  }

  if (help) std::cout << desc_simple << '\n';

  return ut::unit_test_main(&init_unit_test, argc, argv);
}
