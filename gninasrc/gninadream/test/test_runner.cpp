#include <boost/program_options.hpp>
#include <iostream>
#include <cuda_runtime.h>
#include "parsed_args.h"
#include "test_loss.h"
#include <iostream>
#define N_ITERS 5
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

namespace ut = boost::unit_test;
namespace po = boost::program_options;

parsed_args p_args;

void boost_loop_test(void (*func)());

BOOST_AUTO_TEST_SUITE(loss)

BOOST_AUTO_TEST_CASE(cpu_loss) {
  test_cpu_l2();
}

BOOST_AUTO_TEST_CASE(gpu_loss) {
  test_gpu_l2();
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
  int gpu;
  po::positional_options_description positional;
  po::options_description inputs("Input");
  inputs.add_options()("seed", po::value<unsigned>(&p_args.seed),
      "seed for random number generator")("n_iters",
      po::value<unsigned>(&p_args.n_iters),
      "number of iterations to repeat relevant tests")("log",
      po::value<std::string>(&logname), "specify logfile, default is test.log")("gpu", 
      po::value<int>(&gpu));
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
  std::string _dumvar3;
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
        po::value<int>(&_dumvar4));
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
