#include <boost/program_options.hpp>
#include <iostream>
#include "parsed_args.h"
#include "device_buffer.h"
#include "test_gpucode.h"
#include "test_tree.h"
#include "test_cache.h"
#include "test_utils.h"
#define N_ITERS 5
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

namespace ut = boost::unit_test;
namespace po = boost::program_options;

parsed_args p_args;
bool run_on_gpu = true;
int cuda_dev_id = 0;

void boost_loop_test(void (*func)());

//TODO: when we are running with a boost version > 1.58, start using UTF
//datasets with BOOST_PARAM_TEST_CASE

BOOST_AUTO_TEST_SUITE(gpucode)

BOOST_AUTO_TEST_CASE(interaction_energy) {
  boost_loop_test(&test_interaction_energy);
}

BOOST_AUTO_TEST_CASE(eval_intra) {
  boost_loop_test(&test_eval_intra);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_tree_gpu)

BOOST_AUTO_TEST_CASE(set_conf) {
  boost_loop_test(&test_set_conf);
}

BOOST_AUTO_TEST_CASE(derivative) {
  boost_loop_test(&test_derivative);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(cache_gpu)

BOOST_AUTO_TEST_CASE(eval_deriv) {
  boost_loop_test(&test_cache_eval_deriv);
}

BOOST_AUTO_TEST_SUITE_END()



bool init_unit_test() {
  // initializeCUDA(0);
  // TODO: multithread running tests
  thread_buffer.init(available_mem(1));
  std::string logname;
  unsigned seed;
  po::positional_options_description positional;
  po::options_description inputs("Input");
  inputs.add_options()("seed", po::value<unsigned>(&p_args.seed),
      "seed for random number generator")("n_iters",
      po::value<unsigned>(&p_args.n_iters),
      "number of iterations to repeat relevant tests")("log",
      po::value<std::string>(&logname), "specify logfile, default is test.log");
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
  bool help = false;
  po::positional_options_description positional;
  po::options_description inputs("Input");
  inputs.add_options()("seed", po::value<unsigned>(&_dumvar1),
      "seed for random number generator")("n_iters",
      po::value<unsigned>(&_dumvar2),
      "number of iterations to repeat relevant tests")("log",
      po::value<std::string>(&_dumvar3),
      "specify logfile, default is test.log");
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
                ^ po::command_line_style::allow_guessing).positional(positional).allow_unregistered().run(),
        vm);
    notify(vm);
  } catch (po::error& e) {
  }


  if (help) std::cout << desc_simple << '\n';

  return ut::unit_test_main(&init_unit_test, argc, argv);
}
