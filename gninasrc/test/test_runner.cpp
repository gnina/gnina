#include <boost/program_options.hpp>
#include <iostream>
#include "parsed_args.h"
#include "test_gpucode.h"
#define N_ITERS 10
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

namespace ut = boost::unit_test;
namespace po = boost::program_options;

parsed_args p_args;

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

bool init_unit_test()
{
    std::string logname;
    unsigned seed;
    po::positional_options_description positional;
    po::options_description inputs("Input");
    inputs.add_options()
        ("seed", po::value<unsigned>(&p_args.seed), "seed for random number generator")
        ("log", po::value<std::string>(&logname), "specify logfile, default is test.log");
    po::options_description desc, desc_simple;
    desc.add(inputs);
    desc_simple.add(inputs);
    po::variables_map vm;
    try {
        po::store(
                po::command_line_parser(ut::framework::master_test_suite().argc, 
                    ut::framework::master_test_suite().argv).options(desc)
                .style(
                    po::command_line_style::default_style
                    ^ po::command_line_style::allow_guessing)
                .positional(positional).run(), vm);
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
    if (!vm.count("log"))
        logname = "test.log";

    p_args.log.init(logname);

    p_args.params = {p_args.seed};
    if (p_args.many_iters) 
        for (size_t i=0; i<N_ITERS; ++i)
            p_args.params.push_back(std::random_device()());

    return true;
}

int main(int argc, char* argv[]) {
    //The following exists so that passing --help prints both the UTF help and our
    //specific program options - I can't see any way of doing this without
    //parsing the args twice, because UTF chomps args it thinks it owns
    unsigned _dumvar1;
    std::string _dumvar2;
    bool help = false;
    po::positional_options_description positional;
    po::options_description inputs("Input");
    inputs.add_options()
        ("seed,s", po::value<unsigned>(&_dumvar1), "seed for random number generator")
        ("log", po::value<std::string>(&_dumvar2), "specify logfile, default is test.log");
    po::options_description info("Information");
    info.add_options()
        ("help", po::bool_switch(&help), "print usage information");
    po::options_description desc, desc_simple;
    desc.add(inputs).add(info);
    desc_simple.add(inputs).add(info);
    po::variables_map vm;
    try {
        po::store(
                po::command_line_parser(argc, argv).options(desc)
                .style(
                    po::command_line_style::default_style
                    ^ po::command_line_style::allow_guessing)
                .positional(positional).run(), vm);
        notify(vm);
    } catch (po::error& e) {
    }
    
    if (help)
        std::cout << desc_simple << '\n';

    return ut::unit_test_main( &init_unit_test, argc, argv );
}
