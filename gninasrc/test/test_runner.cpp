#include <boost/program_options.hpp>
#include <iostream>
#include "parsed_args.h"
#include "test_gpucode.h"
#include "tee.h"
#define N_ITERS 1
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

namespace ut = boost::unit_test;
namespace po = boost::program_options;

parsed_args p_args;

BOOST_AUTO_TEST_SUITE(gpu_code)

// BOOST_DATA_TEST_CASE(interaction_energy, ut::data::random(0, INT_MAX) 
        // ^ ut::data::xrange(p_args.many_iters * N_ITERS + 1), random_seed, index) {
    // if (p_args.many_iters)
        // p_args.seed = random_seed;
BOOST_AUTO_TEST_CASE(interaction_energy) {
    fl c_out;
    fl g_out;
    std::cout << p_args.seed << std::endl;
    test_interaction_energy(c_out, g_out);
    BOOST_CHECK_CLOSE(c_out, g_out, 0.0001);
}

BOOST_AUTO_TEST_SUITE_END()

bool init_unit_test()
{
    std::string logname;
    unsigned seed;
    po::positional_options_description positional;
    po::options_description inputs("Input");
    inputs.add_options()
        ("seed,s", po::value<unsigned>(&p_args.seed), "seed for random number generator")
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
        std::cout << p_args.seed << std::endl;
    }
    if (!vm.count("log"))
        logname = "test.log";

    p_args.log.init(logname);
    return true;
}

int main(int argc, char* argv[]) {
    unsigned dummy1;
    std::string dummy2;
    bool dummy3 = true;
    po::positional_options_description positional;
    po::options_description _inputs("Input");
    _inputs.add_options()
        ("seed,s", po::value<unsigned>(&dummy1), "seed for random number generator")
        ("log", po::value<std::string>(&dummy2), "specify logfile, default is test.log");
    po::options_description _info("Information");
    _info.add_options()
        ("help", po::bool_switch(&dummy3), "print usage information");
    po::options_description desc, desc_simple;
    desc.add(_inputs).add(_info);
    desc_simple.add(_inputs).add(_info);
    std::vector<std::string> raw_input(argv, argv + argc);
    for (auto& arg : raw_input) {
        if (arg.find("help") != std::string::npos) {
            	std::cout << desc_simple << '\n';
            	return true;
            }
    }
    return ut::unit_test_main( &init_unit_test, argc, argv );
}
