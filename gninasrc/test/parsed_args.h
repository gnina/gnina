#ifndef PARSED_ARGS_H
#define PARSED_ARGS_H
#include "tee.h"
#include <boost/filesystem.hpp>

struct parsed_args {
    unsigned seed;
    bool many_iters;
    unsigned n_iters;
    unsigned iter_count;
    tee log;
    std::vector<unsigned> params;
    std::string rec;

    parsed_args(bool quiet = true)
        : many_iters(false), iter_count(0), log(quiet) {
      rec = boost::filesystem::current_path().string() + "/10gs_pocket.pdb";
    }
};

#endif
