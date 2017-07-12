#ifndef PARSED_ARGS_H
#define PARSED_ARGS_H
#include "tee.h"

struct parsed_args {
    unsigned seed;
    bool many_iters;
    tee log;
    std::vector<unsigned> params;

    parsed_args(bool quiet = true) : many_iters(false), log(quiet) {}
};

#endif
