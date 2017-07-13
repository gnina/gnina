#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <string>
#include "tee.h"
#include "gpucode.h"
#include "parsed_args.h"
#include "atom_constants.h"

extern parsed_args p_args;

//pretty print molecule info for logging
void print_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, tee& log) {
    std::string pad = "    ";
    log << "\n";
    for (size_t i=0; i<atoms.size(); ++i) {
        log << i << pad << types[i] << " " << atoms[i].coords[0] << " " << atoms[i].coords[1]
            << " " << atoms[i].coords[2] << pad << atoms[i].charge << "\n";
    }
    log << "\n";
}

//loop boost test case for energy/force calculations
void boost_loop_test(void (*func)())
{
    p_args.iter_count = 0;
    for (auto& param : p_args.params) {
        p_args.seed = param;
        func();
        p_args.iter_count++;
    }
}

#endif
