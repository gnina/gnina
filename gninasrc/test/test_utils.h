#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <string>
#include "tee.h"
#include "parsed_args.h"

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
void boost_multi_energy_test(model (*func)(fl, fl, std::vector<force_energy_tup>))
{
    //maybe don't need to manually manage memory for prec/non_cache/model
    //objects in test functions. if true, could exploit copy elision or move
    //semantics to get minus_forces owned by the model without copying
    p_args.iter_count = 0;
    for (auto& param : p_args.params) {
        p_args.seed = param;
        fl c_out;
        fl g_out;
        std::vector<force_energy_tup> g_forces;
        std::vector<vec> c_forces = func(c_out, g_out, g_forces);
        BOOST_CHECK_CLOSE(c_out, g_out, 0.0001);
        //improved boost collections support in more recent versions will
        //improve concision for this check
        for (size_t i=0; i<c_forces.size(); ++i)
            BOOST_CHECK_CLOSE(c_forces[i], g_forces[i], 0.0001);
        p_args.iter_count++;
    }
}

#endif
