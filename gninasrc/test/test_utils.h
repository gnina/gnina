#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <string>
#include "tee.h"
#include "gpucode.h"
#include "parsed_args.h"
#include "atom_constants.h"

extern parsed_args p_args;

//TODO: doesn't explicitly prevent/check atoms from overlapping, which could
//theoretically lead to runtime errors later
void make_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, 
             std::mt19937 engine,
             size_t natoms, size_t min_atoms, size_t max_atoms, 
             float max_x, float max_y, float max_z) {

    if (!natoms) {
    //if not provided, randomly generate the number of atoms
        std::uniform_int_distribution<int> natoms_dist(min_atoms, max_atoms+1);
        natoms = natoms_dist(engine);
    }

    //randomly seed reasonable-ish coordinates and types
    //TODO: get charge from type?
    std::uniform_real_distribution<float> coords_dists[3];
    coords_dists[0] = std::uniform_real_distribution<float>(-25, std::nextafter(max_x, FLT_MAX));
    coords_dists[1] = std::uniform_real_distribution<float>(-25, std::nextafter(max_y, FLT_MAX));
    coords_dists[2] = std::uniform_real_distribution<float>(-25, std::nextafter(max_z, FLT_MAX));
    std::uniform_int_distribution<int> charge_dist(-2, 3);
    std::uniform_int_distribution<int> type_dist(0, smina_atom_type::NumTypes-1);

    //set up vector of atoms as well as types
    for (size_t i=0; i<natoms; ++i) {
        atom_params atom;
        atom.charge = charge_dist(engine);
        for (size_t j=0; j<3; ++j) 
            atom.coords[j] = coords_dists[j](engine);
        atoms.push_back(atom);
        atoms[i].charge = charge_dist(engine);
        types.push_back(static_cast<smt>(type_dist(engine)));
    }
}

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

//pretty print tree info for logging
void print_tree(atom_params* atoms, unsigned coords_size, 
        std::vector<unsigned>& ranges, tee& log) {
    unsigned node = 0;
    for (size_t i=0; i<coords_size, ++i) {
        for (size_t j=0; j<ranges.size(); ++j) {
            if (i < ranges[j]) {
                node = j;
                break;
            }
        }
        log << "atom" << i << " " << "node" << node << " " << atoms[i].coords[0] << 
            " " << atoms[i].coords[1] << " " << atoms[i].coords[2] << "\n";
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
