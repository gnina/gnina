#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include <string>
#include "tee.h"

//pretty print molecule info for logging
void print_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, tee& log) {
    std::string pad = "    ";
    log << "\n";
    for (unsigned i=0; i<atoms.size(); ++i) {
        log << i << pad << types[i] << " " << atoms[i].coords[0] << " " << atoms[i].coords[1]
            << " " << atoms[i].coords[2] << pad << atoms[i].charge << "\n";
    }
    log << "\n";
}

#endif
