/*
 * test_utils.cpp
 *
 *  Created on: Jun 16, 2022
 *      Author: dkoes
 */
#include "test_utils.h"


//TODO: doesn't explicitly prevent/check atoms from overlapping, which could
//theoretically lead to runtime errors later
template<typename atomT>
void make_mol(std::vector<atom_params>& atoms, std::vector<atomT>& types,
    std::mt19937 engine, size_t natoms, size_t min_atoms,
    size_t max_atoms, float max_x, float max_y, float max_z) {

  if (!natoms) {
    //if not provided, randomly generate the number of atoms
    std::uniform_int_distribution<int> natoms_dist(min_atoms, max_atoms + 1);
    natoms = natoms_dist(engine);
  }
  //randomly seed reasonable-ish coordinates and types
  //TODO: get charge from type?
  std::uniform_real_distribution<float> coords_dists[3];
  coords_dists[0] = std::uniform_real_distribution<float>(-max_x,
      std::nextafter(max_x, FLT_MAX));
  coords_dists[1] = std::uniform_real_distribution<float>(-max_y,
      std::nextafter(max_y, FLT_MAX));
  coords_dists[2] = std::uniform_real_distribution<float>(-max_z,
      std::nextafter(max_z, FLT_MAX));
  std::uniform_int_distribution<int> charge_dist(-2, 3);
  std::uniform_int_distribution<int> type_dist(0, atomT::NumTypes - 1);

  //set up vector of atoms as well as types
  for (size_t i = 0; i < natoms; ++i) {
    atom_params atom;
    atom.charge = charge_dist(engine);
    for (size_t j = 0; j < 3; ++j)
      atom.coords[j] = coords_dists[j](engine);
    atoms.push_back(atom);
    atoms[i].charge = charge_dist(engine);
    types.push_back(static_cast<atomT>(type_dist(engine)));
  }
}

template
void make_mol(std::vector<atom_params>& atoms, std::vector<smt>& types,
    std::mt19937 engine, size_t natoms, size_t min_atoms,
    size_t max_atoms, float max_x, float max_y, float max_z);

//pretty print molecule info for logging
template<typename atomT>
void print_mol(std::vector<atom_params>& atoms,
    std::vector<atomT>& types, tee& log) {
  std::string pad = "    ";
  log << "\n";
  for (size_t i = 0; i < atoms.size(); ++i) {
    log << i << pad << types[i] << " " << atoms[i].coords[0] << " "
        << atoms[i].coords[1] << " " << atoms[i].coords[2] << pad
        << atoms[i].charge << "\n";
  }
  log << "\n";
}

template
void print_mol(std::vector<atom_params>& atoms,
    std::vector<smt>& types, tee& log);
