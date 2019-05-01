#ifndef TEST_UTILS_H
#define TEST_UTILS_H
#include "caffe/util/rng.hpp"
#include "caffe/layers/molgrid_data_layer.hpp"
#include <string>
#include "tee.h"
#include "gpucode.h"
#include <boost/math/quaternion.hpp>
#include "cnn_scorer.h"
#include "parsed_args.h"
#include "atom_constants.h"
#include "atom.h"
#include "device_buffer.h"

extern parsed_args p_args;
extern bool run_on_gpu;
extern int cuda_dev_id;
typedef boost::math::quaternion<float> quaternion;

//TODO: doesn't explicitly prevent/check atoms from overlapping, which could
//theoretically lead to runtime errors later
template<typename atomT>
inline void make_mol(std::vector<atom_params>& atoms, std::vector<atomT>& types,
    std::mt19937 engine, size_t natoms = 0, size_t min_atoms = 1,
    size_t max_atoms = 200, float max_x = 25, float max_y = 25,
    float max_z = 25) {

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

inline void make_mol(atomv& atoms, std::mt19937 engine, std::vector<int>& types, size_t natoms = 0,
    size_t min_atoms = 1, size_t max_atoms = 200, float max_x = 25, float max_y = 25, 
    float max_z = 25) {
  if (!natoms) {
    //if not provided, randomly generate the number of atoms
    std::uniform_int_distribution<int> natoms_dist(min_atoms, max_atoms + 1);
    natoms = natoms_dist(engine);
  }

  //randomly seed reasonable-ish coordinates and types
  std::uniform_real_distribution<float> coords_dists[3];
  coords_dists[0] = std::uniform_real_distribution<float>(-max_x,
      std::nextafter(max_x, FLT_MAX));
  coords_dists[1] = std::uniform_real_distribution<float>(-max_y,
      std::nextafter(max_y, FLT_MAX));
  coords_dists[2] = std::uniform_real_distribution<float>(-max_z,
      std::nextafter(max_z, FLT_MAX));
  std::uniform_int_distribution<int> type_dist(0, types.size()-1);

  //set up vector of atoms as well as types
  for (size_t i = 0; i < natoms; ++i) {
    atoms.push_back(atom());
    for (size_t j = 0; j < 3; ++j)
      atoms[i].coords[j] = coords_dists[j](engine);
    atoms[i].sm = static_cast<smt>(types[type_dist(engine)]);
  }
}

//pretty print molecule info for logging
template<typename atomT>
inline void print_mol(std::vector<atom_params>& atoms,
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

//pretty print tree info for logging
inline void print_tree(atom_params* atoms, unsigned coords_size, tee& log) {
  for (size_t i = 0; i < coords_size; ++i) {
    log << "atom" << i << " " << atoms[i].coords[0] << " " << atoms[i].coords[1]
        << " " << atoms[i].coords[2] << "\n";
  }
  log << "\n";
}

//set up CNN grids from randomly-generated mol
template <typename atomT, typename MGridT, typename GridMakerT>
inline void set_cnn_grids(MGridT* mgrid, 
  GridMakerT& gmaker, std::vector<atom_params>& mol_atoms, std::vector<atomT>& mol_types) {
  //first set up mgrid
  mgrid->batch_transform.resize(1);
  mgrid->batch_transform[0] = caffe::BaseMolGridDataLayer<CNNScorer::Dtype, GridMaker>::mol_transform();
  caffe::BaseMolGridDataLayer<CNNScorer::Dtype, GridMaker>::mol_transform& transform = mgrid->batch_transform[0];
  vec center(0,0,0);
  for (size_t i=0; i<mol_atoms.size(); ++i) {
    atom_params& ainfo = mol_atoms[i];
    transform.mol.atoms.push_back(make_float4(ainfo.coords.x,
                ainfo.coords.y, ainfo.coords.z, ainfo.charge));
    transform.mol.whichGrid.push_back(mol_types[i]);
    transform.mol.gradient.push_back(make_float3(0,0,0));
    center += vec(ainfo.coords.x, ainfo.coords.y, ainfo.coords.z);
  }
  center /= mol_atoms.size();
  transform.center = center;
  transform.Q = quaternion(1,0,0,0);
  transform.set_random_quaternion(caffe::caffe_rng());

  //then set up gridmaker
  bool spherize = false;
  double dim = round(mgrid->dimension/mgrid->resolution)+1;
  gmaker.initialize(mgrid->resolution, mgrid->dimension, mgrid->radiusmultiple, 
          mgrid->binary, spherize);
  gmaker.setCenter(center[0], center[1], center[2]);
}

//loop boost test case for energy/force calculations
inline void boost_loop_test(void (*func)()) {
  p_args.iter_count = 0;
  for (auto& param : p_args.params) {
    p_args.seed = param;
    thread_buffer.reinitialize();
    func();
    p_args.iter_count++;
  }
}

#endif
