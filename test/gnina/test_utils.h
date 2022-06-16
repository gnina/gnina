#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <boost/math/quaternion.hpp>
#include <string>

#include "caffe/util/rng.hpp"

#include "gpucode.h"
#include "parsed_args.h"
#include "atom_constants.h"
#include "device_buffer.h"

extern parsed_args p_args;
extern bool run_on_gpu;
extern int cuda_dev_id;
typedef boost::math::quaternion<float> quaternion;

//TODO: doesn't explicitly prevent/check atoms from overlapping, which could
//theoretically lead to runtime errors later
template<typename atomT> void make_mol(std::vector<atom_params>& atoms, std::vector<atomT>& types,
    std::mt19937 engine, size_t natoms = 0, size_t min_atoms = 200,
    size_t max_atoms = 500, float max_x = 25, float max_y = 25,
    float max_z = 25);

//pretty print molecule info for logging
template<typename atomT>
void print_mol(std::vector<atom_params>& atoms,
    std::vector<atomT>& types, tee& log);

//pretty print tree info for logging
inline void print_tree(atom_params* atoms, unsigned coords_size, tee& log) {
  for (size_t i = 0; i < coords_size; ++i) {
    log << "atom" << i << " " << atoms[i].coords[0] << " " << atoms[i].coords[1]
        << " " << atoms[i].coords[2] << "\n";
  }
  log << "\n";
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
