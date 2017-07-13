#ifndef TEST_GPUCODE_H
#define TEST_GPUCODE_H
#include "common.h"
#include "gpucode.h"

void make_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, 
             std::mt19937 engine,
             size_t natoms=0, size_t min_atoms=1, size_t max_atoms=200, 
             float max_x=25, float max_y=25, float max_z=25);

void test_interaction_energy();
void test_eval_intra();

#endif
