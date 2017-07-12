#ifndef TEST_GPUCODE_H
#define TEST_GPUCODE_H
#include <numeric>
#include <cmath>
#include <random>
#include "common.h"
#include "parsed_args.h"
#include "gpu_debug.h"
#include "tee.h"
#include "non_cache.h"
#include "non_cache_gpu.h"
#include "tree_gpu.h"
#include "model.h"
#include "curl.h"
#include "weighted_terms.h"
#include "custom_terms.h"
#include "precalculate_gpu.h"
#include "gpucode.h"

void make_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, 
             std::mt19937 engine,
             size_t natoms=0, size_t min_atoms=1, size_t max_atoms=200, 
             float max_x=25, float max_y=25, float max_z=25);

void test_interaction_energy(fl& c_out, fl& g_out);

void test_eval_intra(fl& c_out, fl& g_out, size_t natoms=0, size_t min_atoms=1, 
                     size_t max_atoms=200);
#endif
