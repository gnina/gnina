#include <numeric>
#include <cmath>
#include <random>
#include "common.h"
#include "cache_gpu.h"
#include "weighted_terms.h"
#include "custom_terms.h"
#include "precalculate_gpu.h"
#include "szv_grid.h"
#include "test_cache.h"
#include "parsed_args.h"
#include "test_utils.h"
#include <cuda_runtime.h>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

void test_cache_eval_deriv() {
  p_args.log << "Cache Eval Test \n";
  p_args.log << "Using random seed: " << p_args.seed << "\n";
  p_args.log << "Iteration " << p_args.iter_count;
  p_args.log.endl();
  //set up c++11 random number engine
  std::mt19937 engine(p_args.seed);

  //set up scoring function
  custom_terms t;
  t.add("gauss(o=0,_w=0.5,_c=8)", -0.035579);
  t.add("gauss(o=3,_w=2,_c=8)", -0.005156);
  t.add("repulsion(o=0,_c=8)", 0.840245);
  t.add("hydrophobic(g=0.5,_b=1.5,_c=8)", -0.035069);
  t.add("non_dir_h_bond(g=-0.7,_b=0,_c=8)", -0.587439);
  t.add("num_tors_div", 5 * 0.05846 / 0.1 - 1);

  //set up a bunch of constants
  const fl approx_factor = 10;
  const fl v = 10;
  const fl granularity = 0.375;
  const fl slope = 10;

  weighted_terms wt(&t, t.weights());

  //set up splines
  std::unique_ptr<precalculate_gpu> gprec(
      new precalculate_gpu(wt, approx_factor));
  std::unique_ptr<precalculate_splines> prec(
      new precalculate_splines(wt, approx_factor));

  //set up lig
  std::vector<atom_params> lig_atoms;
  std::vector<smt> lig_types;
  fl max_x = -HUGE_VALF, max_y = -HUGE_VALF, max_z = -HUGE_VALF;
  fl min_x = HUGE_VALF, min_y = HUGE_VALF, min_z = HUGE_VALF;
  make_mol(lig_atoms, lig_types, engine, 0);

  //set up grid
  for (auto& atom : lig_atoms) {
    min_x = std::min(min_x, atom.coords[0]);
    min_y = std::min(min_y, atom.coords[1]);
    min_z = std::min(min_z, atom.coords[2]);
    max_x = std::max(max_x, atom.coords[0]);
    max_y = std::max(max_y, atom.coords[1]);
    max_z = std::max(max_z, atom.coords[2]);
  }

  fl center_x = (max_x + min_x) / 2.0;
  fl center_y = (max_y + min_y) / 2.0;
  fl center_z = (max_z + min_z) / 2.0;
  fl size_x = max_x - min_x;
  fl size_y = max_y - min_y;
  fl size_z = max_z - min_z;

  vec span(size_x, size_y, size_z);
  vec center(center_x, center_y, center_z);
  grid_dims gd;

  for (size_t i = 0; i < 3; ++i) {
    gd[i].n = sz(std::ceil(span[i] / granularity));
    fl real_span = granularity * gd[i].n;
    gd[i].begin = center[i] - real_span / 2;
    gd[i].end = gd[i].begin + real_span;
  }
  grid user_grid;

  //set up rec
  std::vector<atom_params> rec_atoms;
  std::vector<smt> rec_types;
  const float cutoff_sqr = prec->cutoff_sqr();
  const float cutoff = std::sqrt(cutoff_sqr);
  make_mol(rec_atoms, rec_types, engine, 0, 10, 500, max_x + cutoff,
      max_y + cutoff, max_z + cutoff);

  //manually initialize model object
  std::unique_ptr<model> m(new model);
  m->m_num_movable_atoms = lig_atoms.size();
  m->minus_forces = std::vector<vec>(m->m_num_movable_atoms);

  for (size_t i = 0; i < lig_atoms.size(); ++i) {
    m->coords.push_back(*(vec*) &lig_atoms[i]);
    m->atoms.push_back(atom());
    m->atoms[i].sm = lig_types[i];
    m->atoms[i].charge = lig_atoms[i].charge;
    m->atoms[i].coords = *(vec*) &lig_atoms[i];
  }

  for (size_t i = 0; i < rec_atoms.size(); ++i) {
    m->grid_atoms.push_back(atom());
    m->grid_atoms[i].sm = rec_types[i];
    m->grid_atoms[i].charge = rec_atoms[i].charge;
    m->grid_atoms[i].coords = *(vec*) &rec_atoms[i];
  }

  szv_grid_cache gridcache(*m, cutoff_sqr);

  //set up cache
  std::unique_ptr<cache> c(new cache("scoring_function_version001", gd, slope));
  std::unique_ptr<cache_gpu> cg(
      new cache_gpu("scoring_function_version001", gd, slope, &(*gprec)));

  std::vector<smt> atom_types_needed;
  m->get_movable_atom_types(atom_types_needed);
  c->populate(*m, *prec, atom_types_needed, user_grid);
  cg->populate(*m, *gprec, atom_types_needed, user_grid);

  //set up GPU data
  m->initialize_gpu();
  gpu_data& gdat = m->gdata;
  cudaMemset(gdat.minus_forces, 0,
      m->minus_forces.size() * sizeof(gdat.minus_forces[0]));

  //get intermolecular energy and forces
  fl c_out = c->eval_deriv(*m, v, user_grid);
  fl g_out = single_point_calc(cg->get_info(), gdat.coords, gdat.minus_forces,
      v);

  vec g_forces[m->minus_forces.size()];
  cudaMemcpy(g_forces, gdat.minus_forces,
      m->minus_forces.size() * sizeof(gdat.minus_forces[0]),
      cudaMemcpyDeviceToHost);

  //log the mols and check results
  print_mol(rec_atoms, rec_types, p_args.log);
  print_mol(lig_atoms, lig_types, p_args.log);
  p_args.log << "CPU energy: " << c_out << " GPU energy: " << g_out << "\n\n";

  BOOST_REQUIRE_SMALL(c_out - g_out, (float )0.01);
  //TODO: enhanced boost collections support in more recent versions will
  //improve concision for this check
  for (size_t i = 0; i < m->minus_forces.size(); ++i)
    for (size_t j = 0; j < 3; ++j)
      BOOST_REQUIRE_SMALL(m->minus_forces[i][j] - g_forces[i][j], (float )0.01);
}
