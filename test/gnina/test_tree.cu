#include <random>
#include "model.h"
#include "test_tree.h"
#include "test_utils.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

void make_tree(model* m) {
  p_args.log << "Tree Set Conf Test \n";
  p_args.log << "Using random seed: " << p_args.seed << "\n";
  p_args.log << "Iteration " << p_args.iter_count;
  p_args.log.endl();
  //set up c++11 random number engine
  std::mt19937 engine(p_args.seed);

  //set up fake mol
  std::vector<atom_params> mol_atoms;
  std::vector<smt> mol_types;
  make_mol(mol_atoms, mol_types, engine, 0);
  std::uniform_int_distribution<unsigned> natoms_dist(1, mol_atoms.size());
  //subdivide into disjoint atom_ranges; ranges.size() == num_nodes
  unsigned remaining_range = mol_atoms.size() - 1;
  std::vector<unsigned> ranges;
  while (remaining_range) {
    unsigned natoms = natoms_dist(engine);
    natoms = natoms > remaining_range ? (natoms % remaining_range) + 1 : natoms;
    unsigned next_val =
        ranges.size() > 0 ? ranges[ranges.size() - 1] + natoms : natoms;
    ranges.push_back(next_val);
    remaining_range -= natoms;
  }

  //set up model
  m->m_num_movable_atoms = mol_atoms.size();
  m->minus_forces = std::vector<vec>(m->num_movable_atoms());

  for (size_t i = 0; i < mol_atoms.size(); ++i) {
    m->coords.push_back(*(vec*) &mol_atoms[i]);
    m->atoms.push_back(atom());
    m->atoms[i].sm = mol_types[i];
    m->atoms[i].charge = mol_atoms[i].charge;
    m->atoms[i].coords = *(vec*) &mol_atoms[i];
  }

  //set up ligand (derived from CPU heterotree<rigid_body> and atom_range)
  rigid_body root(*(vec*) (&mol_atoms[0]), static_cast<sz>(0),
      static_cast<sz>(ranges[0]) + 1);
  flexible_body flex(root);
  m->ligands.push_back(ligand(flex, ranges.size() - 1));
  ligand& lig = m->ligands[0];
  frame parent = root;
  //TODO: built in a bunch of assumptions here - enough to make this a bad test?
  for (size_t i = 1; i < ranges.size(); i++) {
    unsigned p_node = ranges[i - 1];
    unsigned begin = ranges[i - 1] + 1;
    unsigned end = ranges[i] + 1;
    segment next(*(vec*) (&mol_atoms[begin]), begin, end,
        *(vec*) (&mol_atoms[p_node]), parent);
    parent = next;
    lig.children.push_back(next);
  }

  m->initialize_gpu();
}

__global__
void increment_kernel(conf_gpu c, const change_gpu g, fl factor,
    gpu_data* gdata) {
  c.increment(g, factor, gdata);
}

__global__
void set_conf_kernel(gpu_data* gdata, const conf_gpu c) {
  gdata->treegpu->set_conf(gdata->atom_coords, (vec*) gdata->coords, &c);
}

__global__
void derivatives_kernel(tree_gpu *t, const vec * coords, const vec* forces,
    change_gpu c) {
  t->derivative(coords, forces, &c);
}

void test_set_conf() {
  p_args.log << "Set Coords Test\n";
  p_args.log << "Using random seed; " << p_args.seed << "\n";
  p_args.log << "Iteration " << p_args.iter_count;
  p_args.log.endl();
  model* m = new model;
  make_tree(m);
  conf x_cpu = m->get_initial_conf(false);
  conf_gpu x_gpu(x_cpu, m->gdata, thread_buffer);
  change g_cpu(m->get_size(), false);
  fl factor = 1;

  //function to be tested takes an incremented conf object and updates the Cartesian
  //atom coords accordingly; start by generating randomly incremented conf
  std::mt19937 engine(p_args.seed);
  std::uniform_real_distribution<float> change_dist(-10, 10);
  for (size_t i = 0; i < 3; ++i) {
    g_cpu.ligands[0].rigid.position[i] = change_dist(engine);
    g_cpu.ligands[0].rigid.orientation[i] = change_dist(engine);
  }

  for (auto& torsion : g_cpu.ligands[0].torsions)
    torsion = change_dist(engine);

  change_gpu g_gpu(g_cpu, m->gdata, thread_buffer);
  gpu_data* gpu_gdata;
  CUDA_CHECK_GNINA(cudaMalloc(&gpu_gdata, sizeof(gpu_data)));
  CUDA_CHECK_GNINA(
      cudaMemcpy(gpu_gdata, &m->gdata, sizeof(gpu_data),
          cudaMemcpyHostToDevice));
  x_cpu.increment(g_cpu, factor);
  increment_kernel<<<1, x_gpu.n>>>(x_gpu, g_gpu, factor, gpu_gdata);

  //now set the coords
  m->set(x_cpu);
  set_conf_kernel<<<1, m->gdata.coords_size>>>(gpu_gdata, x_gpu);

  //log the tree and check correctness
  std::vector<atom_params> gpu_coords;
  gpu_coords.resize(m->gdata.coords_size);
  CUDA_CHECK_GNINA(
      cudaMemcpy(&gpu_coords[0], m->gdata.coords,
          sizeof(gpu_coords[0]) * m->gdata.coords_size,
          cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();
  p_args.log << "CPU tree \n";
  print_tree((atom_params*) (&m->coords[0]), m->coords.size(), p_args.log);
  p_args.log << "GPU tree \n";
  print_tree(&gpu_coords[0], m->gdata.coords_size, p_args.log);

  for (size_t i = 0; i < m->coords.size(); ++i)
    for (size_t j = 0; j < 3; ++j)
      BOOST_REQUIRE_SMALL(m->coords[i][j] - gpu_coords[i].coords[j],
          (float )0.01);

  delete m;
}

void test_derivative() {
  model* m = new model;
  make_tree(m);
  //for this one we need to generate forces and provide a change object to be set
  change g_cpu(m->get_size(), false);
  change_gpu g_gpu(g_cpu, m->gdata, thread_buffer);
  std::mt19937 engine(p_args.seed);
  std::uniform_real_distribution<float> forces_dist(-10, 10);
  for (size_t i = 0; i < m->coords.size(); ++i)
    for (size_t j = 0; j < 3; ++j)
      m->minus_forces[i][j] = forces_dist(engine);

  CUDA_CHECK_GNINA(
      cudaMemcpy(m->gdata.minus_forces, &m->minus_forces[0],
          sizeof(m->minus_forces[0]) * m->gdata.forces_size,
          cudaMemcpyHostToDevice));
  m->ligands.derivative(m->coords, m->minus_forces, g_cpu.ligands);
  derivatives_kernel<<<1, m->coords.size()>>>(m->gdata.treegpu,
      (vec*) m->gdata.coords, (vec*) m->gdata.minus_forces, g_gpu);
  std::vector<float> gpu_out;
  gpu_out.resize(g_gpu.n);
  CUDA_CHECK_GNINA(
      cudaMemcpy(&gpu_out[0], g_gpu.values, sizeof(gpu_out[0]) * g_gpu.n,
          cudaMemcpyDeviceToHost));
  cudaDeviceSynchronize();

  //log change and check results
  p_args.log << "Derivative Test\n";
  p_args.log << "Using random seed; " << p_args.seed << "\n";
  p_args.log << "Iteration " << p_args.iter_count;
  p_args.log.endl();
  for (size_t i = 0; i < g_gpu.n; ++i) {
    if (i < 3) {
      BOOST_REQUIRE_SMALL(g_cpu.ligands[0].rigid.position[i] - gpu_out[i],
          (float )0.01);
      p_args.log << "CPU " << g_cpu.ligands[0].rigid.position[i] << " GPU "
          << gpu_out[i] << "\n";
    } else
      if (i < 6) {
        BOOST_REQUIRE_SMALL(
            g_cpu.ligands[0].rigid.orientation[i - 3] - gpu_out[i],
            (float )0.01);
        p_args.log << "CPU " << g_cpu.ligands[0].rigid.orientation[i - 3]
            << " GPU " << gpu_out[i] << "\n";
      } else {
        BOOST_REQUIRE_SMALL(g_cpu.ligands[0].torsions[i - 6] - gpu_out[i],
            (float )0.01);
        p_args.log << "CPU " << g_cpu.ligands[0].torsions[i - 6] << " GPU "
            << gpu_out[i] << "\n";
      }
  }

  delete m;
}
