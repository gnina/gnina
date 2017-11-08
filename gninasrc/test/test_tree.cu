#include "model.h"
#include "parsed_args.h"
#include <random>
#include "test_tree.h"
#include "test_utils.h"
#include "tree.h"
#include "tree_gpu.h"

extern parsed_args p_args;
thread_local float_buffer buffer;

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
    unsigned remaining_range = mol_atoms.size();
    std::vector<unsigned> ranges;
    while (remaining_range) {
        unsigned natoms = natoms_dist(engine);
        natoms = natoms > remaining_range ? (natoms % remaining_range) + 1: natoms;
        ranges[i] = i > 0 ? ranges[i-1] + natoms : natoms;
        remaining_range -= ranges[i];
    }

    //set up model
    m->m_num_movable_atoms = mol_atoms.size();
    m->minus_forces = std::vector<vec>(m->num_movable_atoms);

    for (size_t i=0; i <mol_atoms.size(); ++i) {
        m->coords.push_back(*(vec*)&mol_atoms[i]);
        m->atoms.push_back(atom());
        m->atoms[i].sm = mol_types[i];
        m->atoms[i].charge = mol_atoms[i].charge;
        m->atoms[i].coords = *(vec*)&mol_atoms[i];
    }

    //set up ligand (CPU heterotree<rigid_body>)
    rigid_body root(mol_atoms[0].coords, 0, ranges[0]);
    m->ligands[0] = ligand(root);
    ligand& lig = m->ligands[0];
    frame parent = root;
    //TODO: built in a bunch of assumptions here - enough to make this a bad test?
    for (size_t i=1; i < ranges.size(); i++) {
        unsigned p_node = range[i-1];
        unsigned begin = ranges[i-1] + 1;
        unsigned end = ranges[i];
        segment next(mol_atoms[begin].coords, begin, end, mol_atoms[p_node].coords, parent);
        parent = next;
        lig.children.push_back(next);
    }

    m->initialize_gpu();
}

__global__ 
void increment_kernel(conf_gpu c, const change_gpu g, fl factor, gpu_data* gdata) { 
    c.increment(g, factor, gdata);
}

__global__ 
void set_conf_kernel(gpu_data* gdata, const conf_gpu c) {
    gdata->treegpu->set_conf(gdata.atom_coords, gdata.coords, &c);
}

void test_set_conf() {
    model* m = new model;
    make_tree(m);
    conf x_cpu = m->get_initial_conf();
    conf_gpu x_gpu(c_cpu, m->gdata, buffer);
    change g_cpu(m->get_size());
    fl factor = 1;

    //function to be tested takes an incremented conf object and updates the Cartesian
    //atom coords accordingly; start by generating randomly incremented conf
    std::mt19937 engine(p_args.seed);
    std::uniform_real_distribution<float> change_dist(-10, 10);
    for (size_t i=0; i<3; ++i) {
        g_cpu.ligands[0].rigid_change.position[i] = change_dist(engine);
        g_cpu.ligands[0].rigid_change.orientation[i] = change_dist(engine);
    }

    for (auto& torsion : g_cpu.ligands[0].torsions)
        torsion = change_dist(engine);

    change_gpu g_gpu(g_cpu, m->gdata, buffer);
    gpu_data* gpu_gdata;
    CUDA_CHECK_GNINA(cudaMalloc(&gpu_gdata, sizeof(gpu_data)));
    CUDA_CHECK_GNINA(cudaMemcpy(&gpu_gdata, m->gdata, sizeof(gpu_data), cudaMemcpyHostToDevice));
    x_cpu.increment(g_cpu, factor);
    increment_kernel<<<1, x_gpu.n>>>(x_gpu, g_gpu, factor, gpu_gdata);

    //now set the coords
    m->set(x_cpu);
    set_conf_kernel<<<1, mol_atoms.size()>>>(gpu_gdata, x_gpu);

    //log the tree and check correctness
    CUDA_CHECK_GNINA(cudaMemcpy(m->gdata->coords, &gpu_gdata->coords, 
                sizeof(m->gdata->coords) * m->gdata->coords.size, 
                cudaMemcpyDeviceToHost));
    p_args.log << "CPU tree \n";
    print_tree(static_cast<atom_params>(&m->coords[0]), m->coords.size(), ranges, p_args.log);
    p_args.log << "GPU tree \n";
    print_tree(m->gdata->coords, m->gdata->coords.size, ranges, p_args.log);

    for (size_t i=0; i<m->coords.size(); ++i) 
        for (size_t j=0; j<3; ++j)
            BOOST_REQUIRE_SMALL(m->coords[i][j] - m->gdata->coords[i].coords[j], 
                    (float)0.01);

    m->deallocate_gpu();
    delete m;
}
