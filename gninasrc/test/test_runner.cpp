#include <cuda_runtime.h>
#include <stdio.h>
#include <boost/program_options.hpp>
#include <cmath>
#include <random>
#include "common.h"
#include "gpu_debug.h"
#include "tee.h"
#include "non_cache.h"
#include "non_cache_gpu.h"
#include "model.h"
#include "curl.h"
#include "weighted_terms.h"
#include "custom_terms.h"
#include "precalculate_gpu.h"
#include "gpucode.h"
//TODO: doesn't explicitly prevent/check atoms from overlapping, which could
//theoretically lead to runtime errors later
void make_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, 
             std::mt19937 engine,
             size_t natoms=0, size_t min_atoms=1, size_t max_atoms=200, 
             float max_x=25, float max_y=25, float max_z=25) {

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

int main(int argc, char* argv[]) {
    //TODO: include progress bar?
    //set up program options
    std::string logname;
    unsigned seed;
    namespace po = boost::program_options;
	po::positional_options_description positional; // remains empty
    po::options_description inputs("Input");
    inputs.add_options()
        ("seed,s", po::value<unsigned>(&seed), "seed for random number generator")
        ("log", po::value<std::string>(&logname), "specify logfile, default is test.log");
    po::options_description desc, desc_simple;
    desc.add(inputs);
    desc_simple.add(inputs);
    po::variables_map vm;
    try
    {
        po::store(
                po::command_line_parser(argc, argv).options(desc)
                .style(
                    po::command_line_style::default_style
                    ^ po::command_line_style::allow_guessing)
                .positional(positional).run(), vm);
        notify(vm);
    } catch (po::error& e)
    {
		std::cerr << "Command line parse error: " << e.what() << '\n'
				<< "\nCorrect usage:\n" << desc_simple << '\n';
		return 1;
	}
    if (!vm.count("seed"))
        seed = std::random_device()();
    if (!vm.count("logname"))
        logname = "test.log";
    //set up c++11 random number engine
    std::mt19937 engine(seed);

    //set up logging
    bool quiet = true;
    tee log(quiet);
    log.init(logname);

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
    const fl granularity = 0.375;
    const vec v = vec(10, 10, 10);
    const fl slope = 10;

    weighted_terms wt(&t, t.weights());

    //set up splines
    const precalculate_gpu* gprec = new precalculate_gpu(wt, approx_factor);
    const precalculate_splines* prec = new precalculate_splines(wt, approx_factor);

    //set up fake lig and rec for testing
    std::vector<atom_params> lig_atoms;
    std::vector<smt> lig_types;
    fl max_x = -HUGE_VALF, max_y = -HUGE_VALF, max_z = -HUGE_VALF;
    fl min_x = HUGE_VALF, min_y = HUGE_VALF, min_z = HUGE_VALF;
    make_mol(lig_atoms, lig_types, engine, 10);

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

    for (size_t i; i < 3; ++i) {
        gd[i].n = sz(std::ceil(span[i] / granularity));
        fl real_span = granularity * gd[i].n;
        gd[i].begin = center[i] - real_span / 2;
        gd[i].end = gd[i].begin + real_span;
    }

    //make rec
    std::vector<atom_params> rec_atoms;
    std::vector<smt> rec_types;
    const float cutoff_sqr = prec->cutoff_sqr();
    const float cutoff = std::sqrt(cutoff_sqr);
    make_mol(rec_atoms, rec_types, engine, 1377, 10, 2500, max_x + cutoff, max_y + 
            cutoff, max_z + cutoff);

    //manually initialize model object
    model* m = new model;
    m->m_num_movable_atoms = lig_atoms.size();
    m->minus_forces = std::vector<vec>(m->m_num_movable_atoms);

    for (size_t i=0; i <lig_atoms.size(); ++i) {
        m->coords.push_back(*(vec*)&lig_atoms[i]);
        m->atoms.push_back(atom());
        m->atoms[i].sm = lig_types[i];
        m->atoms[i].charge = lig_atoms[i].charge;
        m->atoms[i].coords = *(vec*)&lig_atoms[i];
    }

    for (size_t i=0; i<rec_atoms.size(); ++i) {
        m->grid_atoms.push_back(atom());
        m->grid_atoms[i].sm = rec_types[i];
        m->grid_atoms[i].charge = rec_atoms[i].charge;
        m->grid_atoms[i].coords = *(vec*)&rec_atoms[i];
    }

    szv_grid_cache gridcache(*m, cutoff_sqr);
    non_cache* nc = new non_cache(gridcache, gd, prec);
    non_cache_gpu* nc_gpu = new non_cache_gpu(gridcache, gd, gprec, slope);
    m->initialize_gpu();
    grid user_grid;
    gpu_data& gdat = m->gdata;
    cudaMemset(gdat.minus_forces, 0, m->minus_forces.size()*sizeof(gdat.minus_forces[0]));
    std::cout << rec_atoms.size() << "\n";
    std::cout << lig_atoms.size() << "\n";

    //get intermolecular energy, check agreement
    float g_out = single_point_calc(nc_gpu->info, gdat.coords, gdat.minus_forces, v[0]);
    float c_out = nc->eval_deriv(*m, v[0], user_grid);
    vec g_forces[m->minus_forces.size()];
    cudaMemcpy(g_forces, gdat.minus_forces, m->minus_forces.size()*sizeof(gdat.minus_forces[0]), cudaMemcpyDeviceToHost);
    pretty_print_vec_array(g_forces, m->minus_forces.size(), "gpu forces");
    pretty_print_vec_array(&m->minus_forces[0], m->minus_forces.size(), "cpu forces");

    std::cout << g_out << "\n";
    std::cout << c_out << "\n";

    m->deallocate_gpu();
    delete gprec;
    delete prec;
    delete nc;
    delete nc_gpu;
    delete m;
}
