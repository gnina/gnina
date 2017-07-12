#include "test_gpucode.h"

extern parsed_args p_args;

//TODO: doesn't explicitly prevent/check atoms from overlapping, which could
//theoretically lead to runtime errors later
void make_mol(std::vector<atom_params>& atoms, std::vector<smt>& types, 
             std::mt19937 engine,
             size_t natoms, size_t min_atoms, size_t max_atoms, 
             float max_x, float max_y, float max_z) {

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

void test_interaction_energy(fl& c_out, fl& g_out) {
    p_args.log << "Using random seed: " << p_args.seed;
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
    const precalculate_gpu* gprec = new precalculate_gpu(wt, approx_factor);
    const precalculate_splines* prec = new precalculate_splines(wt, approx_factor);

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

    for (size_t i=0; i < 3; ++i) {
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
    make_mol(rec_atoms, rec_types, engine, 0, 10, 2500, max_x + cutoff, max_y + 
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

    //make non_cache
    non_cache* nc = new non_cache(gridcache, gd, prec, slope);
    non_cache_gpu* nc_gpu = new non_cache_gpu(gridcache, gd, gprec, slope);

    //set up GPU data
    m->initialize_gpu();
    gpu_data& gdat = m->gdata;
    cudaMemset(gdat.minus_forces, 0, m->minus_forces.size()*sizeof(gdat.minus_forces[0]));

    //get intermolecular energy, check agreement
    g_out = single_point_calc(nc_gpu->info, gdat.coords, gdat.minus_forces, v);
    c_out = nc->eval_deriv(*m, v, user_grid);
    //TODO: check whether the forces match, too
    vec g_forces[m->minus_forces.size()];
    cudaMemcpy(g_forces, gdat.minus_forces, m->minus_forces.size()*sizeof(gdat.minus_forces[0]), cudaMemcpyDeviceToHost);

    //clean up after yourself
    delete nc;
    delete nc_gpu;
    delete gprec;
    delete prec;
    m->deallocate_gpu();
    delete m;
}

void test_eval_intra(fl& c_out, fl& g_out, size_t natoms, size_t min_atoms, 
                     size_t max_atoms) {

    p_args.log << "Using random seed: " << p_args.seed;
    p_args.log.endl();

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
    const fl slope = 10;
    
    weighted_terms wt(&t, t.weights());

    //set up splines
    const precalculate_gpu* gprec = new precalculate_gpu(wt, approx_factor);
    const precalculate_splines* prec = new precalculate_splines(wt, approx_factor);

    //generate the mol
    std::mt19937 engine(p_args.seed);
    std::vector<atom_params> atoms;
    std::vector<smt> types;
    make_mol(atoms, types, engine, natoms, min_atoms, max_atoms);

    //generate pairs vector consisting of every combination that isn't super close
    model* m = new model;
    for (size_t i=0; i<atoms.size(); ++i) {
        for (size_t j=i+1; j<atoms.size(); ++j) {
            float r2 = vec_distance_sqr(*(vec*)&atoms[i], *(vec*)&atoms[j]);
            //TODO? threshold for "closeness" is pretty arbitrary here...
            if (r2 > 4) {
                interacting_pair ip;
                ip.t1 = types[i];
                ip.t2 = types[j];
                ip.a = i;
                ip.b = j;
                m->other_pairs.push_back(ip);
            }
        }
    }

    //set up model
    m->minus_forces = std::vector<vec>(atoms.size());
    for (size_t i=0; i <atoms.size(); ++i) {
        m->coords.push_back(*(vec*)&atoms[i]);
        m->atoms.push_back(atom());
        m->atoms[i].sm = types[i];
        m->atoms[i].charge = atoms[i].charge;
        m->atoms[i].coords = *(vec*)&atoms[i];
    }

    //set up GPU data structures
    gpu_data& gdat = m->gdata;
    m->initialize_gpu();
    cudaMemset(gdat.minus_forces, 0, m->minus_forces.size()*sizeof(gdat.minus_forces[0]));
    GPUNonCacheInfo dinfo;
    dinfo.cutoff_sq = prec->cutoff_sqr();
    dinfo.splineInfo = gprec->deviceData;

    //compute intra energy and compare
    c_out = m->eval_interacting_pairs_deriv(*prec, v, m->other_pairs, m->coords, m->minus_forces);
    g_out = gdat.eval_interacting_pairs_deriv_gpu(dinfo, v, gdat.other_pairs, 
            m->other_pairs.size());

    m->deallocate_gpu();
    delete m;
    delete gprec;
    delete prec;
}
