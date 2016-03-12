#include "non_cache_gpu.h"
#include "loop_timer.h"
#include "gpu_math.h"

force_energy_tup::force_energy_tup(void){};

non_cache_gpu::non_cache_gpu(szv_grid_cache& gcache,
                             const grid_dims& gd_,
                             const precalculate_gpu* p_,
                             fl slope_) :
	non_cache(gcache, gd_, p_, slope_)
{
    const model& m = gcache.getModel();
    info.cutoff_sq = p->cutoff_sqr();

    unsigned nlig_atoms = m.num_movable_atoms();
    info.nlig_atoms = nlig_atoms;
    //allocate memory for positions, partial charges, and atom types of movable atoms
    cudaMalloc(&info.lig_atoms, sizeof(atom_params[nlig_atoms]));
    cudaMalloc(&info.lig_penalties, sizeof(force_energy_tup[nlig_atoms]));
    cudaMalloc(&info.result, sizeof(force_energy_tup[nlig_atoms]));
    cudaMalloc(&info.types, sizeof(unsigned[nlig_atoms]));

    //initialize atom types and partial charges
    std::vector<unsigned> htypes(nlig_atoms);

    VINA_FOR(i, nlig_atoms)
    {
        htypes[i] = m.atoms[i].get();
    }
    cudaMemcpy(info.types, &htypes[0], sizeof(unsigned[nlig_atoms]),
               cudaMemcpyHostToDevice);

    info.gridbegins = float3(gd[0].begin, gd[1].begin, gd[2].begin);
    info.gridends = float3(gd[0].end, gd[1].end, gd[2].end);

    assert(gd[0].n > 0);
    assert(gd[1].n > 0);
    assert(gd[2].n > 0);

    //figure out all possibly relevant receptor atoms
    szv recatomids;
    gcache.compute_relevant(gd_, recatomids);
    unsigned nrec_atoms = recatomids.size();
    info.nrec_atoms = nrec_atoms;

    //allocate memory for positions, atom types, and partial charges of all
    //possibly relevant receptor atoms
    cudaMalloc(&info.rec_atoms, sizeof(atom_params[nrec_atoms]));
    cudaMalloc(&info.rectypes, sizeof(unsigned[nrec_atoms]));

    //initialize
    std::vector<atom_params> hrec_atoms(nrec_atoms);
    std::vector<unsigned> hrectypes(nrec_atoms);
    for (unsigned i = 0; i < nrec_atoms; i++)
    {
        unsigned index = recatomids[i];
        const vec& c = m.grid_atoms[index].coords;
        atom_params *a = &hrec_atoms[i];
        a->coords = c;
        a->charge = m.grid_atoms[index].charge;
        
        hrectypes[i] = m.grid_atoms[index].get();
    }
    cudaMemcpy(info.rec_atoms, &hrec_atoms[0], sizeof(atom_params[nrec_atoms]),
               cudaMemcpyHostToDevice);
    cudaMemcpy(info.rectypes, &hrectypes[0], sizeof(unsigned[nrec_atoms]),
               cudaMemcpyHostToDevice);

    info.ntypes = p_->num_types();
    info.splineInfo = p_->getDeviceData();
}

non_cache_gpu::~non_cache_gpu()
{
    //deallocate device memory
    cudaFree(info.lig_atoms);
    cudaFree(info.lig_penalties);
    cudaFree(info.result);
    cudaFree(info.types);
    
    cudaFree(info.rec_atoms);
    cudaFree(info.rectypes);
}

fl non_cache_gpu::eval(const model& m, fl v) const
{
    abort(); //not implemented
}

//evaluate the model on the gpu, v is the curl amount
//sets m.minus_forces and returns total energy
fl non_cache_gpu::eval_deriv(model& m, fl v, const grid& user_grid) const
{
    static loop_timer t;
    t.resume();
    
    //clear energies
    if(user_grid.initialized())
    {
        std::cerr << "usergrid not supported in gpu code yet\n";
        exit(-1);
    }
    
    unsigned nlig_atoms = m.num_movable_atoms();
    cudaMemset(info.result, 0, sizeof(force_energy_tup[nlig_atoms]));
    
    /* TODO:charges */
    atom_params hlig_atoms[nlig_atoms];

    //update coordinates
    for (unsigned i = 0; i < nlig_atoms; i++)
    {
        const vec& c = m.coords[i];
        atom_params *a = &hlig_atoms[i];
        for (unsigned j = 0; j < 3; j++)
        {
            get(a->coords, j) = c[j];
        }
    }
    cudaMemcpy(info.lig_atoms, hlig_atoms, sizeof(hlig_atoms), cudaMemcpyHostToDevice);

    //this will calculate the per-atom energies and forces; curl ignored
    double e = single_point_calc(&info, info.result, slope,
                                 info.nlig_atoms, info.nrec_atoms, v);

    //get forces
    force_energy_tup out[nlig_atoms];
    cudaMemcpy(out, info.result, sizeof(out), cudaMemcpyDeviceToHost);

    for(unsigned i = 0; i < nlig_atoms; i++) {
        force_energy_tup *t = &out[i];
        for(unsigned j = 0; j < 3; j++) {
            m.minus_forces[i][j] = get(t->minus_force, j);
        }
    }

    t.stop();

    return e;
}
