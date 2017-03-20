#include "non_cache_gpu.h"
#include "loop_timer.h"
#include "gpu_math.h"


non_cache_gpu::non_cache_gpu(szv_grid_cache& gcache,
                             const grid_dims& gd_,
                             const precalculate_gpu* p_,
                             fl slope_) :
	non_cache(gcache, gd_, p_, slope_)
{
  const model& m = gcache.getModel();
  info.cutoff_sq = p->cutoff_sqr();
  info.slope = slope;

  m.print_counts(info.nrec_atoms);

  unsigned nlig_atoms = m.num_movable_atoms();
  info.nlig_atoms = nlig_atoms;
  //allocate memory for positions, partial charges, and atom types of movable atoms
  cudaMalloc(&info.lig_penalties, sizeof(force_energy_tup[nlig_atoms]));
  cudaMalloc(&info.types, sizeof(unsigned[nlig_atoms]));

  //initialize atom types and partial charges
  std::vector<unsigned> htypes(nlig_atoms);

  VINA_FOR(i, nlig_atoms)
  {
    htypes[i] = m.atoms[i].get();
    /* TODO breaking const */
    ((atom_params *) &m.coords[0])[i].charge = m.atoms[i].charge;
    /* lig_atoms_scratch[i].charge = 101010; */
  }
  cudaMemcpy(info.types, &htypes[0], sizeof(unsigned[nlig_atoms]),
             cudaMemcpyHostToDevice);

  info.gridbegins = float3(gd[0].begin, gd[1].begin, gd[2].begin);
  info.gridends = float3(gd[0].end, gd[1].end, gd[2].end);

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
    a->coords.x = c[0];
    a->coords.y = c[1];
    a->coords.z = c[2];
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
  cudaFree(info.lig_penalties);
  cudaFree(info.types);
    
  cudaFree(info.rec_atoms);
  cudaFree(info.rectypes);

  /* print_hits(); */
}

fl non_cache_gpu::eval(const model& m, fl v) const
{
  abort(); //not implemented
}

void non_cache_gpu::setSlope(fl sl) {
    slope = sl;
    info.slope = sl;
}

