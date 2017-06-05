#include <cuda_runtime.h>
#include <stdio.h>
#include <random>
#include "common.h"
#include "curl.h"
#include "weighted_terms.h"
#include "custom_terms.h"
#include "precalculate_gpu.h"
#include "gpucode.h"

//TODO: logging, user-provided random seed

void make_lig(std::vector<atom_params>& lig_atoms,
             std::mt19937 engine,
             size_t nlig_atoms=0, size_t min_ligatoms=1, size_t max_ligatoms=200) {

    if (!nlig_atoms) {
    //if not provided, randomly generate the number of ligand atoms
        std::uniform_int_distribution<int> natoms_dist(min_ligatoms, max_ligatoms+1);
        nlig_atoms = natoms_dist(engine);
    }

    //randomly seed reasonable-ish coordinates and types
    //TODO: get charge from type?
    std::uniform_real_distribution<float> coords_dist(0, std::nextafter(50, FLT_MAX));
    std::uniform_int_distribution<int> charge_dist(-2, 3);
    std::uniform_int_distribution<int> type_dist(0, 28); //what is NumTypes?

    //set up vector of lig atoms
    for (size_t i=0; i<nlig_atoms; ++i) {
        atom_params atom;
        atom.charge = charge_dist(engine);
        for (size_t j=0; j<3; ++j) 
            atom.coords[j] = coords_dist(engine);
        lig_atoms.push_back(atom);
    }
}

void make_dinfo(GPUNonCacheInfo& dinfo, 
        std::mt19937 engine, size_t nlig_atoms, std::vector<unsigned>& lig_types, 
        std::vector<atom_params>& rec_atoms, std::vector<unsigned>& rec_types, 
        size_t nrec_atoms=0, size_t min_recatoms=10, size_t max_recatoms=2500) {

    if (!nrec_atoms) {
        std::uniform_int_distribution<int> natoms_dist(min_recatoms, max_recatoms+1);
        nrec_atoms = natoms_dist(engine);
    }

    std::uniform_real_distribution<float> coords_dist(0, std::nextafter(150, FLT_MAX));
    std::uniform_int_distribution<int> charge_dist(-2, 3);
    std::uniform_int_distribution<int> type_dist(0, 28); //what is NumTypes?

    //set up vector of rec atoms
    for (size_t i=0; i<nrec_atoms; ++i) {
        atom_params atom;
        atom.charge = charge_dist(engine);
        for (size_t j=0; j<3; ++j) 
            atom.coords[j] = coords_dist(engine);
        rec_atoms.push_back(atom);
        rec_types.push_back(type_dist(engine));
    }

    CUDA_CHECK_GNINA(cudaMalloc(&dinfo.rec_atoms, sizeof(atom_params)*rec_atoms.size()));
    CUDA_CHECK_GNINA(cudaMemcpy(dinfo.rec_atoms, &rec_atoms[0], sizeof(atom_params)*rec_atoms.size(), cudaMemcpyHostToDevice));
    CUDA_CHECK_GNINA(cudaMalloc(&dinfo.rectypes, sizeof(unsigned)*rec_types.size()));
    CUDA_CHECK_GNINA(cudaMemcpy(dinfo.rectypes, &rec_types[0], sizeof(unsigned)*rec_types.size(), cudaMemcpyHostToDevice));

    //TODO: test receptor flexibility too
    dinfo.num_movable_atoms = nlig_atoms;
    dinfo.nrec_atoms = nrec_atoms;
    dinfo.cutoff_sq = 100;
    dinfo.slope = 10;
    
    for (size_t i=0; i<nlig_atoms; ++i)
        lig_types.push_back(type_dist(engine));

    unsigned* glig_types;
    CUDA_CHECK_GNINA(cudaMalloc(&glig_types, sizeof(unsigned)*lig_types.size()));
    CUDA_CHECK_GNINA(cudaMemcpy(glig_types, &lig_types[0], sizeof(unsigned)*lig_types.size(), cudaMemcpyHostToDevice));
    dinfo.types = glig_types;
}

fl check_bounds_deriv(const vec& a_coords, vec& adjusted_a_coords, vec& out_of_bounds_deriv, grid_dims& gd, float slope) 
{
	fl out_of_bounds_penalty = 0;
	adjusted_a_coords = a_coords;
	VINA_FOR_IN(j, gd)
	{
		if (gd[j].n > 0)
		{
			if (a_coords[j] < gd[j].begin)
			{
				adjusted_a_coords[j] = gd[j].begin;
				out_of_bounds_deriv[j] = -1;
				out_of_bounds_penalty += std::abs(
						a_coords[j] - gd[j].begin);
			}
			else if (a_coords[j] > gd[j].end)
			{
				adjusted_a_coords[j] = gd[j].end;
				out_of_bounds_deriv[j] = 1;
				out_of_bounds_penalty += std::abs(a_coords[j] - gd[j].end);
			}
		}
	}
	out_of_bounds_penalty *= slope;
	out_of_bounds_deriv *= slope;
	return out_of_bounds_penalty;
}

float eval_inter_cpu(const precalculate_splines* p, unsigned num_movable_atoms, const std::vector<atom_params>& lig_atoms, const std::vector<atom_params>& rec_atoms, const std::vector<unsigned>& lig_types, const std::vector<unsigned>& rec_types, std::vector<force_energy_tup>& minus_forces, grid_dims& gd, float slope, float v) {
    //derived from non_cache::eval_deriv but modified such that it does not
    //depend on the existence of a model object. might be desirable to put this
    //somewhere else, but currently it's only used here anyway
	fl e = 0;
	const fl cutoff_sqr = p->cutoff_sqr();

    //TODO: generalize
	sz n = 27;

	VINA_FOR(i, num_movable_atoms)
	{
		const auto& a = lig_atoms[i];
		auto t1 = lig_types[i];
		if (t1 >= n || is_hydrogen((smina_atom_type::type)t1))
		{
			minus_forces[i] = force_energy_tup(0,0,0,0);
			continue;
		}

		const vec& a_coords = vec(a.coords[0], a.coords[1], a.coords[2]);
		vec adjusted_a_coords;
		vec out_of_bounds_deriv(0, 0, 0);
		fl out_of_bounds_penalty = check_bounds_deriv(a_coords, adjusted_a_coords, out_of_bounds_deriv, gd, slope);
		
		fl this_e = 0;
		vec deriv(0, 0, 0);
		VINA_FOR_IN(j, rec_atoms)
		{
			const auto& b = rec_atoms[j];
			auto t2 = rec_types[j];
			vec r_ba;
			r_ba = adjusted_a_coords - vec(b.coords[0], b.coords[1], b.coords[2]);
			fl r2 = sqr(r_ba);
			if (r2 < cutoff_sqr)
			{
				if (r2 < epsilon_fl) {
					throw std::runtime_error(
							"Ligand atom exactly overlaps receptor atom.  I can't deal with this.");
				}
                //construct correct atom_base objects
                atom_base a_ab;
                a_ab.sm = (smina_atom_type::type)lig_types[i];
                a_ab.charge = lig_atoms[i].charge;
                atom_base b_ab;
                b_ab.sm = (smina_atom_type::type)rec_types[j];
                b_ab.charge = rec_atoms[j].charge;

				pr e_dor = p->eval_deriv(a_ab, b_ab, r2);
				this_e += e_dor.first;
				deriv += e_dor.second * r_ba;
			}
		}
		curl(this_e, deriv, v);
        vec tmp = deriv + out_of_bounds_deriv;
		minus_forces[i] = *(force_energy_tup*)(&tmp);
		e += this_e + out_of_bounds_penalty;
	}
	return e;
}

int main() {
    //TODO: include progress bar?
    //set up c++11 random number engine
    auto const seed = std::random_device()();
    std::mt19937 engine(seed);

    //set up logging 
    int ntests = 1;
    int* failed_tests[ntests];

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

    weighted_terms wt(&t, t.weights());

    //set up splines
    const precalculate_gpu* gprec = new precalculate_gpu(wt, approx_factor);
    const precalculate_splines* prec = new precalculate_splines(wt, approx_factor);

    //set up fake lig for testing
    std::vector<atom_params> lig_atoms;
    fl max_x = -HUGE_VALF, max_y = -HUGE_VALF, max_z = -HUGE_VALF;
    fl min_x = HUGE_VALF, min_y = HUGE_VALF, min_z = HUGE_VALF;
    make_lig(lig_atoms, engine);

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

    //set up dinfo, mostly involves randomly generating some receptor atoms.
    //to my knowledge lig_penalties is never used, so I don't even bother
    //mallocing it here; if I'm wrong a segfault is a good way to find out
    GPUNonCacheInfo dinfo;
    dinfo.gridends = float3(gd[0].begin, gd[1].begin, gd[2].begin);
    dinfo.gridbegins = float3(gd[0].end, gd[1].end, gd[2].end);
    dinfo.ntypes = gprec->num_types();
    dinfo.splineInfo = gprec->getDeviceData();

    std::vector<unsigned> lig_types;
    std::vector<atom_params> rec_atoms;
    std::vector<unsigned> rec_types;
    make_dinfo(dinfo, engine, lig_atoms.size(), lig_types, rec_atoms, rec_types);
    std::vector<force_energy_tup> minus_forces(lig_atoms.size());
    force_energy_tup* d_forces;
    cudaMalloc(&d_forces, sizeof(force_energy_tup)*minus_forces.size());
    cudaMemcpy(d_forces, &minus_forces[0], sizeof(force_energy_tup)*minus_forces.size(), cudaMemcpyHostToDevice);

    //malloc and copy lig_atoms
    atom_params* glig_atoms;
    cudaMalloc(&glig_atoms, sizeof(atom_params) * lig_atoms.size());
    cudaMemcpy(glig_atoms, &lig_atoms[0], sizeof(atom_params) * lig_atoms.size(), cudaMemcpyHostToDevice);

    //get intermolecular energy, check agreement
    float g_out = single_point_calc(dinfo, glig_atoms, d_forces, v[0]);
    // float c_out = eval_inter_cpu(prec, dinfo.num_movable_atoms, lig_atoms, rec_atoms, lig_types, rec_types, minus_forces, gd, dinfo.slope, v[1]);

    printf("%f\n", g_out);
    // printf("%f\n", c_out);

    cudaFree(glig_atoms);
    cudaFree(d_forces);
    cudaFree(dinfo.types);
    delete gprec;
    delete prec;
}
