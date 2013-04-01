/*

 Copyright (c) 2006-2010, The Scripps Research Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 Author: Dr. Oleg Trott <ot14@columbia.edu>,
 The Olson Lab,
 The Scripps Research Institute

 */

#ifndef VINA_NON_CACHE_GPU_H
#define VINA_NON_CACHE_GPU_H

#ifdef SMINA_GPU

#include "non_cache.h"
#include "precalculate_gpu.h"

struct non_cache_gpu: public non_cache
{
private:
	GPUNonCacheInfo info; //host struct of device pointers
	GPUNonCacheInfo *dinfo;//device pointer to info;

public:
	non_cache_gpu(szv_grid_cache& gcache, const grid_dims& gd_,
			const precalculate_gpu* p_, fl slope_) :
	non_cache(gcache, gd_, p_, slope_)
	{
		const model& m = gcache.getModel();
		info.cutoff_sq = p->cutoff_sqr();

		unsigned natoms = m.num_movable_atoms();
		info.natoms = natoms;
		//allocate memory for positions, partial charges, and atom types of movable atoms
		cudaMalloc(&info.coords, sizeof(float) * natoms * 3);
		cudaMalloc(&info.charges, sizeof(float) * natoms);
		cudaMalloc(&info.types, sizeof(unsigned) * natoms);
		cudaMalloc(&info.energies, sizeof(float) * natoms);

		//allocate memory for minus forces of movable atoms, init to zero
		cudaMalloc(&info.minus_forces, sizeof(unsigned) * natoms * 3);
		cudaMemset(info.minus_forces, 0, sizeof(unsigned) * natoms * 3);
		//initialize atom types and partial charges
		std::vector<unsigned> htypes(natoms);
		std::vector<float> hcharges(natoms);

		VINA_FOR(i, natoms)
		{
			htypes[i] = m.atoms[i].get();
			hcharges[i] = m.atoms[i].charge;
		}
		cudaMemcpy(info.charges, &hcharges[0], sizeof(float) * natoms,
				cudaMemcpyHostToDevice);
		cudaMemcpy(info.types, &htypes[0], sizeof(unsigned) * natoms,
				cudaMemcpyHostToDevice);

		//allocate memory and initialize grid dimensions
		cudaMalloc(&info.gridends, sizeof(float) * 3);
		cudaMalloc(&info.gridbegins, sizeof(float) * 3);

		float hbegins[3] =
		{	gd[0].begin, gd[1].begin, gd[2].begin};
		float hends[3] =
		{	gd[0].end, gd[1].end, gd[2].end};

		assert(gd[0].n > 0);
		assert(gd[1].n > 0);
		assert(gd[2].n > 0);

		cudaMemcpy(info.gridbegins, hbegins, sizeof(float) * 3,
				cudaMemcpyHostToDevice);
		cudaMemcpy(info.gridends, hends, sizeof(float) * 3, cudaMemcpyHostToDevice);

		//figure out all possibly relevant receptor atoms
		szv recatomids;
		gcache.compute_relevant(gd_, recatomids);
		unsigned nrecatoms = recatomids.size();
		info.nrecatoms = nrecatoms;

		//allocate memory for positions, atom types, and partial charges of all
		//possibly relevant receptor atoms
		cudaMalloc(&info.recoords, nrecatoms * 3 * sizeof(float));
		cudaMalloc(&info.reccharges, nrecatoms * sizeof(float));
		cudaMalloc(&info.rectypes, nrecatoms * sizeof(unsigned));

		//initialize
		std::vector<float> hrecoords(3 * nrecatoms);
		std::vector<unsigned> hrectypes(nrecatoms);
		std::vector<float> hreccharges(nrecatoms);
		for (unsigned i = 0; i < nrecatoms; i++)
		{
			unsigned index = recatomids[i];
			const vec& c = m.grid_atoms[index].coords;
			for (unsigned j = 0; j < 3; j++)
			{
				hrecoords[3 * i + j] = c[j];
			}
			hreccharges[i] = m.grid_atoms[index].charge;
			hrectypes[i] = m.grid_atoms[index].get();
		}
		cudaMemcpy(info.recoords, &hrecoords[0], nrecatoms * sizeof(float) * 3,
				cudaMemcpyHostToDevice);
		cudaMemcpy(info.reccharges, &hreccharges[0], nrecatoms * sizeof(float),
				cudaMemcpyHostToDevice);
		cudaMemcpy(info.rectypes, &hrectypes[0], nrecatoms * sizeof(float),
				cudaMemcpyHostToDevice);

		info.ntypes = p_->num_types();
		info.splineInfo = p_->getDeviceData();
		//create struct
		cudaMalloc(&dinfo, sizeof(GPUNonCacheInfo));
		cudaMemcpy(dinfo, &info, sizeof(GPUNonCacheInfo), cudaMemcpyHostToDevice);
	}

	virtual ~non_cache_gpu()
	{
		//deallocate device memory
		cudaFree(info.coords);
		cudaFree(info.charges);
		cudaFree(info.types);
		cudaFree(info.minus_forces);
		cudaFree(info.energies);

		cudaFree(info.gridbegins);
		cudaFree(info.gridends);

		cudaFree(info.recoords);
		cudaFree(info.reccharges);
		cudaFree(info.rectypes);

		cudaFree(dinfo);
	}

	fl eval(const model& m, fl v) const
	{
		abort(); //not implemented
	}

	//evaluate the model on the gpu, v is the curl amount
	//sets m.minus_forces and returns total energy
	virtual fl eval_deriv(model& m, fl v) const
	{
		//clear energies
		unsigned natoms = m.num_movable_atoms();
		cudaMemset(info.energies, 0, sizeof(float) * natoms);

		float hcoords[natoms * 3];

		//update coordinates
		for (unsigned i = 0; i < natoms; i++)
		{
			const vec& c = m.coords[i];
			for (unsigned j = 0; j < 3; j++)
			{
				hcoords[3 * i + j] = c[j];
			}
		}
		cudaMemcpy(info.coords, hcoords, sizeof(float) * natoms * 3,
				cudaMemcpyHostToDevice);
		cudaMemset(info.minus_forces, 0, natoms*3*sizeof(float));
		cudaMemset(info.energies, 0, natoms*sizeof(float));

		//this will calculate the per-atom energies and forces; curl ignored
		double e = single_point_calc(dinfo, info.energies, slope, info.natoms, info.nrecatoms, v);

		//get forces
		float forces[natoms*3];
		cudaMemcpy(forces, info.minus_forces, natoms*3*sizeof(float), cudaMemcpyDeviceToHost);

		for(unsigned i = 0; i < natoms; i++) {
			for(unsigned j = 0; j < 3; j++) {
				m.minus_forces[i][j] = forces[3*i+j];
			}
		}

		return e;
	}

};

#endif
#endif

