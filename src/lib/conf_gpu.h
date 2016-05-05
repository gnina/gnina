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

#ifndef VINA_CONF_GPU_H
#define VINA_CONF_GPU_H

#include <boost/ptr_container/ptr_vector.hpp> // typedef output_container

#include "quaternion.h"
#include "random.h"
#include "conf.h"
#include <iostream>
#include <string>

struct change_gpu {
	std::vector<ligand_change> ligands;
	std::vector<residue_change> flex;
	change_gpu(const change& c) {
		cudaMalloc(this, sizeof(change_gpu));
		cudaMemcpy(&this, &c, sizeof(change), cudaMemcpyHostToDevice);
	}	
	change_gpu(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0);
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0);
	}
	//dkoes - zeros out all differences
	void clear()
	{
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid = rigid_change();
			ligands[i].torsions.assign(ligands[i].torsions.size(),0);
		}
		VINA_FOR_IN(i, flex)
			flex[i].torsions.assign(flex[i].torsions.size(), 0);
	}

	//dkoes - multiply by -1
	void invert()
	{
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid.position *= -1;
			ligands[i].rigid.orientation *= -1;
			for(unsigned j = 0, n = ligands[i].torsions.size(); j < n; j++)
				ligands[i].torsions[j] *= -1;
		}
		VINA_FOR_IN(i, flex) {
			for(unsigned j = 0, n = flex[i].torsions.size(); j < n; j++)
				flex[i].torsions[j] *= -1;
		}
	}
	fl operator()(sz index) const { // returns by value
		VINA_FOR_IN(i, ligands) {
			const ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			const residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return 0; // shouldn't happen, placating the compiler
	}
	fl& operator()(sz index) {
		VINA_FOR_IN(i, ligands) {
			ligand_change& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3) return lig.rigid.orientation[index];
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			residue_change& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false); 
		return ligands[0].rigid.position[0]; // shouldn't happen, placating the compiler
	}
	sz num_floats() const {
		sz tmp = 0;
		VINA_FOR_IN(i, ligands)
			tmp += 6 + ligands[i].torsions.size();
		VINA_FOR_IN(i, flex)
			tmp += flex[i].torsions.size();
		return tmp;
	}
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
};

struct conf_gpu {
	std::vector<ligand_conf> ligands;
	std::vector<residue_conf> flex;
	conf_gpu(const conf& c) {
		for (int i=0; i<c.ligands.size(); i++) {
			ligands.push_back(c.ligands[i]);
		}
		for (int i=0; i<c.flex.size(); i++) {
			flex.push_back(c.flex[i]);
		}
	}
	conf_gpu() {}
	conf_gpu(const conf_size& s) : ligands(s.ligands.size()), flex(s.flex.size()) {
		VINA_FOR_IN(i, ligands)
			ligands[i].torsions.resize(s.ligands[i], 0); // FIXME?
		VINA_FOR_IN(i, flex)
			flex[i].torsions.resize(s.flex[i], 0); // FIXME?
	}
	conf to_conf() {
		conf c;
		c.ligands = ligands;
		c.flex = flex;
		return c;
	}
	void set_to_null() {
		VINA_FOR_IN(i, ligands)
			ligands[i].set_to_null();
		VINA_FOR_IN(i, flex)
			flex[i].set_to_null();
	}
	void increment(const change_gpu& c, fl factor) { // torsions get normalized, orientations do not
		VINA_FOR_IN(i, ligands)
			ligands[i].increment(c.ligands[i], factor);
		VINA_FOR_IN(i, flex)
			flex[i]   .increment(c.flex[i],    factor);
	}
	bool internal_too_close(const conf_gpu& c, fl torsions_cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!torsions_too_close(ligands[i].torsions, c.ligands[i].torsions, torsions_cutoff))
				return false;
		return true;
	}
	bool external_too_close(const conf_gpu& c, const scale& cutoff) const {
		assert(ligands.size() == c.ligands.size());
		VINA_FOR_IN(i, ligands)
			if(!ligands[i].rigid.too_close(c.ligands[i].rigid, cutoff.position, cutoff.orientation))
				return false;
		assert(flex.size() == c.flex.size());
		VINA_FOR_IN(i, flex)
			if(!torsions_too_close(flex[i].torsions, c.flex[i].torsions, cutoff.torsion))
				return false;
		return true;
	}
	bool too_close(const conf_gpu& c, const scale& cutoff) const {
		return internal_too_close(c, cutoff.torsion) &&
			   external_too_close(c, cutoff); // a more efficient implementation is possible, probably
	}
	void generate_internal(fl torsion_spread, fl rp, const conf* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			ligands[i].rigid.position.assign(0);
			ligands[i].rigid.orientation = qt_identity;
			const flv* torsions_rs = rs ? (&rs->ligands[i].torsions) : NULL;
			torsions_generate(ligands[i].torsions, torsion_spread, rp, torsions_rs, generator);
		}
	}
	void generate_external(const scale& spread, fl rp, const conf_gpu* rs, rng& generator) { // torsions are not normalized after this
		VINA_FOR_IN(i, ligands) {
			const rigid_conf* rigid_conf_rs = rs ? (&rs->ligands[i].rigid) : NULL;
			ligands[i].rigid.generate(spread.position, spread.orientation, rp, rigid_conf_rs, generator);
		}
		VINA_FOR_IN(i, flex) {
			const flv* torsions_rs = rs ? (&rs->flex[i].torsions) : NULL;
			torsions_generate(flex[i].torsions, spread.torsion, rp, torsions_rs, generator);
		}
	}
	void randomize(const vec& corner1, const vec& corner2, rng& generator) {
		VINA_FOR_IN(i, ligands)
			ligands[i].randomize(corner1, corner2, generator);
		VINA_FOR_IN(i, flex)
			flex[i].randomize(generator);
	}
	void print() const {
		VINA_FOR_IN(i, ligands)
			ligands[i].print();
		VINA_FOR_IN(i, flex)
			flex[i].print();
	}
	//dkoes - index into position values; corresponds to change indexing
	//read only because of quaternions
	fl operator()(sz index) const { // returns by value
		VINA_FOR_IN(i, ligands) {
			const ligand_conf& lig = ligands[i];
			if(index < 3) return lig.rigid.position[index];
			index -= 3;
			if(index < 3)
			{
				vec ang = quaternion_to_angle(lig.rigid.orientation);
				return ang[index];
			}
			index -= 3;
			if(index < lig.torsions.size()) return lig.torsions[index];
			index -= lig.torsions.size();
		}
		VINA_FOR_IN(i, flex) {
			const residue_conf& res = flex[i];
			if(index < res.torsions.size()) return res.torsions[index];
			index -= res.torsions.size();
		}
		VINA_CHECK(false);
		return 0; // shouldn't happen, placating the compiler
	}

private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & ligands;
		ar & flex;
	}
};

struct output_type_gpu {
	conf_gpu c;
	fl e;
	vecv coords;
	output_type_gpu(const conf_gpu& c_, fl e_) : c(c_), e(e_) {}
	output_type_gpu(output_type& ot) : c(ot.c), e(ot.e), coords(ot.coords) {
		cudaMalloc(&outgpu, sizeof(output_type_gpu));
		cudaMemcpy(outgpu, &out, sizeof(output_type), cudaMemcpyHostToDevice);
	}
};

typedef boost::ptr_vector<output_type_gpu> output_container_gpu;

inline bool operator<(const output_type_gpu& a, const output_type_gpu& b) { // for sorting output_container
	return a.e < b.e;
}

#endif
