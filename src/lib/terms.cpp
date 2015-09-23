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

#include "terms.h"
#include "brick.h"

conf_independent_inputs::operator flv() const
{
	flv tmp;
	tmp.push_back(num_tors);
	tmp.push_back(num_rotors);
	tmp.push_back(num_heavy_atoms);
	tmp.push_back(num_hydrophobic_atoms);
	tmp.push_back(ligand_max_num_h_bonds);
	tmp.push_back(num_ligands);
	tmp.push_back(ligand_lengths_sum);
	return tmp;
}

unsigned conf_independent_inputs::num_bonded_heavy_atoms(const model& m,
		const atom_index& i) const
		{ // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned acc = 0;
	const std::vector<bond>& bonds = m.get_atom(i).bonds;
	VINA_FOR_IN(j, bonds)
	{
		const bond& b = bonds[j];
		const atom& a = m.get_atom(b.connected_atom_index);
		if (!a.is_hydrogen())
			++acc;
	}
	return acc;
}

// FIXME? - could be static, but I don't feel like declaring function friends
unsigned conf_independent_inputs::atom_rotors(const model& m,
		const atom_index& i) const
		{ // the number of rotatable bonds to heavy ligand atoms
	unsigned acc = 0;
	const atom& base = m.get_atom(i);
	const std::vector<bond>& bonds = base.bonds;

	VINA_FOR_IN(j, bonds)
	{
		const bond& b = bonds[j];
		const atom& a = m.get_atom(b.connected_atom_index);
		if (b.rotatable && !a.is_hydrogen()
				&& num_bonded_heavy_atoms(m, b.connected_atom_index) > 1)
		{

			if (num_bonded_heavy_atoms(m, i) > 1
					|| !get_fixed_rotable_hydrogens())
			{
				// not counting CH_3, etc
				// dkoes - only check one way if hbonds aren't fixed to
				// be bug-compatible with Vina
				++acc;
			}
		}
	}
	return acc;
}

conf_independent_inputs::conf_independent_inputs(const model& m)
{
	num_tors = 0;
	num_rotors = 0;
	num_heavy_atoms = 0;
	num_hydrophobic_atoms = 0;
	ligand_max_num_h_bonds = 0;
	num_ligands = m.ligands.size();
	ligand_lengths_sum = 0;

	VINA_FOR_IN(ligand_i, m.ligands)
	{
		const ligand& lig = m.ligands[ligand_i];
		ligand_lengths_sum += m.ligand_length(ligand_i);
		VINA_RANGE(i, lig.begin, lig.end)
		{
			const atom& a = m.atoms[i];
			if (!a.is_hydrogen())
			{
				unsigned ar = atom_rotors(m, atom_index(i, false));

				num_tors += 0.5 * ar;

				if (ar > 2)
					num_rotors += 0.5;
				else
					num_rotors += 0.5 * ar;

				++num_heavy_atoms;
				if (xs_is_hydrophobic(a.get()))
					++num_hydrophobic_atoms;

				if (xs_is_acceptor(a.get()) || xs_is_donor(a.get()))
					++ligand_max_num_h_bonds;
			}
		}
	}
}

std::vector<std::string> conf_independent_inputs::get_names() const
{ // FIXME should probably be static
	std::vector<std::string> tmp;
	tmp.push_back("num_tors");
	tmp.push_back("num_rotors");
	tmp.push_back("num_heavy_atoms");
	tmp.push_back("num_hydrophobic_atoms");
	tmp.push_back("ligand_max_num_h_bonds");
	tmp.push_back("num_ligands");
	tmp.push_back("ligand_lengths_sum");
	VINA_CHECK(static_cast<flv>(*this).size() == tmp.size());
	// FIXME?
	return tmp;
}

conf_independent_inputs::conf_independent_inputs() :
		num_tors(0), num_rotors(0), num_heavy_atoms(0),
				num_hydrophobic_atoms(0), ligand_max_num_h_bonds(0), num_ligands(
						0),
				ligand_lengths_sum(0)
{
}

inline fl inner_product_shortest(const flv& a, const flv& b)
{
	sz n = (std::min)(a.size(), b.size());
	fl acc = 0;
	VINA_FOR(i, n)
		acc += a[i] * b[i];
	return acc;
}

fl factors::eval(const flv& weights, bool include_internal) const
		{
	fl tmp = inner_product_shortest(e, weights);
	if (include_internal)
		tmp += inner_product_shortest(i, weights);
	return tmp;
}

// terms

std::vector<std::string> terms::get_names(bool enabled_only) const
		{ // does not include conf-independent
	std::vector<std::string> tmp;

	charge_independent_terms.get_names(enabled_only, tmp);
	charge_dependent_terms.get_names(enabled_only, tmp);
	distance_additive_terms.get_names(enabled_only, tmp);
	additive_terms.get_names(enabled_only, tmp);
	intermolecular_terms.get_names(enabled_only, tmp);
	return tmp;
}

sz terms::size_internal() const
{
	return distance_additive_terms.size() + charge_independent_terms.size()
			+
			charge_dependent_terms.size() + additive_terms.size();
}

sz terms::size_conf_independent(bool enabled_only) const
		{ // number of parameters does not necessarily equal the number of operators
	sz acc = 0;
	VINA_FOR_IN(i, conf_independent_terms)
		if (!enabled_only || conf_independent_terms.enabled[i])
			acc += conf_independent_terms[i].size();
	return acc;
}

fl terms::max_r_cutoff() const
{
	fl tmp = 0;
	tmp = (std::max)(tmp, charge_independent_terms.max_cutoff());
	tmp = (std::max)(tmp, charge_dependent_terms.max_cutoff());
	tmp = (std::max)(tmp, distance_additive_terms.max_cutoff());
	tmp = (std::max)(tmp, additive_terms.max_cutoff());
	return tmp;
}

void terms::eval_additive_aux(const model& m, const atom_index& i,
		const atom_index& j, fl r, flv& out) const
		{ // out is added to
	const atom& a = m.get_atom(i);
	const atom& b = m.get_atom(j);

	sz offset = 0;

	VINA_FOR_IN(k, charge_independent_terms)
		if (r < charge_independent_terms[k].cutoff)
			out[offset + k] += charge_independent_terms[k].eval(a, b, r);
	offset += charge_independent_terms.size();

	VINA_FOR_IN(k, charge_dependent_terms)
		if (r < charge_dependent_terms[k].cutoff)
			out[offset + k] += charge_dependent_terms[k].eval(a, b, r);
	offset += charge_dependent_terms.size();

	VINA_FOR_IN(k, distance_additive_terms)
		if (r < distance_additive_terms[k].cutoff)
			out[offset + k] += distance_additive_terms[k].eval(a, b, r);
	offset += distance_additive_terms.size();

	VINA_FOR_IN(k, additive_terms)
		if (r < additive_terms[k].cutoff)
			out[offset + k] += additive_terms[k].eval(m, i, j);

	VINA_CHECK(offset + additive_terms.size() == size_internal());
}

//dkoes - evaluate from terms directly with no precomputation
//also, fill out each atom's contribution
//this routine is not intended to be efficient (see precalc)
flv terms::evale_robust(const model& m, std::vector<flv>& per_atom) const
{
	flv tmp(size(), 0);
	// only single-ligand systems are supported by this procedure
	if(m.ligands.size() == 0) //nolig
	  return tmp;

  VINA_CHECK(m.ligands.size() == 1);
	fl max_r_cutoff_sqr = sqr(max_r_cutoff());

	grid_dims box = m.movable_atoms_box(0); // add nothing
	vec box_begin = grid_dims_begin(box);
	vec box_end = grid_dims_end(box);

	const sz n = num_atom_types();

	std::vector<atom_index> relevant_atoms;

	VINA_FOR_IN(j, m.grid_atoms)
	{
		const atom& a = m.grid_atoms[j];
		const smt t = a.get();
		if (brick_distance_sqr(box_begin, box_end, a.coords)
				< max_r_cutoff_sqr && t < n && !is_hydrogen(t)) // exclude, say, Hydrogens
			relevant_atoms.push_back(atom_index(j, true));
	}

	per_atom.clear();
	per_atom.resize(m.atoms.size(), flv(size(), 0));
	VINA_FOR_IN(j, m.atoms)
	{
		const atom& a = m.atoms[j];
		const vec& a_coords = m.coords[j];
		if (m.find_ligand(j) < m.ligands.size())
			continue; // skip ligand atoms, add only flex/inflex
		const smt t = a.get();
		if (brick_distance_sqr(box_begin, box_end, a_coords)
				< max_r_cutoff_sqr && t < n && !is_hydrogen(t)) // exclude, say, Hydrogens
			relevant_atoms.push_back(atom_index(j, false));
	}

	VINA_FOR_IN(lig_i, m.ligands)
	{
		const ligand& lig = m.ligands[lig_i];
		VINA_RANGE(i, lig.begin, lig.end)
		{
			const vec& coords = m.coords[i];
			const atom& a = m.atoms[i];
			const smt t = a.get();

			if (t < n && !is_hydrogen(t))
			{ // exclude, say, Hydrogens
				VINA_FOR_IN(relevant_j, relevant_atoms)
				{
					const atom_index& j = relevant_atoms[relevant_j];
					fl d2 = vec_distance_sqr(coords, m.atom_coords(j));
					if (d2 > max_r_cutoff_sqr)
						continue; // most likely scenario
					fl d = std::sqrt(d2);
					eval_additive_aux(m, atom_index(i, false), j, d, per_atom[i]);
				}

				for(sz j = 0, m = tmp.size(); j < m; j++)
					tmp[j] += per_atom[i][j];
			}
		}
	}

	sz offset = size_internal();
	VINA_CHECK(intermolecular_terms.size() == 0);
	VINA_CHECK(offset + intermolecular_terms.size() == tmp.size());

	return tmp;
}


fl terms::eval_conf_independent(const conf_independent_inputs& in, fl x,
		flv::const_iterator& it) const
		{ // evaluates enabled only
	VINA_FOR_IN(i, conf_independent_terms)
		if (conf_independent_terms.enabled[i])
			x = conf_independent_terms[i].eval(in, x, it);
	return x;
}

flv terms::filter_external(const flv& v) const
		{
	flv tmp;
	flv::const_iterator i = v.begin();
	charge_independent_terms.filter(i, tmp);
	charge_dependent_terms.filter(i, tmp);
	distance_additive_terms.filter(i, tmp);
	additive_terms.filter(i, tmp);
	intermolecular_terms.filter(i, tmp);
	VINA_CHECK(i == v.end());
	return tmp;
}

flv terms::filter_internal(const flv& v) const
		{
	flv tmp;
	flv::const_iterator i = v.begin();
	charge_independent_terms.filter(i, tmp);
	charge_dependent_terms.filter(i, tmp);
	distance_additive_terms.filter(i, tmp);
	additive_terms.filter(i, tmp);
	VINA_CHECK(i == v.end());
	return tmp;
}

factors terms::filter(const factors& f) const
		{
	factors tmp;
	tmp.e = filter_external(f.e);
	tmp.i = filter_internal(f.i);
	return tmp;
}

void terms::display_info() const
{
	std::vector<std::string> enabled_names = get_names(true);
	std::cout << "Enabled terms: \n";
	VINA_FOR_IN(i, enabled_names)
		std::cout << enabled_names[i] << '\n';
	std::cout << '\n';

	std::vector<std::string> enabled_operators;
	conf_independent_terms.get_names(true, enabled_operators);
	std::cout << "Enabled conf-independent operators: \n";
	VINA_FOR_IN(i, enabled_operators)
		std::cout << enabled_operators[i] << '\n';
	std::cout << '\n';

	VINA_SHOW(size());
	VINA_SHOW(size_internal());
	VINA_SHOW(max_r_cutoff());
}
