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

/*
 * SMINA NOTICE
 * dkoes - if you add a term here, in order to use it with --custom_scoring
 * you must add a parser to custom_terms.{h,cpp}
 *
 */

#ifndef VINA_EVERYTHING_H
#define VINA_EVERYTHING_H

#include "terms.h"
#include "int_pow.h"

inline fl gaussian(fl x, fl width)
{
	return std::exp(-sqr(x / width));
}

// distance_additive terms

template<unsigned i>
struct electrostatic: public distance_additive
{
	fl cap;
	electrostatic(fl cap_, fl cutoff_) :
			distance_additive(cutoff_), cap(cap_)
	{
		name = std::string("electrostatic(i=") + to_string(i) + ",_^="
				+ to_string(cap) + ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom_base& a, const atom_base& b, fl r) const
	{
		fl tmp = int_pow<i>(r);
		fl q1q2 = a.charge * b.charge;
		if (tmp < epsilon_fl)
			return q1q2 * cap;
		else
			return q1q2 * (std::min)(cap, 1 / int_pow<i>(r));
	}
};

struct ad4_solvation: public distance_additive
{
	fl desolvation_sigma;
	fl solvation_q;
	bool charge_dependent;
	ad4_solvation(fl desolvation_sigma_, fl solvation_q_,
			bool charge_dependent_, fl cutoff_) :
			distance_additive(cutoff_), solvation_q(solvation_q_), charge_dependent(
					charge_dependent_), desolvation_sigma(desolvation_sigma_)
	{
		name = std::string("ad4_solvation(d-sigma=")
				+ to_string(desolvation_sigma) + ",_s/q="
				+ to_string(solvation_q) + ",_q=" + to_string(charge_dependent)
				+ ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(const atom_base& a, const atom_base& b, fl r) const;
};

inline fl optimal_distance(sz xs_t1, sz xs_t2)
{
	return xs_radius(xs_t1) + xs_radius(xs_t2);
}

struct gauss: public usable
{
	fl offset; // added to optimal distance
	fl width;
	gauss(fl offset_, fl width_, fl cutoff_) :
			usable(cutoff_), offset(offset_), width(width_)
	{
		name = std::string("gauss(o=") + to_string(offset) + ",_w="
				+ to_string(width) + ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		return gaussian(r - (optimal_distance(t1, t2) + offset), width);
	}
};

struct repulsion: public usable
{
	fl offset; // added to vdw
	repulsion(fl offset_, fl cutoff_) :
			usable(cutoff_), offset(offset_)
	{
		name = std::string("repulsion(o=") + to_string(offset) + ",_c="
				+ to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		fl d = r - (optimal_distance(t1, t2) + offset);
		if (d > 0)
			return 0;
		return d * d;
	}
};

inline fl slope_step(fl x_bad, fl x_good, fl x)
{
	if (x_bad < x_good)
	{
		if (x <= x_bad)
			return 0;
		if (x >= x_good)
			return 1;
	}
	else
	{
		if (x >= x_bad)
			return 0;
		if (x <= x_good)
			return 1;
	}
	return (x - x_bad) / (x_good - x_bad);
}

struct hydrophobic: public usable
{
	fl good;
	fl bad;
	hydrophobic(fl good_, fl bad_, fl cutoff_) :
			usable(cutoff_), good(good_), bad(bad_)
	{
		name = "hydrophobic(g=" + to_string(good) + ",_b=" + to_string(bad)
				+ ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		if (xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		else
			return 0;
	}
};

struct non_hydrophobic: public usable
{
	fl good;
	fl bad;
	non_hydrophobic(fl good_, fl bad_, fl cutoff_) :
			usable(cutoff_), good(good_), bad(bad_)
	{
		name = "non_hydrophobic(g=" + to_string(good) + ",_b=" + to_string(bad)
				+ ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		if (!xs_is_hydrophobic(t1) && !xs_is_hydrophobic(t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		else
			return 0;
	}
};

template<unsigned n, unsigned m>
void find_vdw_coefficients(fl position, fl depth, fl& c_n, fl& c_m)
{
	BOOST_STATIC_ASSERT(n != m);
	c_n = int_pow<n>(position) * depth * m / (fl(n) - fl(m));
	c_m = int_pow<m>(position) * depth * n / (fl(m) - fl(n));
}

template<unsigned i, unsigned j>
struct vdw: public usable
{
	fl smoothing;
	fl cap;
	vdw(fl smoothing_, fl cap_, fl cutoff_)
	:
			usable(cutoff_), smoothing(smoothing_), cap(cap_)
	{
		name = "vdw(i=" + to_string(i) + ",_j=" + to_string(j) + ",_s="
				+ to_string(smoothing) + ",_^=" + to_string(cap) + ",_c="
				+ to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		fl d0 = optimal_distance(t1, t2);
		fl depth = 1;
		fl c_i = 0;
		fl c_j = 0;
		find_vdw_coefficients<i, j>(d0, depth, c_i, c_j);
		if (r > d0 + smoothing)
			r -= smoothing;
		else if (r < d0 - smoothing)
			r += smoothing;
		else
			r = d0;

		fl r_i = int_pow<i>(r);
		fl r_j = int_pow<j>(r);
		if (r_i > epsilon_fl && r_j > epsilon_fl)
			return (std::min)(cap, c_i / r_i + c_j / r_j);
		else
			return cap;
	}
};

/* A 10-12 LJ potential */
struct non_dir_h_bond_lj: public usable
{
	fl offset;
	fl cap;
	non_dir_h_bond_lj(fl offset_, fl cap_, fl cutoff_) :
			usable(cutoff_), offset(offset_), cap(cap_)
	{
		name = std::string("non_dir_h_bond_lj(o=") + to_string(offset)
				+ ",_^=" + to_string(cap)
				+ ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		if (xs_h_bond_possible(t1, t2))
		{
			fl d0 = optimal_distance(t1, t2)+offset;
			fl depth = 5;
			fl c_i = 0;
			fl c_j = 0;
			find_vdw_coefficients<10, 12>(d0, depth, c_i, c_j);

			fl r_i = int_pow<10>(r);
			fl r_j = int_pow<12>(r);
			if (r_i > epsilon_fl && r_j > epsilon_fl)
				return (std::min)(cap, c_i / r_i + c_j / r_j);
			else
				return cap;
		}
		return 0;
	}
};

/* This mimics repulsion, but only between possible bonders. More for testing */
struct non_dir_h_bond_quadratic: public usable
{
	fl offset;
	non_dir_h_bond_quadratic(fl offset_, fl cutoff_) :
			usable(cutoff_), offset(offset_)
	{
		name = std::string("non_dir_h_bond_quadratic(o=") + to_string(offset)
				+ ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		if (xs_h_bond_possible(t1, t2))
		{
			fl d = r - (optimal_distance(t1, t2) + offset);
			if (d > 0)
				return 0;
			return d * d;
		}
		return 0;
	}
};

//classic Vina hbond term
struct non_dir_h_bond: public usable
{
	fl good;
	fl bad;
	non_dir_h_bond(fl good_, fl bad_, fl cutoff_) :
			usable(cutoff_), good(good_), bad(bad_)
	{
		name = std::string("non_dir_h_bond(g=") + to_string(good) + ",_b="
				+ to_string(bad) + ",_c=" + to_string(cutoff) + ")";
	}
	fl eval(sz t1, sz t2, fl r) const
	{
		if (xs_h_bond_possible(t1, t2))
			return slope_step(bad, good, r - optimal_distance(t1, t2));
		return 0;
	}
};

inline fl read_iterator(flv::const_iterator& i)
{
	fl x = *i;
	++i;
	return x;
}

fl smooth_div(fl x, fl y);

struct num_tors_add: public conf_independent
{
	num_tors_add()
	{
		name = "num_tors_add";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		//fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		fl w = read_iterator(i); // FIXME?
		return x + w * in.num_tors;
	}
};

struct num_tors_sqr: public conf_independent
{
	num_tors_sqr()
	{
		name = "num_tors_sqr";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		fl add = w * sqr(fl(in.num_tors)) / 5;
		fl ret = x + add;
		return ret;
	}
};

struct num_tors_sqrt: public conf_independent
{
	num_tors_sqrt()
	{
		name = "num_tors_sqrt";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 0.1 * read_iterator(i); // [-1 .. 1]
		return x + w * std::sqrt(fl(in.num_tors)) / sqrt(5.0);
	}
};

struct num_tors_div: public conf_independent
{
	num_tors_div()
	{
		name = "num_tors_div";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 0.1 * (read_iterator(i) + 1); // w is in [0..0.2]
		return smooth_div(x, 1 + w * in.num_tors / 5.0);
	}
};

struct ligand_length: public conf_independent
{
	ligand_length()
	{
		name = "ligand_length";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = read_iterator(i);
		return x + w * in.ligand_lengths_sum;
	}
};

struct num_ligands: public conf_independent
{
	num_ligands()
	{
		name = "num_ligands";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 1 * read_iterator(i); // w is in [-1.. 1]
		return x + w * in.num_ligands;
	}
};

struct num_heavy_atoms_div: public conf_independent
{
	num_heavy_atoms_div()
	{
		name = "num_heavy_atoms_div";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 0.05 * read_iterator(i);
		return smooth_div(x, 1 + w * in.num_heavy_atoms);
	}
};

struct num_heavy_atoms: public conf_independent
{
	num_heavy_atoms()
	{
		name = "num_heavy_atoms";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 0.05 * read_iterator(i);
		return x + w * in.num_heavy_atoms;
	}
};

struct num_hydrophobic_atoms: public conf_independent
{
	num_hydrophobic_atoms()
	{
		name = "num_hydrophobic_atoms";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = 0.05 * read_iterator(i);
		return x + w * in.num_hydrophobic_atoms;
	}
};

struct constant_term: public conf_independent
{
	constant_term()
	{
		name = "constant_term";
	}
	sz size() const
	{
		return 1;
	}
	fl eval(const conf_independent_inputs& in, fl x,
			flv::const_iterator& i) const
	{
		fl w = read_iterator(i);
		return x + w;
	}
};

#endif
