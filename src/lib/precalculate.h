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

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include "scoring_function.h"
#include "matrix.h"
#include "splines.h"

//base class for precaluting classes
class precalculate
{
public:
	//return just the fast evaluation of types, no derivative
	virtual result_components eval_fast(smt t1, smt t2, fl r2) const = 0;

	//return value and derivative
	//IMPORTANT: derivative is scaled by sqrt(r2) so that when
	//multiplied by the direction vector the result is normalized
	virtual pr eval_deriv(const atom_base& a, const atom_base& b,
			fl r2) const = 0;

	precalculate(const scoring_function& sf,
			const minimization_params& minparms, fl factor_) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			m_cutoff(sf.cutoff()),
					m_cutoff_sqr(sqr(sf.cutoff())),
					cutoff_smoothing(minparms.cutoff_smoothing),
					factor(factor_),
					scoring(sf)
	{

		VINA_CHECK(factor > epsilon_fl);
	}

	virtual ~precalculate()
	{

	}

	fl cutoff_sqr() const
	{
		return m_cutoff_sqr;
	}
	bool has_slow() const
	{
		return scoring.has_slow();
	} //dkoes

	fl eval_slow(const atom_base& a, const atom_base& b, fl r2) const
			{
		//dkoes - un-precalculable terms - widening isn't supported here
		if (scoring.has_slow())
		{ //dkoes - this check is just to avoid the sqrt..
			fl r = sqrt(r2);
			return scoring.eval_slow(a, b, r);
		}
		return 0;
	}
protected:
	fl m_cutoff;
	fl m_cutoff_sqr;
	fl factor;
	fl cutoff_smoothing;
	const scoring_function& scoring;

};

typedef std::vector< std::pair<result_components, result_components> > compvec;
class precalculate_linear_element
{
	friend class precalculate_linear;
	precalculate_linear_element(sz n, fl factor_) :
			fast(n, flv(result_components::size(), 0)),
			smooth(n), //index by control point then by component
			factor(factor_)
	{
	}

	result_components eval_fast(fl r2) const
	{
		assert(r2 * factor < fast.size());
		sz i = sz(factor * r2); // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i < fast.size());
		return fast[i];
	}

	pr eval_deriv(const atom_base& a, const atom_base& b,fl r2) const
	{
		fl r2_factored = factor * r2;
		assert(r2_factored + 1 < smooth.size());
		sz i1 = sz(r2_factored);
		sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i1 < smooth.size());
		assert(i2 < smooth.size());
		fl rem = r2_factored - i1;
		assert(rem >= -epsilon_fl);
		assert(rem < 1 + epsilon_fl);
		fl e1 = smooth[i1].first.eval(a,b);
		fl e2 = smooth[i2].first.eval(a,b);
		fl d1 = smooth[i1].second.eval(a,b);
		fl d2 = smooth[i2].second.eval(a,b);

		fl e = e1 + rem * (e2 - e1);
		fl dor = d1 + rem * (d2 - d1);
		return pr(e, dor);
	}

	void init_from_smooth_fst(const flv& rs)
	{
		sz n = smooth.size();
		VINA_CHECK(rs.size() >= n);
		VINA_CHECK(fast.size() == n);
		VINA_FOR(i, n)
		{
			for(sz j = 0, m = result_components::size(); j < m; j++)
			{
			// calculate dor's
				fl& dor = smooth[i].second[j];
				if (i == 0 || i == n - 1)
					dor = 0;
				else
				{
					fl delta = rs[i + 1] - rs[i - 1];
					fl r = rs[i];
					dor = (smooth[i + 1].first[j] - smooth[i - 1].first[j]) / (delta * r);
				}
				// calculate fast's from smooth.first's
				fl f1 = smooth[i].first[j];
				fl f2 = (i + 1 >= n) ? 0 : smooth[i + 1].first[j];
				fast[i][j] = (f2 + f1) / 2;
			}
		}
	}

	std::vector<result_components> fast;
	compvec smooth; // [(e, dor)] for each component
	fl factor;
};

class precalculate_linear: public precalculate
{
public:
	precalculate_linear(const scoring_function& sf,
			const minimization_params& minparms, fl factor_) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf, minparms, factor_),
					n(sz(factor_ * m_cutoff_sqr) + 3), // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
					data(num_atom_types(),
							precalculate_linear_element(n, factor_))
	{
		VINA_CHECK(sz(m_cutoff_sqr*factor) + 1 < n);
		// cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
		VINA_CHECK(m_cutoff_sqr*factor + 1 < n);

		calculate_rs();

		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim())
			{
				precalculate_linear_element& p = data(t1, t2);
				// init smooth[].first
				VINA_FOR_IN(i, p.smooth)
				{
					p.smooth[i].first = sf.eval_fast((smt)t1, (smt)t2, rs[i]);
				}
				// init the rest
				p.init_from_smooth_fst(rs);
			}
	}
	result_components eval_fast(smt t1, smt t2, fl r2) const
			{
		if (t1 > t2)
			std::swap(t1, t2);
		assert(r2 <= m_cutoff_sqr);
		return data(t1, t2).eval_fast(r2);
	}

	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		sz type_pair_index = get_type_pair_index(a, b);
		pr ret = data(type_pair_index).eval_deriv(a,b,r2);
		if (scoring.has_slow())
		{
			//dkoes - recompute "derivative" computation on the fly,
			//I am attempting to exactly mimic the precomputation, including
			//the discretization
			fl r2_factored = factor * r2;
			//compute rounded positions w,x,y,z, with r between x and y
			sz x = sz(r2_factored);
			if (x > 0)
			{
				sz w = x - 1;
				sz y = x + 1;
				sz z = x + 2;
				//value at positions w,x,y,z
				fl W = scoring.eval_slow(a, b, rs[w]);
				fl X = scoring.eval_slow(a, b, rs[x]);
				fl Y = scoring.eval_slow(a, b, rs[y]);
				fl Z = scoring.eval_slow(a, b, rs[z]);

				fl rem = r2_factored - x; //how much beyond y we are

				fl e = X + rem * (Y - X); //linearly interpolate

				//calc derivitives
				fl delta0 = rs[y] - rs[w];
				fl dor0 = (Y - W) / (delta0 * rs[x]);

				fl delta1 = rs[z] - rs[x];
				fl dor1 = (Z - X) / (delta1 * rs[y]);

				fl dor = dor0 + rem * (dor1 - dor0);

				ret.first += e;
				ret.second += dor;
			}
		}
		return ret;
	}

private:
	sz n;
	triangular_matrix<precalculate_linear_element> data;
	flv rs; //actual distance of index locations

	void calculate_rs() //calculate square roots of control points once
	{
		rs = flv(n + 2, 0); //dkoes - so I don't have to be careful with eval slow
		VINA_FOR(i, n+2)
			rs[i] = std::sqrt(i / factor);
	}
};

//evaluates spline between two smina atom types as needed
//will decompose charge dependent terms
class spline_cache
{
	const scoring_function* sf;
	fl cutoff;
	sz n;
	sz numcut;
	smt t1, t2;
	//the following are only computed when needed
	mutable Spline<spline_cache> spline;
	mutable bool valid;
	public:

	spline_cache() :
			sf(NULL), cutoff(0), n(0), numcut(0), t1(smina_atom_type::NumTypes), t2(smina_atom_type::NumTypes), valid(false)
	{
	}

	//intialize values to approprate types etc - do not compute spline
	void set(const scoring_function& sf_, smt t1_, smt t2_, fl cut, sz n_,
			sz ncut)
	{
		sf = &sf_;
		cutoff = cut;
		n = n_;
		numcut = ncut;
		t1 = t1_;
		t2 = t2_;
		valid = false;

		if (t1 > t2)
			std::swap(t1, t2);
	}

	//function called by spline
	fl operator()(fl r) const
	{
		return sf->eval_fast(t1, t2, r);
	}

	pr eval(fl r) const
	{
		if (!valid)
		{
			//create spline
			spline.initialize(*this, cutoff, n, numcut);
			valid = true;
		}

		return spline(r);
	}
};

// dkoes - using cubic spline interpolation instead of linear for nice
// smooth gradients
class precalculate_splines: public precalculate
{
public:
	precalculate_splines(const scoring_function& sf,
			const minimization_params& minparms, fl factor_, fl v = max_fl) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf, minparms, factor_),
					data(num_atom_types(), spline_cache()),
					delta(0.000005)
	{
		unsigned n = factor * m_cutoff;
		unsigned numcut = cutoff_smoothing * factor;
		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim())
			{
				//initialize spline cache - this doesn't create the splines
				data(t1, t2).set(sf,(smt)t1,(smt)t2, m_cutoff, n, numcut);
			}
	}

	fl eval_fast(smt t1, smt t2, fl r2) const
	{
		assert(r2 <= m_cutoff_sqr);
		fl r = sqrt(r2);
		if (t1 > t2)
			std::swap(t1, t2);
		return data(t1, t2).eval(r).first;
	}

	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const
	{
		assert(r2 <= m_cutoff_sqr);
		smt t1 = a.get();
		smt t2 = b.get();
		if (t1 > t2)
			std::swap(t1, t2);
		fl r = sqrt(r2);
		pr ret = data(t1, t2).eval(r);

		if (scoring.has_slow())
		{
			//compute value and numerical derivative directly from function

			fl X = scoring.eval_slow(a, b, r);
			ret.first += X;

			fl rhi = r + delta;
			fl rlo = r - delta;
			if (rlo < 0)
				rlo = 0;
			if (rhi > m_cutoff)
				rhi = m_cutoff;

			fl W = scoring.eval_slow(a, b, rlo);
			fl Y = 0;
			if (rhi < m_cutoff)
				Y = scoring.eval_slow(a, b, rhi);

			fl dx = (Y - W) / (rhi - rlo);
			ret.second += dx;
		}
		ret.second /= r;
		return ret;
	}

private:

	triangular_matrix<spline_cache> data;
	fl delta;
};

// dkoes - do a full function recomputation (no precalculation)
class precalculate_exact: public precalculate
{
public:
	precalculate_exact(const scoring_function& sf,
			const minimization_params& minparms, fl factor_, fl v = max_fl) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf, minparms, factor_),
					delta(0.000005)
	{
	}
	result_components eval_fast(smt t1, smt t2, fl r2) const
	{
		assert(r2 <= m_cutoff_sqr);
		fl r = sqrt(r2);
		return scoring.eval_fast(t1, t2, r);
	}

	//numerical exact derivative - ignore cutoff for now
	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const
			{
		smt ta = a.get();
		smt tb = b.get();
		fl r = sqrt(r2);
		result_components res = scoring.eval_fast(ta, tb, r);

		fl X = res.eval(a, b);
		fl rhi = r + delta;
		fl rlo = r - delta;
		if (rlo < 0)
			rlo = 0;

		fl W = scoring.eval_fast(ta, tb, rlo).eval(a, b);
		fl Y = scoring.eval_fast(ta, tb, rhi).eval(a, b);
		if (scoring.has_slow())
		{
			X += scoring.eval_slow(a, b, r);

			W += scoring.eval_slow(a, b, rlo);
			Y += scoring.eval_slow(a, b, rhi);
		}

		fl dx = (Y - W) / (rhi - rlo);
		return pr(X, dx / r);
	}

private:

	fl delta;
	sz num_types;
};

#endif
