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

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
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

	precalculate(const scoring_function& sf) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			m_cutoff(sf.cutoff()),
					m_cutoff_sqr(sqr(sf.cutoff())),
					scoring(sf)
	{

	}

	virtual ~precalculate()
	{

	}

	fl cutoff_sqr() const
	{
		return m_cutoff_sqr;
	}
	bool has_components() const
	{
		return scoring.num_used_components() > 1;
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

	fl eval(const atom_base& a, const atom_base& b, fl r2) const
	{
		fl ret = eval_fast(a.get(), b.get(), r2).eval(a, b);
		return ret + eval_slow(a, b, r2);
	}
protected:
	fl m_cutoff;
	fl m_cutoff_sqr;
	const scoring_function& scoring;

};

typedef std::vector<prv> prvv; //index by component, then point
class precalculate_linear_element
{
	friend class precalculate_linear;
	precalculate_linear_element(sz n, sz num_components, fl factor_) :
			fast(n, flv(num_components, 0)),
					smooth(num_components, prv(n, pr(0, 0))), //index by control point then by component
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

	pr eval_deriv(sz num_components, const atom_base& a, const atom_base& b,
			fl r2) const
			{
		fl r2_factored = factor * r2;
		assert(smooth.size() == num_components);
		assert(r2_factored + 1 < smooth[0].size());
		sz i1 = sz(r2_factored);
		sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i1 < smooth[0].size());
		assert(i2 < smooth[0].size());
		fl rem = r2_factored - i1;
		assert(rem >= -epsilon_fl);
		assert(rem < 1 + epsilon_fl);
		fl e1, e2, d1, d2;
		if (num_components == 1) //very slight speedup here
		{
			e1 = smooth[0][i1].first;
			e2 = smooth[0][i2].first;
			d1 = smooth[0][i1].second;
			d2 = smooth[0][i2].second;
		}
		else
		{
			result_components e1comp, e2comp, d1comp, d2comp;
			for (sz c = 0; c < num_components; c++)
			{
				e1comp[c] = smooth[c][i1].first;
				e2comp[c] = smooth[c][i2].first;
				d1comp[c] = smooth[c][i1].second;
				d2comp[c] = smooth[c][i2].second;
			}
			e1 = e1comp.eval(a, b);
			e2 = e2comp.eval(a, b);
			d1 = d1comp.eval(a, b);
			d2 = d2comp.eval(a, b);
		}

		fl e = e1 + rem * (e2 - e1);
		fl dor = d1 + rem * (d2 - d1);
		return pr(e, dor);
	}

	void init_from_smooth_fst(sz num_components, const flv& rs)
	{
		sz n = smooth[0].size();
		VINA_CHECK(rs.size() >= n);
		VINA_CHECK(fast.size() == n);
		for (sz c = 0; c < num_components; c++)
		{
			VINA_FOR(i, n)
			{
				// calculate dor's
				fl& dor = smooth[c][i].second;
				if (i == 0 || i == n - 1)
					dor = 0;
				else
				{
					fl delta = rs[i + 1] - rs[i - 1];
					fl r = rs[i];
					dor = (smooth[c][i + 1].first - smooth[c][i - 1].first)
							/ (delta * r);
				}
				// calculate fast's from smooth.first's
				fl f1 = smooth[c][i].first;
				fl f2 = (i + 1 >= n) ? 0 : smooth[c][i + 1].first;
				fast[i][c] = (f2 + f1) / 2;
			}
		}

	}

	std::vector<result_components> fast;
	prvv smooth; // [(e, dor)] for each component, indexed first by component
	fl factor;
};

class precalculate_linear: public precalculate
{
	//evaluate data while properly swapping types
	result_components eval_fast_data(smt t1, smt t2, fl r2) const
	{
		if (t1 <= t2)
		{
			return data(t1, t2).eval_fast(r2);
		}
		else
		{
			result_components ret = data(t2, t1).eval_fast(r2);
			ret.swapOrder();
			return ret;
		}
	}

public:
	precalculate_linear(const scoring_function& sf, fl factor_) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf),
					n(sz(factor_ * m_cutoff_sqr) + 3), // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
					data(num_atom_types(),
							precalculate_linear_element(n,
									sf.num_used_components(), factor_)),
					num_components(sf.num_used_components()),
					factor(factor_)
	{
		VINA_CHECK(factor > epsilon_fl);
		VINA_CHECK(sz(m_cutoff_sqr*factor) + 1 < n);
		// cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
		VINA_CHECK(m_cutoff_sqr*factor + 1 < n);

		calculate_rs();

		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim())
			{
				precalculate_linear_element& p = data(t1, t2);
				// init smooth[].first
				VINA_FOR(i, n)
				{
					result_components res = sf.eval_fast((smt) t1, (smt) t2,
							rs[i]);
					for (sz c = 0; c < num_components; c++)
						p.smooth[c][i].first = res[c];
				}
				// init the rest
				p.init_from_smooth_fst(num_components, rs);
			}
	}

	result_components eval_fast(smt t1, smt t2, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		return eval_fast_data(t1, t2, r2);
	}

	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		smt t1 = a.get();
		smt t2 = b.get();
		pr ret;
		if (t1 <= t2)
			ret = data(t1, t2).eval_deriv(num_components, a, b, r2);
		else
			ret = data(t2, t1).eval_deriv(num_components, b, a, r2);

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
	sz num_components;
	fl factor;

	void calculate_rs() //calculate square roots of control points once
	{
		rs = flv(n + 2, 0); //dkoes - so I don't have to be careful with eval slow
		VINA_FOR(i, n+2)
			rs[i] = std::sqrt(i / factor);
	}
};

typedef std::pair<result_components, result_components> component_pair;
//evaluates spline between two smina atom types as needed
//will decompose charge dependent terms
class spline_cache
{
	const scoring_function* sf;
	fl cutoff;
	sz n;
	smt t1, t2;
	//the following are only computed when needed
	mutable Spline *splines; //one for each component
	mutable unsigned num_splines; //size of splines
	mutable bool valid;
	mutable boost::mutex lock; //make thread safe

	//create control points for spline
	//poitns indexed by component first; nonzero indexec by component
	void setup_points(std::vector<std::vector<pr> >& points,
			std::vector<bool>& nonzero) const
			{
		assert(n >= 2);
		fl fraction = cutoff / (fl) n;
		sz numc = sf->num_used_components();

		//clear out arguments
		nonzero.assign(numc, false);
		points.resize(numc);
		for (sz i = 0; i < numc; i++)
		{
			points[i].clear();
			points[i].reserve(numc + 1);
		}

		//compute points
		for (unsigned i = 0; i < n; i++)
		{
			fl xval = i * fraction;
			result_components res = sf->eval_fast(t1, t2, xval);
			for (unsigned c = 0; c < numc; c++)
			{
				points[c].push_back(pr(xval, res[c]));
				if (res[c] != 0)
					nonzero[c] = true;
			}
		}
		//last point at cutoff is zero
		for (unsigned c = 0; c < numc; c++)
		{
			points[c].push_back(pr(cutoff, 0));
		}
	}
public:

	spline_cache() :
			sf(NULL), cutoff(0), n(0), t1(smina_atom_type::NumTypes), t2(
					smina_atom_type::NumTypes), splines(NULL), num_splines(0), valid(false)
	{
	}

	//need explicit copy constructor to deal with mutex vairable
	spline_cache(const spline_cache& rhs) :
			sf(rhs.sf), cutoff(rhs.cutoff), n(rhs.n), t1(rhs.t1), t2(rhs.t2), splines(NULL), num_splines(0), valid(false)
	{ //mutable member do not get copied
	}

	~spline_cache() {
		if(splines) delete [] splines;
	}

	//intialize values to approprate types etc - do not compute spline
	void set(const scoring_function& sf_, smt t1_, smt t2_, fl cut, sz n_)
	{
		sf = &sf_;
		cutoff = cut;
		n = n_;
		t1 = t1_;
		t2 = t2_;

		if (t1 > t2)
			std::swap(t1, t2);
	}

	component_pair eval(fl r) const
			{
		if (splines == NULL)
		{
			//create spline, thread safe
			boost::lock_guard<boost::mutex> L(lock);

			if (splines == NULL) //another thread didn't fix it for us
			{
				std::vector<std::vector<pr> > points;
				std::vector<bool> nonzero;
				setup_points(points, nonzero);

				num_splines = sf->num_used_components();
				Spline *tmpsplines = new Spline[num_splines];
				for (sz i = 0, n = num_splines; i < n; i++)
				{
					if (nonzero[i]) //worth interpolating
						tmpsplines[i].initialize(points[i]);
				}
				splines = tmpsplines; //this is assumed atomic
			}
		}

		result_components val, deriv;
		for (sz i = 0, n = num_splines; i < n; i++)
		{
			pr ret = splines[i].eval_deriv(r);
			val[i] = ret.first;
			deriv[i] = ret.second;
		}
		return component_pair(val, deriv);
	}
};

// dkoes - using cubic spline interpolation instead of linear for nice
// smooth gradients
class precalculate_splines: public precalculate
{
	//evaluates splines at t1/t2 and r, properly swaping result
	component_pair evaldata(smt t1, smt t2, fl r) const
			{
		if (t1 <= t2)
		{
			return data(t1, t2).eval(r);
		}
		else
		{
			component_pair ret = data(t2, t1).eval(r);
			ret.first.swapOrder();
			ret.second.swapOrder();
			return ret;
		}
	}
public:
	precalculate_splines(const scoring_function& sf, fl factor_) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf),
					data(num_atom_types(), spline_cache()),
					delta(0.000005),
					factor(factor_)
	{
		VINA_CHECK(factor > epsilon_fl);
		unsigned n = factor * m_cutoff;
		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim())
			{
				//initialize spline cache - this doesn't create the splines
				data(t1, t2).set(sf, (smt) t1, (smt) t2, m_cutoff, n);
			}
	}

	result_components eval_fast(smt t1, smt t2, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		fl r = sqrt(r2);
		return evaldata(t1, t2, r).first;
	}

	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		smt t1 = a.get();
		smt t2 = b.get();
		fl r = sqrt(r2);

		component_pair rets = evaldata(t1, t2, r);

		pr ret(rets.first.eval(a, b), rets.second.eval(a, b));

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
	fl factor;
};

// dkoes - do a full function recomputation (no precalculation)
class precalculate_exact: public precalculate
{
public:
	precalculate_exact(const scoring_function& sf) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf), delta(0.000005)
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
};

#endif
