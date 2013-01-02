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

struct precalculate_element {
	precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) {}
	fl eval_fast(fl r2) const {
		assert(r2 * factor < fast.size());
		sz i = sz(factor * r2);  // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i < fast.size()); 
		return fast[i];
	}
	pr eval_deriv(fl r2) const {
		fl r2_factored = factor * r2;
		assert(r2_factored + 1 < smooth.size());
		sz i1 = sz(r2_factored); 
		sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i1 < smooth.size());
		assert(i2 < smooth.size());
		fl rem = r2_factored - i1;
		assert(rem >= -epsilon_fl);
		assert(rem < 1 + epsilon_fl);
		const pr& p1 = smooth[i1];
		const pr& p2 = smooth[i2];
		fl e   = p1.first  + rem * (p2.first  - p1.first);
		fl dor = p1.second + rem * (p2.second - p1.second); 
		return pr(e, dor);
	}
	void init_from_smooth_fst(const flv& rs) {
		sz n = smooth.size();
		VINA_CHECK(rs.size() >= n);
		VINA_CHECK(fast.size() == n);
		VINA_FOR(i, n) {
			// calculate dor's
			fl& dor = smooth[i].second;
			if(i == 0 || i == n-1)
				dor = 0;
			else {
				fl delta = rs[i+1] - rs[i-1];
				fl r = rs[i];
				dor = (smooth[i+1].first - smooth[i-1].first) / (delta * r);
			}
			// calculate fast's from smooth.first's
			fl f1 = smooth[i].first;
			fl f2 = (i+1 >= n) ? 0 : smooth[i+1].first;
			fast[i] = (f2 + f1) / 2;
		}
	}

	flv fast;
	prv smooth; // [(e, dor)]
	fl factor;
};

struct precalculate {
	precalculate(const scoring_function& sf, const minimization_params& minparms, fl v = max_fl, fl factor_ = 32) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
		m_cutoff(sf.cutoff()),
		m_cutoff_sqr(sqr(sf.cutoff())),
		n(sz(factor_ * m_cutoff_sqr) + 3),  // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
		cutoff_smoothing(minparms.cutoff_smoothing),
		factor(factor_),
		scoring(sf),
		data(num_atom_types(sf.atom_typing_used()), precalculate_element(n, factor_)),
		m_atom_typing_used(sf.atom_typing_used()) {

		VINA_CHECK(factor > epsilon_fl);
		VINA_CHECK(sz(m_cutoff_sqr*factor) + 1 < n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
		VINA_CHECK(m_cutoff_sqr*factor + 1 < n);

		calculate_rs();

		VINA_FOR(t1, data.dim())
			VINA_RANGE(t2, t1, data.dim()) {
				precalculate_element& p = data(t1, t2);
				// init smooth[].first
				VINA_FOR_IN(i, p.smooth)
					p.smooth[i].first = (std::min)(v, sf.eval_fast(t1, t2, rs[i]));

				//dkoes - need to smooth the transition into the cutoff or
				//artifacts appear in the energy landscape
				if (cutoff_smoothing > 0)
				{
					fl smoothrange = minparms.cutoff_smoothing;
					fl start = std::max(m_cutoff - smoothrange, 0.0);
					fl startsq = sqr(start);
					sz rstart = startsq * factor;
					smoothrange = rs[rstart]; //adjust to boundary
					sz end = factor * m_cutoff_sqr;
					//get the endpoint value at the cutoff
					if (rstart < end)
					{
						fl endval = p.smooth[end].first;
						//zero the very end
						for (unsigned i = end; i < p.smooth.size(); i++)
						{
							p.smooth[i].first = 0;
						}
						//subtract off increments of endval to end up at zero
						for (unsigned i = rstart; i < end; i++)
						{
								p.smooth[i].first -= endval
									* (rs[i] + smoothrange - m_cutoff)
									/ smoothrange;
						}
					}
				}
				// init the rest
				p.init_from_smooth_fst(rs);
			}
	}
	fl eval_fast(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_fast(r2);
	}

	fl eval_slow(const atom_base& a, const atom_base& b, fl r2) const {
		//dkoes - un-precalculable terms - widening isn't supported here
		if(scoring.has_slow())
		{ //dkoes - this check is just to avoid the sqrt..
			fl r = sqrt(r2);
			return scoring.eval_slow(a, b, r);
		}
		return 0;
	}

	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		sz type_pair_index = get_type_pair_index(atom_typing_used(), a, b);
		pr ret = data(type_pair_index).eval_deriv(r2);
		if(scoring.has_slow())
		{
			//dkoes - recompute "derivative" computation on the fly,
			//I am attempting to exactly mimic the precomputation, including
			//the discretization
			fl r2_factored = factor * r2;
			//compute rounded positions w,x,y,z, with r between x and y
			sz x = sz(r2_factored);
			if(x > 0)
			{
				sz w = x-1;
				sz y = x+1;
				sz z = x+2;
				//value at positions w,x,y,z
				fl W = scoring.eval_slow(a,b,rs[w]);
				fl X = scoring.eval_slow(a,b,rs[x]);
				fl Y = scoring.eval_slow(a,b,rs[y]);
				fl Z = scoring.eval_slow(a,b,rs[z]);

				if(cutoff_smoothing > 0)
				{
					sz rstart = factor*(m_cutoff_sqr-sqr(cutoff_smoothing));

					if(z > rstart) //have to smooth at least some values
					{
						fl smoothrange = m_cutoff-rs[rstart];
						fl endval = scoring.eval_slow(a,b,m_cutoff);
						Z -= endval* (rs[z] + smoothrange - m_cutoff)
								/ smoothrange;
						if(y > rstart)
							Y -= endval* (rs[y] + smoothrange - m_cutoff)/ smoothrange;
						if(x > rstart)
							X -= endval* (rs[x] + smoothrange - m_cutoff)/ smoothrange;
						if(w > rstart)
							W -= endval* (rs[w] + smoothrange - m_cutoff)/ smoothrange;
					}
				}

				fl rem = r2_factored - x; //how much beyond y we are

				fl e   = X  + rem * (Y  - X); //linearly interpolate

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
	sz index_permissive(sz t1, sz t2) const { return data.index_permissive(t1, t2); }
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }
	fl cutoff_sqr() const { return m_cutoff_sqr; }
	bool has_slow() const { return scoring.has_slow(); } //dkoes
private:
	void calculate_rs()  {
		rs = flv(n+2, 0); //dkoes - so I don't have to be careful with eval slow
		VINA_FOR(i, n+2)
			rs[i] = std::sqrt(i / factor);
	}
	fl m_cutoff;
	fl m_cutoff_sqr;
	sz n;
	fl factor;
	fl cutoff_smoothing;
	atom_type::t m_atom_typing_used;
	const scoring_function& scoring;
	triangular_matrix<precalculate_element> data;
	flv rs; //actual distance of index locations
};

#endif
