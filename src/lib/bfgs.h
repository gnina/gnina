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

#ifndef VINA_BFGS_H
#define VINA_BFGS_H

#include "matrix.h"
#include <numeric>
typedef triangular_matrix<fl> flmat;

template<typename Change>
void minus_mat_vec_product(const flmat& m, const Change& in, Change& out)
{
	sz n = m.dim();
	VINA_FOR(i, n)
	{
		fl sum = 0;
		VINA_FOR(j, n)
			sum += m(m.index_permissive(i, j)) * in(j);
		out(i) = -sum;
	}
}

template<typename Change>
inline fl scalar_product(const Change& a, const Change& b, sz n)
{
	fl tmp = 0;
	VINA_FOR(i, n)
		tmp += a(i) * b(i);
	return tmp;
}

template<typename Change>
inline bool bfgs_update(flmat& h, const Change& p, const Change& y,
		const fl alpha)
{
	const fl yp = scalar_product(y, p, h.dim());
	if (alpha * yp < epsilon_fl)
		return false; // FIXME?
	Change minus_hy(y);
	minus_mat_vec_product(h, y, minus_hy);
	const fl yhy = -scalar_product(y, minus_hy, h.dim());
	const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
	const sz n = p.num_floats();
	VINA_FOR(i, n)
		VINA_RANGE(j, i, n) // includes i
			h(i, j) += alpha * r * (minus_hy(i) * p(j)
					+ minus_hy(j) * p(i)) +
					+alpha * alpha * (r * r * yhy + r) * p(i) * p(j); // s * s == alpha * alpha * p * p
	return true;
}

//dkoes - this is the line search method used by vina,
//it is simple and fast, but may return an inappropriately large alpha
template<typename F, typename Conf, typename Change>
fl fast_line_search(F& f, sz n, const Conf& x, const Change& g, const fl f0,
		const Change& p, Conf& x_new, Change& g_new, fl& f1)
{ // returns alpha
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	const fl pg = scalar_product(p, g, n);

	VINA_U_FOR(trial, max_trials)
	{
		x_new = x;
		x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
		if (f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
			break;
		alpha *= multiplier;
	}
	return alpha;
}

//dkoes - this line search is modeled after lnsrch in numerical recipes, it puts
//a bit of effort into calculating a good scaling factor, and ensures that alpha
//will actually result in a smaller value
template<typename F, typename Conf, typename Change>
fl accurate_line_search(F& f, sz n, const Conf& x, const Change& g, const fl f0,
		const Change& p, Conf& x_new, Change& g_new, fl& f1)
{ // returns alpha
	fl a, alpha, alpha2 = 0, alamin, b, disc, f2 = 0;
	fl rhs1, rhs2, slope = 0, sum = 0, temp, test, tmplam;
	int i;
	const fl ALF = 1.0e-4;
	const fl FIRST = 1.0;
	sum = scalar_product(p, p, n);
	sum = sqrt(sum);

	slope = scalar_product(g, p, n);
	if (slope >= 0)
	{
		//gradient isn't actually in a decreasing direction
		x_new = x;
		g_new.clear(); //dkoes - set gradient to zero
		return 0;
	}
	test = 0;
	//compue lambdamin
	for (i = 0; i < n; i++)
	{
		temp = fabs(p(i)) / std::max(fabs(x(i)), 1.0f);
		if (temp > test)
			test = temp;
	}

	alamin = std::numeric_limits<float>::epsilon() / test;
	alpha = FIRST; //single newton step
	for (;;) //always try full newton step first
	{
		x_new = x;
		x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
		//std::cout << "alpha " << alpha << "  f " << f1 << "\tslope " << slope << "\n";
		if (alpha < alamin) //convergence
		{
			x_new = x;
			g_new.clear(); //dkoes - set gradient to zero
			return 0;
		}
		else if (f1 <= f0 + ALF * alpha * slope)
		{
			//sufficient function decrease, stop searching
			return alpha;
		}
		else //have to backtrack
		{
			if (alpha == FIRST)
			{
				//first time
				tmplam = -slope / (2.0 * (f1 - f0 - slope));
			}
			else //subsequent backtracks
			{
				rhs1 = f1 - f0 - alpha * slope;
				rhs2 = f2 - f0 - alpha2 * slope;
				a = (rhs1 / (alpha * alpha) - rhs2 / (alpha2 * alpha2))
						/ (alpha - alpha2);
				b = (-alpha2 * rhs1 / (alpha * alpha)
						+ alpha * rhs2 / (alpha2 * alpha2)) / (alpha - alpha2);
				if (a == 0.0)
					tmplam = -slope / (2.0 * b);
				else
				{
					disc = b * b - 3.0 * a * slope;
					if (disc < 0)
						tmplam = 0.5 * alpha;
					else if (b <= 0)
						tmplam = (-b + sqrt(disc)) / (3.0 * a);
					else
						tmplam = -slope / (b + sqrt(disc));
				}
				if (tmplam > .5 * alpha)
					tmplam = .5 * alpha; //always at least cut in half
			}
		}
		alpha2 = alpha;
		f2 = f1;
		//std::cout << "TMPLAM " << tmplam << "\n";
		//considered slowing things down with f1 > 0, but it was slow without actually improving scores
		alpha = std::max(tmplam, (fl)0.1 * alpha); //never smaller than a tenth
	}
}

inline void set_diagonal(flmat& m, fl x)
{
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

template<typename Change>
void subtract_change(Change& b, const Change& a, sz n)
{ // b -= a
	VINA_FOR(i, n)
		b(i) -= a(i);
}

template<typename F, typename Conf, typename Change>
fl bfgs(F& f, Conf& x, Change& g, const fl average_required_improvement,
		const minimization_params& params)
{ // x is I/O, final value is returned
	sz n = g.num_floats();
	flmat h(n, 0);
	set_diagonal(h, 1);

	Change g_new(g);
	Conf x_new(x);
	fl f0 = f(x, g);

	fl f_orig = f0;
	Change g_orig(g);
	Conf x_orig(x);

	Change p(g);

//	std::ofstream fout("minout.sdf");
	VINA_U_FOR(step, params.maxiters)
	{
		minus_mat_vec_product(h, g, p);
		fl f1 = 0;
		fl alpha;

//		f.m->set(x);
//		f.m->write_sdf(fout);
//		fout << "$$$$\n";

		if (params.type == minimization_params::BFGSAccurateLineSearch)
			alpha = accurate_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
		else
			alpha = fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);

		if(alpha == 0)
			break; //line direction was wrong, give up

		Change y(g_new);
		subtract_change(y, g, n);

		fl prevf0 = f0;
		f0 = f1;
		x = x_new;

		if (params.early_term)
		{
			//dkoes - use the progress in reducing the function value as an indication of when to stop
			fl diff = prevf0 - f0;
			if (fabs(diff) < 1e-5) //arbitrary cutoff
			{
				break;
			}
		}

		g = g_new; // dkoes - check the convergence of the new gradient

		fl gradnormsq = scalar_product(g, g, n);
		//std::cout << "step " << step << " " << f0 << " " << gradnormsq << " " << alpha << "\n";

		if (!(gradnormsq >= 1e-4)) //slightly arbitrary cutoff - works with fp
		{
			break; // breaks for nans too // FIXME !!??
		}

		if (step == 0)
		{
			const fl yy = scalar_product(y, y, n);
			if (std::abs(yy) > epsilon_fl)
				set_diagonal(h, alpha * scalar_product(y, p, n) / yy);
		}

		bool h_updated = bfgs_update(h, p, y, alpha);
	}

	if (!(f0 <= f_orig))
	{ // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
	return f0;
}

//set g = g_new + B*g
template<typename Change>
void conjugate_update(Change& s, fl B, const Change& g_new, sz n)
{
	VINA_FOR(i, n)
	{
		s(i) = -g_new(i) + B*s(i);
	}
}

//dkoes - conjugate gradient method, from
//http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
//This does not do nearly as well either in terms of convergence or performance.
//This may be partly due to an inadequate line search method or a buggy
//implementation, but I don't feel compelled to invest any more time into it.
template<typename F, typename Conf, typename Change>
fl conjgrad(F& f, Conf& x, Change& g, const fl average_required_improvement,
		const minimization_params& params)
{ // x is I/O, final value is returned
	sz n = g.num_floats();

	Change g_new(g);
	Conf x_new(x);
	fl f0 = f(x, g);

	fl f_orig = f0;
	Change g_orig(g);
	Conf x_orig(x);
	Change s(g);
	s.invert();

	VINA_U_FOR(step, params.maxiters)
	{
		fl f1 = 0;
		fl alpha;

		//update position with line search and calculate new gradient
		//WARNING: NR says this method isn't accurate enought
		alpha = accurate_line_search(f, n, x, g, f0, s, x_new, g_new, f1);

		//compute B, multiplier of previous gradient

		fl denom = scalar_product(g, g, n);

		if(!(denom >= 1e-5*1e-5))
		{
			break; //very small gradient, consider ourselves converged
		}

		if(alpha == 0) {
			//line direction wrong, giveup
			break;
		}

		//use Polak-Ribiere:
		/*
		Change y(g_new);
		subtract_change(y, g, n);
		fl numerator = scalar_product(g_new, y, n);
		fl B = numerator/denom;
		*/
		//use Fletcher-REeves
		fl numerator = scalar_product(g_new, g_new, n);
		fl B = numerator/denom;
		B = std::max((fl)0.0, B);

		//update direction
		//s = -g_new + s*B
		conjugate_update(s, B, g_new, n);

		fl prevf0 = f0;
		f0 = f1;
		x = x_new;
		g = g_new;

		if (params.early_term)
		{
			//dkoes - use the
			//progress in reducing the function value as an indication of when to stop
			fl diff = prevf0 - f0;
			if (fabs(diff) < 1e-5) //arbitrary cutoff
			{
				break;
			}
		}

		//check convergence of new gradient
		fl gradnormsq = scalar_product(g, g, n);
		if (!(gradnormsq >= 1e-5 * 1e-5))
		{
			break; // breaks for nans too // FIXME !!??
		}

	}

	if (!(f0 <= f_orig))
	{ // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
	return f0;
}

#endif
