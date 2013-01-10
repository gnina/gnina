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
		temp = fabs(p(i)) / std::max(fabs(x(i)), 1.0);
		if (temp > test)
			test = temp;
	}

	alamin = std::numeric_limits<fl>::epsilon() / test;
	alpha = FIRST; //single newton step
	for (;;) //always try full newton step first
	{
		x_new = x;
		x_new.increment(p, alpha);
		f1 = f(x_new, g_new);
//std::cout << "alpha " << alpha << "  f " << f1 << "\n";
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
		alpha = std::max(tmplam, 0.1 * alpha); //never smaller than a tenth
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

	VINA_U_FOR(step, params.maxiters)
	{
		minus_mat_vec_product(h, g, p);
		fl f1 = 0;
		fl alpha;

		if (params.type == minimization_params::BFGSAccurateLineSearch)
			alpha = accurate_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
		else
			alpha = fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);

		Change y(g_new);
		subtract_change(y, g, n);

		fl prevf0 = f0;
		f0 = f1;
		x = x_new;

		if (params.early_term)
		{
			//dkoes - since the gradient check doesn't actually work, use the
			//progress in reducing the function value as an indication of when to stop
			fl diff = prevf0 - f0;
			if (fabs(diff) < 1e-5) //arbitrary cutoff
			{
				break;
			}
		}

		g = g_new; // dkoes - check the convergence of the new gradient

		//dkoes - in practice the gradient never seems to converge to zero
		//(it is set to zero when accurate linesearch fails though)
		fl gradnormsq = scalar_product(g, g, n);
//std::cout.precision(10);
//std::cout << "\n" << step << " f " << f0 << " grad " << sqrt(gradnormsq) << " alpha " << alpha << "\n";
//g.print();
//x.print();
		if (!(gradnormsq >= 1e-5 * 1e-5))
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

#endif
