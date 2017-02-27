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
#include "conf_gpu.h"
#include <numeric>

inline void minus_mat_vec_product(const flmat& m, const change& in, change& out)
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


inline void minus_mat_vec_product(const flmat_gpu& m, const change_gpu& in, change_gpu& out)
{
	in.minus_mat_vec_product(m, out);
}

inline fl scalar_product(const change& a, const change& b, sz n)
{
	fl tmp = 0;
	VINA_FOR(i, n)
		tmp += a(i) * b(i);
	return tmp;
}

inline fl scalar_product(const change_gpu& a, const change_gpu& b, sz n)
{
	return a.dot(b);
}


inline bool bfgs_update(flmat& h, const change& p, const change& y,
		const fl alpha)
{
	const fl yp = scalar_product(y, p, h.dim());
	if (alpha * yp < epsilon_fl)
		return false; // FIXME?
	change minus_hy(y);
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

void bfgs_update(const flmat_gpu& h, const change_gpu& p, const change_gpu& y,
		const fl alpha);

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

inline fl compute_lambdamin(const change& p, const conf& x, sz n)
{
	fl test = 0;
	//compute lambdamin
	for (sz i = 0; i < n; i++)
	{
		//static_assert(std::is_same<decltype(std::fabs(1.0f)),float>::value,"Not a float.\n");
		fl temp = std::fabs(p(i)) / std::max(std::fabs(x(i)), 1.0f);
		if (temp > test)
			test = temp;
	}
	return test;
}

fl compute_lambdamin(const change_gpu& p, const conf_gpu& x, sz n);

//dkoes - this line search is modeled after lnsrch in numerical recipes, it puts
//a bit of effort into calculating a good scaling factor, and ensures that alpha
//will actually result in a smaller value
template<typename F, typename Conf, typename Change>
fl accurate_line_search(F& f, sz n, const Conf& x, const Change& g, const fl f0,
		const Change& p, Conf& x_new, Change& g_new, fl& f1)
{ // returns alpha
	fl a, alpha, alpha2 = 0, alamin, b, disc, f2 = 0;
	fl rhs1, rhs2, slope = 0, sum = 0, test, tmplam;
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
	test = compute_lambdamin(p, x, n);

	alamin = std::numeric_limits<float>::epsilon() / test;
	alpha = FIRST; //single newton step
	for (;;) //always try full newton step first
	{
		x_new = x;
		x_new.increment(p, alpha);

		f1 = f(x_new, g_new);

		//std::cout << "alpha " << alpha << "  f " << f1 << "\tslope " << slope << " f0ALF " << f0 + ALF * alpha * slope << "\n";
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

void set_diagonal(const flmat_gpu& m, fl x);

inline void subtract_change(change& b, const change& a, sz n)
{ // b -= a
	VINA_FOR(i, n)
		b(i) -= a(i);
}

inline void subtract_change(change_gpu& b, const change_gpu& a, sz n)
{ // b -= a
	b.sub(a);
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
	std::cout << std::setprecision(8);
	std::cout << "f0 " << f0 << "\n";
	//std::ofstream fout("minout.sdf");
	VINA_U_FOR(step, params.maxiters)
	{
		minus_mat_vec_product(h, g, p);
		fl f1 = 0;
		fl alpha;
/*
		f.m->set(x);
		f.m->write_sdf(fout);
		fout << "$$$$\n";
*/
		if (params.type == minimization_params::BFGSAccurateLineSearch)
			alpha = accurate_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
		else
			alpha = fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);

		if(alpha == 0) {
			//std::cout << "alpha 0\n";
			break; //line direction was wrong, give up
		}

		Change y(g_new);
        // Update line direction
		subtract_change(y, g, n);

		fl prevf0 = f0;
		f0 = f1;
		x = x_new;

		if (params.early_term)
		{
			//dkoes - use the progress in reducing the function value as an indication of when to stop
			fl diff = prevf0 - f0;
			if (std::fabs(diff) < 1e-5) //arbitrary cutoff
			{
				break;
			}
		}

		g = g_new; // dkoes - check the convergence of the new gradient

		fl gradnormsq = scalar_product(g, g, n);
		std::cout << "step " << step << " " << f0 << " " << gradnormsq << " " << alpha << "\n";

		if (!(gradnormsq >= 1e-4)) //slightly arbitrary cutoff - works with fp
		{
			//std::cout << "gradnormsq " << gradnormsq << "\n";
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
//	std::cout << "final f0 " << f0 << "\n";

	return f0;
}

template<typename F>
fl bfgs(F& f, conf_gpu& x, change_gpu& g, const fl average_required_improvement,
		const minimization_params& params) {
    sz n = g.num_floats();

    // Initialize and copy Hessian
    flmat_gpu h(n);

    // Initialize and copy additional conf and change objects
	change_gpu g_new(g);
	conf_gpu x_new(x);
    //TODO: change model::eval_deriv_gpu so it doesn't memcpy at the end
	fl f0 = f(x, g);
	fl f_orig = f0;
	change_gpu g_orig(g);
	conf_gpu x_orig(x);

	change_gpu p(g);
    // For now, keep everything on the GPU but control is maintained on CPU
    // which launches the relevant kernels...maybe this will change
	VINA_U_FOR(step, params.maxiters)
	{
		minus_mat_vec_product(h, g, p);
        // f1 is the returned energy for the next iteration of eval_deriv_gpu
		fl f1 = 0;
		fl alpha;

		if (params.type == minimization_params::BFGSAccurateLineSearch)
			alpha = accurate_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
		else
			alpha = fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);

		if(alpha == 0) {
			break;
		}

		change_gpu y(g_new);
        // Update line direction
		subtract_change(y, g, n);

		fl prevf0 = f0;
		f0 = f1;
		x = x_new;

		if (params.early_term)
		{
			fl diff = prevf0 - f0;
			if (std::fabs(diff) < 1e-5) 
			{
				break;
			}
		}

		g = g_new; 

		fl gradnormsq = scalar_product(g, g, n);
//		std::cout << "step " << step << " " << f0 << " " << gradnormsq << " " << alpha << "\n";

		if (!(gradnormsq >= 1e-4)) //slightly arbitrary cutoff - works with fp
		{
			//std::cout << "gradnormsq " << gradnormsq << "\n";
			break; // breaks for nans too // FIXME !!??
		}

		if (step == 0)
		{
			const fl yy = scalar_product(y, y, n);
			if (std::abs(yy) > epsilon_fl)
				set_diagonal(h, alpha * scalar_product(y, p, n) / yy);
		}
        // bfgs_update used to return a bool, but the value of that bool never
        // got checked anyway
		bfgs_update(h, p, y, alpha);
	}

	if (!(f0 <= f_orig))
	{ // succeeds for nans too
		f0 = f_orig;
		x = x_orig;
		g = g_orig;
	}
//	std::cout << "final f0 " << f0 << "\n";
    CUDA_CHECK_GNINA(cudaFree(h.m_data));
	return f0;
}

#endif
