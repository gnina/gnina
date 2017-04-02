#include "quasi_newton.h"
#include "conf_gpu.h"
#include "matrix.h"
#include "bfgs.h"

#include <cuda_runtime.h>

thread_local float_buffer buffer;

__device__ fl compute_lambdamin(const change_gpu& p, const conf_gpu& x, sz n)

{
    fl test = 0;
	for (sz i = 0; i < n; i++)
	{
		fl temp = fabsf(p.values[i]) / fmaxf(fabsf(x.values[i]),
                                                    1.0f);
		if (temp > test)
			test = temp;
	}
    return test;
}

//TODO: operator -=
__device__ inline 
void subtract_change(change_gpu& b, const change_gpu& a, sz n)
{ // b -= a
	b.sub(a);
}

__device__
void set_diagonal(flmat_gpu& m, fl x)
{
    VINA_FOR(i, m.dim())
		m(i, i) = x;
}

__device__ inline fl scalar_product(const change_gpu& a, const change_gpu& b, sz n)
{
	return a.dot(b);
}

__device__ inline void minus_mat_vec_product(const flmat_gpu& m,
                                             const change_gpu& in, change_gpu& out)
{
	in.minus_mat_vec_product(m, out);
}

__device__
fl accurate_line_search_gpu(quasi_newton_aux_gpu& f, sz n, const conf_gpu& x,
                            const change_gpu& g, const fl f0,
                            const change_gpu& p, conf_gpu& x_new,
                            change_gpu& g_new, fl& f1)
{ 
	fl a, alpha2 = 0, b, disc, f2 = 0;
	fl rhs1, rhs2, slope = 0, test, tmplam;
	const fl ALF = 1.0e-4;
	const fl FIRST = 1.0;

    __shared__ fl alpha;
    __shared__ fl alamin;
    int idx = threadIdx.x;
	slope = scalar_product(g, p, n);
	if (slope >= 0)
	{
		//gradient isn't actually in a decreasing direction
        if (idx == 0) {
		    x_new = x;
		    g_new.clear(); //dkoes - set gradient to zero
        }
		return 0;
	}
    if (idx == 0) {
	    test = compute_lambdamin(p, x, n);

	    alamin = epsilon_fl / test;
	    alpha = FIRST; //single newton step
    }
	for (;;) //always try full newton step first
	{
        if (idx == 0) {
		    x_new = x;
        }
        __syncthreads();
        if (idx < x_new.n)
		    x_new.increment(p, alpha, f.gdata.treegpu);

        if (idx == 0)
            f1 = f(x_new, g_new);

        __syncthreads();
		//std::cout << "alpha " << alpha << "  f " << f1 << "\tslope " << slope << " f0ALF " << f0 + ALF * alpha * slope << "\n";
		if (alpha < alamin) //convergence
		{
            if (idx == 0) {
			    x_new = x;
			    g_new.clear(); //dkoes - set gradient to zero
            }
			return 0;
		}
		else if (f1 <= f0 + ALF * alpha * slope)
		{
			//sufficient function decrease, stop searching
            return alpha;
		}
		else //have to backtrack
		{
            if (idx == 0) {
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
		}
        if (idx == 0) {
		    alpha2 = alpha;
		    f2 = f1;
		    //std::cout << "TMPLAM " << tmplam << "\n";
		    alpha = fmaxf(tmplam, (fl)0.1 * alpha); //never smaller than a tenth
        }
	}

	return 0; // absolutely necessary to workaround nvcc compiler bug!!! (only took N days to find..)

}

__device__
void bfgs_update(flmat_gpu& h, const change_gpu& p,
                 const change_gpu& y, const fl alpha,
                 change_gpu &minus_hy) {
	const fl yp = y.dot(p);
	const sz n = p.num_floats();
    int idx = threadIdx.x;
	if (alpha * yp < epsilon_fl)
		return; // FIXME?

    if (idx == 0)
        minus_hy = y;

    if (idx < minus_hy.n)
	    y.minus_mat_vec_product(h, minus_hy);

    __syncthreads();

	const fl yhy = -y.dot(minus_hy);
    if (idx < n) {
	    const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon

	    float coef = +alpha * alpha * (r * r * yhy + r) ;

        float *minus_hyvec = minus_hy.values;
        float *pvec = p.values;
	    VINA_RANGE(j, idx, n) // includes i
            atomicAdd(&h(idx, j), alpha * r *
                       (minus_hyvec[idx] * pvec[j] + minus_hyvec[j] * pvec[idx])
                       + coef * pvec[idx]	* pvec[j]);
    }
    __syncthreads();
    // s * s == alpha * alpha * p * p	} *
}

__global__
void bfgs_gpu(quasi_newton_aux_gpu f,
              conf_gpu x, conf_gpu x_orig, conf_gpu x_new,
              change_gpu g, change_gpu g_orig, change_gpu g_new,
              change_gpu p, change_gpu y, flmat_gpu h, change_gpu minus_hy,
              const fl average_required_improvement,
              const minimization_params params,
              float* out_energy)
{
    sz n = g.n;
    __shared__ fl alpha;
    __shared__ fl diff;
    __shared__ fl f1;
    __shared__ fl f0;
    float f_orig;
    int idx = threadIdx.x;

    if (idx == 0) {
        f0 = f(x, g);
	    f_orig = f0;
    }
	VINA_U_FOR(step, params.maxiters)
	{
        if (idx < g.n) {
		    minus_mat_vec_product(h, g, p);
            // f1 is the returned energy for the next iteration of eval_deriv_gpu
		    f1 = 0;
            //do we even care about the fast_line_search?
		    assert(params.type == minimization_params::BFGSAccurateLineSearch);
        }
        __syncthreads();
		alpha = accurate_line_search_gpu(f, n, x, g, f0,
                                                p, x_new, g_new, f1);
		if(alpha == 0) 
			break;
        fl prevf0;

        if (idx == 0) {
		    y = g_new;

		    prevf0 = f0;
		    f0 = f1;
		    x = x_new;
        }

        // Update line direction
        if (idx < y.n)
		    subtract_change(y, g, n);

		if (params.early_term)
		{
            if (idx == 0)
			    diff = prevf0 - f0;
            __syncthreads();
			if (fabsf(diff) < 1e-5) 
				break;
		}

        if (idx == 0) 
		    g = g_new; 

        __syncthreads();
		fl gradnormsq = scalar_product(g, g, n);
//		std::cout << "step " << step << " " << f0 << " " << gradnormsq << " " << alpha << "\n";

		if (!(gradnormsq >= 1e-4)) //slightly arbitrary cutoff - works with fp
			break; // breaks for nans too // FIXME !!??

		if (step == 0)
		{
			const fl yy = scalar_product(y, y, n);
			if (fabsf(yy) > epsilon_fl) {
                const fl yp = scalar_product(y, p, n);
                if (idx == 0)
				    set_diagonal(h, alpha * yp / yy);
            }
		}
        // bfgs_update used to return a bool, but the value of that bool never
        // got checked anyway
		bfgs_update(h, p, y, alpha, minus_hy);
	}
    if (idx == 0) {
	    if (!(f0 <= f_orig))
	    { // succeeds for nans too
	    	f0 = f_orig;
	    	x = x_orig;
	    	g = g_orig;
	    }
        *out_energy = f0;
    }
}

fl bfgs(quasi_newton_aux_gpu &f, conf_gpu& x,
        change_gpu& g, const fl average_required_improvement,
		const minimization_params& params) {
    sz n = g.num_floats();

    // Initialize and copy Hessian
    flmat_gpu h(n);

    // Initialize and copy additional conf and change objects
	change_gpu g_orig(g, buffer);
	change_gpu g_new(g, buffer);
    
	conf_gpu x_orig(x, buffer);
	conf_gpu x_new(x, buffer);

	change_gpu p(g, buffer);
    change_gpu y(g, buffer);

    change_gpu minus_hy(g, buffer);
    float* f0;
    float out_energy;

    CUDA_CHECK_GNINA(cudaMalloc(&f0, sizeof(float)));
    //TODO: make safe for the case where num_movable_atoms > 1024
    assert(f.ig.num_movable_atoms <= 1024);
    bfgs_gpu<<<1,max(WARPSIZE,f.ig.num_movable_atoms)>>>(f,
                      x, x_orig, x_new,
                      g, g_orig, g_new,
                      p, y, h, minus_hy,
                      average_required_improvement, params, f0);
    sync_and_errcheck();
    CUDA_CHECK_GNINA(cudaFree(h.m_data));
    CUDA_CHECK_GNINA(cudaMemcpy(&out_energy,
                                f0, sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK_GNINA(cudaFree(f0));
	return out_energy;
}
