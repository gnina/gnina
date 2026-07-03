#include "quasi_newton.h"
#include "conf_gpu.h"
#include "matrix.h"
#include "bfgs.h"
#include "device_buffer.h"

#ifdef USE_METAL
  #include "cuda_metal_compat.h"
#else
  #include <cuda_runtime.h>
#endif

static __device__ fl new_lambdamin_val(fl p_val, fl x_val, fl cur_val) {
  // return fabsf(p.values[i]) / fmaxf(fabsf(x.values[i]),
  //                                   1.0f);
  fl v = fabsf(p_val / fmaxf(fabsf(x_val), 1.0f));
  return fmaxf(v, cur_val);
}

__device__ fl compute_lambdamin(const change_gpu& p, const conf_gpu& x, sz n,
    const gpu_data& d)

    {
  const tree_gpu& tree = *d.treegpu;
  fl test = 0;
  for (sz i = 0; i < tree.num_nodes; i++) {
    if (i < tree.nlig_roots) {
      float* conf_start = &x.values[i * 7];
      float* change_start = &p.values[i * 6];

      for (sz i = 0; i < 3; i++)
        test = new_lambdamin_val(change_start[i], conf_start[i], test);

      qt orientation(conf_start[3], conf_start[4], conf_start[5],
          conf_start[6]);
      // qt *orientation = (qt *)(void*) &conf_start[3];
      vec angle = quaternion_to_angle(orientation);
      for (sz i = 0; i < 3; i++)
        test = new_lambdamin_val(change_start[3 + i], angle[i], test);

    } else {
      float* conf_start = &x.values[i + tree.nlig_roots * 6];
      float* change_start = &p.values[i + tree.nlig_roots * 5];
      test = new_lambdamin_val(*change_start, *conf_start, test);
    }

    // fl temp = fabsf(p.values[i]) / fmaxf(fabsf(x.values[i]),
    //                                      1.0f);
    // if (temp > test)
    // 	test = temp;
  }
  return test;
}

//TODO: operator -=
__device__ inline
void subtract_change(change_gpu& b, const change_gpu& a, sz n) { // b -= a
  b.sub(a);
}

__device__ inline
void set_to_neg(change_gpu& b, const change_gpu& a, sz n) { // b = -a
  b = a;
  b.invert();
}

__device__
void set_diagonal(flmat_gpu& m, fl x) {
  VINA_FOR(i, m.dim())
    m(i, i) = x;
}

__device__  inline fl scalar_product(const change_gpu& a, const change_gpu& b,
    sz n) {
  return a.dot(b);
}

__device__ inline void minus_mat_vec_product(const flmat_gpu& m,
    const change_gpu& in, change_gpu& out) {
  in.minus_mat_vec_product(m, out);
}

template<typename infoT>
__device__ fl fast_line_search(quasi_newton_aux_gpu<infoT>& f, sz n,
    const conf_gpu& x, const change_gpu& g, const fl f0, const change_gpu& p,
    conf_gpu& x_new, change_gpu& g_new, fl& f1) { // returns alpha
  const fl c0 = 0.0001;
  const unsigned max_trials = 10;
  const fl multiplier = 0.5;
  fl alpha = 1;
  int idx = threadIdx.x;

  const fl pg = scalar_product(p, g, n);

  VINA_U_FOR(trial, max_trials) {
    if (idx == 0) x_new = x;
    __syncthreads();
    if (idx < x_new.n) x_new.increment(p, alpha, &f.gdata);
    __syncthreads();
    if (idx == 0) f1 = f(x_new, g_new);
    __syncthreads();
    if (f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
    break;
    alpha *= multiplier;
  }
  return alpha;
}

template<typename infoT>
__device__ fl accurate_line_search_gpu(quasi_newton_aux_gpu<infoT>& f, sz n,
    const conf_gpu& x, const change_gpu& g, const fl f0, const change_gpu& p,
    conf_gpu& x_new, change_gpu& g_new, fl& f1) {
  fl a, alpha2 = 0, b, disc, f2 = 0;
  fl rhs1, rhs2, slope = 0, test, tmplam;
  const fl ALF = 1.0e-4;
  const fl FIRST = 1.0;

  __shared__ fl alpha;
  __shared__ fl alamin;
  int idx = threadIdx.x;
  slope = scalar_product(g, p, n);
  if (slope >= 0) {
    //gradient isn't actually in a decreasing direction
    if (idx == 0) {
      x_new = x;
      g_new.clear(); //dkoes - set gradient to zero
    }
    return 0;
  }
  if (idx == 0) {
    test = compute_lambdamin(p, x, n, f.gdata);

    alamin = epsilon_fl / test;
    alpha = FIRST; //single newton step
  }
  for (;;) //always try full newton step first
      {
    if (idx == 0) {
      x_new = x;
    }
    __syncthreads();
    if (idx < x_new.n) x_new.increment(p, alpha, &f.gdata);

    __syncthreads();
    if (idx == 0) f1 = f(x_new, g_new);

    __syncthreads();
    //std::cout << "alpha " << alpha << "  f " << f1 << "\tslope " << slope << " f0ALF " << f0 + ALF * alpha * slope << "\n";
    if (alpha < alamin) //convergence
        {
      if (idx == 0) {
        x_new = x;
        g_new.clear(); //dkoes - set gradient to zero
      }
      return 0;
    } else
      if (f1 <= f0 + ALF * alpha * slope) {
        //sufficient function decrease, stop searching
        return alpha;
      } else //have to backtrack
      {
        if (idx == 0) {
          if (alpha == FIRST) {
            //first time
            tmplam = -slope / (2.0 * (f1 - f0 - slope));
          } else //subsequent backtracks
          {
            rhs1 = f1 - f0 - alpha * slope;
            rhs2 = f2 - f0 - alpha2 * slope;
            a = (rhs1 / (alpha * alpha) - rhs2 / (alpha2 * alpha2))
                / (alpha - alpha2);
            b = (-alpha2 * rhs1 / (alpha * alpha)
                + alpha * rhs2 / (alpha2 * alpha2)) / (alpha - alpha2);
            if (a == 0.0)
              tmplam = -slope / (2.0 * b);
            else {
              disc = b * b - 3.0 * a * slope;
              if (disc < 0)
                tmplam = 0.5 * alpha;
              else
                if (b <= 0)
                  tmplam = (-b + sqrt(disc)) / (3.0 * a);
                else
                  tmplam = -slope / (b + sqrt(disc));
            }
            if (tmplam > .5 * alpha) tmplam = .5 * alpha; //always at least cut in half
          }
        }
      }
    if (idx == 0) {
      alpha2 = alpha;
      f2 = f1;
      //std::cout << "TMPLAM " << tmplam << "\n";
      alpha = fmaxf(tmplam, (fl) 0.1 * alpha); //never smaller than a tenth
    }
  }

  return 0; // absolutely necessary to workaround nvcc compiler bug!!! (only took N days to find..)

}

__device__
void bfgs_update(flmat_gpu& h, const change_gpu& p, const change_gpu& y,
    const fl alpha, change_gpu &minus_hy) {
  const fl yp = y.dot(p);
  const sz n = p.num_floats();
  int idx = threadIdx.x;
  if (alpha * yp < epsilon_fl) return; // FIXME?

  if (idx == 0) minus_hy = y;

  if (idx < minus_hy.n) y.minus_mat_vec_product(h, minus_hy);

  __syncthreads();

  const fl yhy = -y.dot(minus_hy);
  if (idx < n) {
    const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon

    float coef = +alpha * alpha * (r * r * yhy + r);

    float *minus_hyvec = minus_hy.values;
    float *pvec = p.values;
    VINA_RANGE(j, idx, n) // includes i
      atomicAdd(&h(idx, j),
          alpha * r * (minus_hyvec[idx] * pvec[j] + minus_hyvec[j] * pvec[idx])
              + coef * pvec[idx] * pvec[j]);
  }
  __syncthreads();
  // s * s == alpha * alpha * p * p	} *
}

#ifdef USE_METAL
// CPU-serial BFGS helpers for Metal / Apple Silicon.
// conf_gpu / change_gpu value arrays live in MTLStorageModeShared unified
// memory and are CPU-accessible. The __device__ functions use threadIdx.x
// for intra-warp parallelism; these serial replacements are used instead.

static inline fl cpu_dot(const change_gpu& a, const change_gpu& b) {
  fl s = 0;
  for (sz i = 0; i < a.n; i++) s += (fl) a.values[i] * b.values[i];
  return s;
}

static void cpu_sub(change_gpu& b, const change_gpu& a) { // b -= a
  for (sz i = 0; i < a.n; i++) b.values[i] -= a.values[i];
}

static void cpu_minus_mvp(const flmat_gpu& m, const change_gpu& in,
    change_gpu& out) {
  for (sz idx = 0; idx < in.n; idx++) {
    fl sum = 0;
    for (sz j = 0; j < in.n; j++)
      sum += m(m.index_permissive(idx, j)) * in.values[j];
    out.values[idx] = (float) -sum;
  }
}

static void cpu_increment(conf_gpu& x, const change_gpu& c, fl factor,
    gpu_data& gdata) {
  tree_gpu& tree = *gdata.treegpu;
  unsigned lig_roots = tree.nlig_roots;
  for (unsigned idx = 0; idx < lig_roots; idx++) {
    unsigned conf_offset   = idx * 7;
    unsigned change_offset = idx * 6;
    for (int i = 0; i < 3; i++)
      x.values[conf_offset + i] += (float) (c.values[change_offset + i] * factor);
    qt orientation(x.values[conf_offset + 3], x.values[conf_offset + 4],
        x.values[conf_offset + 5], x.values[conf_offset + 6]);
    vec rotation((fl) factor * c.values[change_offset + 3],
        (fl) factor * c.values[change_offset + 4],
        (fl) factor * c.values[change_offset + 5]);
    quaternion_increment(orientation, rotation);
    x.values[conf_offset + 3] = (float) orientation.R_component_1();
    x.values[conf_offset + 4] = (float) orientation.R_component_2();
    x.values[conf_offset + 5] = (float) orientation.R_component_3();
    x.values[conf_offset + 6] = (float) orientation.R_component_4();
  }
  for (unsigned idx = lig_roots; idx < tree.num_nodes; idx++) {
    unsigned conf_offset   = idx + (6 * lig_roots);
    unsigned change_offset = idx + (5 * lig_roots);
    x.values[conf_offset] += (float) normalized_angle(
        (fl) factor * c.values[change_offset]);
    normalize_angle(x.values[conf_offset]);
  }
}

static void cpu_bfgs_update(flmat_gpu& h, const change_gpu& p,
    const change_gpu& y, fl alpha, change_gpu& minus_hy) {
  const fl yp = cpu_dot(y, p);
  const sz n  = p.n;
  if (alpha * yp < epsilon_fl) return;

  minus_hy = y; // copy via operator= (memcpy under Metal)
  cpu_minus_mvp(h, y, minus_hy); // minus_hy = -H * y

  const fl yhy  = -cpu_dot(y, minus_hy); // y^T * H * y
  const fl r    = 1.0 / (alpha * yp);
  const fl coef = alpha * alpha * (r * r * yhy + r);
  for (sz i = 0; i < n; i++) {
    for (sz j = i; j < n; j++) {
      h(i, j) += (fl) (alpha * r
          * (minus_hy.values[i] * p.values[j]
              + minus_hy.values[j] * p.values[i])
          + coef * p.values[i] * p.values[j]);
    }
  }
}

template<typename infoT>
static fl cpu_fast_line_search(quasi_newton_aux_gpu<infoT>& f, sz n,
    const conf_gpu& x, const change_gpu& g, fl f0, const change_gpu& p,
    conf_gpu& x_new, change_gpu& g_new, fl& f1) {
  const fl c0 = 0.0001;
  const unsigned max_trials = 10;
  const fl multiplier = 0.5;
  fl alpha = 1.0;
  const fl pg = cpu_dot(p, g);

  VINA_U_FOR(trial, max_trials) {
    x_new = x;
    cpu_increment(x_new, p, alpha, f.gdata);
    f1 = f(x_new, g_new);
    if (f1 - f0 < c0 * alpha * pg) break;
    alpha *= multiplier;
  }
  return alpha;
}

template<typename infoT>
static fl cpu_accurate_line_search(quasi_newton_aux_gpu<infoT>& f, sz n,
    const conf_gpu& x, const change_gpu& g, fl f0, const change_gpu& p,
    conf_gpu& x_new, change_gpu& g_new, fl& f1) {
  fl alpha2 = 0, f2 = 0;
  fl rhs1, rhs2, tmplam;
  const fl ALF   = 1.0e-4;
  const fl FIRST = 1.0;
  fl alpha = FIRST;
  fl alamin;

  fl slope = cpu_dot(g, p);
  if (slope >= 0) {
    x_new = x;
    g_new.clear();
    return 0;
  }

  fl test = compute_lambdamin(p, x, n, f.gdata);
  alamin = epsilon_fl / test;

  for (;;) {
    x_new = x;
    cpu_increment(x_new, p, alpha, f.gdata);
    f1 = f(x_new, g_new);

    if (alpha < alamin) {
      x_new = x;
      g_new.clear();
      return 0;
    } else if (f1 <= f0 + ALF * alpha * slope) {
      return alpha;
    } else {
      fl a, b, disc;
      if (alpha == FIRST) {
        tmplam = -slope / (2.0 * (f1 - f0 - slope));
      } else {
        rhs1 = f1 - f0 - alpha * slope;
        rhs2 = f2 - f0 - alpha2 * slope;
        a = (rhs1 / (alpha * alpha) - rhs2 / (alpha2 * alpha2))
            / (alpha - alpha2);
        b = (-alpha2 * rhs1 / (alpha * alpha)
            + alpha * rhs2 / (alpha2 * alpha2)) / (alpha - alpha2);
        if (a == 0.0)
          tmplam = -slope / (2.0 * b);
        else {
          disc = b * b - 3.0 * a * slope;
          if (disc < 0)
            tmplam = 0.5 * alpha;
          else if (b <= 0)
            tmplam = (-b + sqrt(disc)) / (3.0 * a);
          else
            tmplam = -slope / (b + sqrt(disc));
        }
        if (tmplam > 0.5 * alpha) tmplam = 0.5 * alpha;
      }
      alpha2 = alpha;
      f2     = f1;
      alpha  = fmaxf((float) tmplam, 0.1f * (float) alpha);
    }
  }
  return 0;
}

template<typename infoT>
static fl bfgs_metal(quasi_newton_aux_gpu<infoT>& f, conf_gpu& x, change_gpu& g,
    const fl /*average_required_improvement*/, const minimization_params& params) {
  sz n = g.num_floats();
  flmat_gpu h(n);

  change_gpu g_orig(g, thread_buffer);
  change_gpu g_new(g, thread_buffer);
  conf_gpu   x_orig(x, thread_buffer);
  conf_gpu   x_new(x, thread_buffer);
  change_gpu p(g, thread_buffer);
  change_gpu y(g, thread_buffer);
  change_gpu minus_hy(g, thread_buffer);

  fl f0     = f(x, g);
  fl f_orig = f0;
  g_orig = g;
  x_orig = x;

  VINA_U_FOR(step, params.maxiters) {
    cpu_minus_mvp(h, g, p); // p = -H * g (search direction)

    fl f1 = 0, alpha;
    if (params.type == minimization_params::BFGSAccurateLineSearch)
      alpha = cpu_accurate_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
    else
      alpha = cpu_fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);

    if (alpha == 0) break;

    y = g_new;
    fl prevf0 = f0;
    f0 = f1;
    x  = x_new;

    cpu_sub(y, g); // y = g_new - g_old (gradient change)

    if (params.early_term && fabsf((float) (prevf0 - f0)) < 1e-5f) break;

    g = g_new;

    fl gradnormsq = cpu_dot(g, g);
    if (!(gradnormsq >= 1e-4)) break; // also breaks for NaN

    if (step == 0) {
      const fl yy = cpu_dot(y, y);
      if (fabsf((float) yy) > epsilon_fl) {
        const fl yp = cpu_dot(y, p);
        set_diagonal(h, alpha * yp / yy);
      }
    }

    cpu_bfgs_update(h, p, y, alpha, minus_hy);
  }

  if (!(f0 <= f_orig)) { // restore best if no improvement (also for NaN)
    f0 = f_orig;
    x  = x_orig;
    g  = g_orig;
  }
  return f0;
}
#endif // USE_METAL

#ifndef USE_METAL  // ── CUDA kernel ──
template<typename infoT>
__global__
void bfgs_gpu(quasi_newton_aux_gpu<infoT> f, conf_gpu x, conf_gpu x_orig,
    conf_gpu x_new, change_gpu g, change_gpu g_orig, change_gpu g_new,
    change_gpu p, change_gpu y, flmat_gpu h, change_gpu minus_hy,
    const fl average_required_improvement, const minimization_params params,
    float* out_energy) {
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
    g_orig = g;
    x_orig = x; 
    p = g;
  }
  VINA_U_FOR(step, params.maxiters) {
    if (idx < g.n) {
      minus_mat_vec_product(h, g, p);
      // f1 is the returned energy for the next iteration of eval_deriv_gpu
      f1 = 0;
    }
    __syncthreads();
    if (params.type == minimization_params::BFGSAccurateLineSearch)
      alpha = accurate_line_search_gpu(f, n, x, g, f0, p, x_new, g_new, f1);
    else
      alpha = fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
    if (alpha == 0) break;
    fl prevf0;

    if (idx == 0) {
      y = g_new;

      prevf0 = f0;
      f0 = f1;
      x = x_new;
    }

    // Update line direction
    if (idx < y.n) subtract_change(y, g, n);

    if (params.early_term) {
      if (idx == 0) diff = prevf0 - f0;
      __syncthreads();
      if (fabsf(diff) < 1e-5) break;
    }

    if (idx == 0) g = g_new;

    __syncthreads();
    fl gradnormsq = scalar_product(g, g, n);
//		std::cout << "step " << step << " " << f0 << " " << gradnormsq << " " << alpha << "\n";

    if (!(gradnormsq >= 1e-4)) //slightly arbitrary cutoff - works with fp
      break;// breaks for nans too // FIXME !!??

    if (step == 0) {
      const fl yy = scalar_product(y, y, n);
      if (fabsf(yy) > epsilon_fl) {
        const fl yp = scalar_product(y, p, n);
        if (idx == 0) set_diagonal(h, alpha * yp / yy);
      }
    }
    // bfgs_update used to return a bool, but the value of that bool never
    // got checked anyway
    bfgs_update(h, p, y, alpha, minus_hy);
  }
  if (idx == 0) {
    if (!(f0 <= f_orig)) { // succeeds for nans too
      f0 = f_orig;
      x = x_orig;
      g = g_orig;
    }
    *out_energy = f0;
  }
}
#endif // USE_METAL

template<typename infoT>
fl bfgs(quasi_newton_aux_gpu<infoT> &f, conf_gpu& x, change_gpu& g,
    const fl average_required_improvement, const minimization_params& params) {
  sz n = g.num_floats();

  // Initialize and copy Hessian
  flmat_gpu h(n);

  // Initialize and copy additional conf and change objects
  change_gpu g_orig(g, thread_buffer);
  change_gpu g_new(g, thread_buffer);

  conf_gpu x_orig(x, thread_buffer);
  conf_gpu x_new(x, thread_buffer);

  change_gpu p(g, thread_buffer);
  change_gpu y(g, thread_buffer);

  change_gpu minus_hy(g, thread_buffer);
  float* f0;
  float out_energy;

#ifndef USE_METAL
  CUDA_CHECK_GNINA(device_malloc(&f0, sizeof(float)));
  //TODO: make safe for the case where num_movable_atoms > 1024
  assert(f.ig.num_movable_atoms <= 1024);
  bfgs_gpu<<<1, ROUND_TO_WARP(max(WARPSIZE, f.ig.num_movable_atoms))>>>(f, x,
      x_orig, x_new, g, g_orig, g_new, p, y, h, minus_hy,
      average_required_improvement, params, f0);
  CUDA_CHECK_GNINA(
      definitelyPinnedMemcpy(&out_energy, f0, sizeof(float),
          cudaMemcpyDeviceToHost));
  CUDA_CHECK_GNINA(device_free(f0));
#else
  (void)f0;
  out_energy = bfgs_metal(f, x, g, average_required_improvement, params);
#endif
  return out_energy;
}

template __device__ fl accurate_line_search_gpu(
    quasi_newton_aux_gpu<GPUNonCacheInfo>&, sz, const conf_gpu&,
    const change_gpu&, const fl, const change_gpu&, conf_gpu&, change_gpu&,
    fl&);
template __device__ fl accurate_line_search_gpu(
    quasi_newton_aux_gpu<GPUCacheInfo>&, sz, const conf_gpu&, const change_gpu&,
    const fl, const change_gpu&, conf_gpu&, change_gpu&, fl&);

template __device__ fl fast_line_search(quasi_newton_aux_gpu<GPUNonCacheInfo>&,
    sz, const conf_gpu&, const change_gpu&, const fl, const change_gpu&,
    conf_gpu&, change_gpu&, fl&);
template __device__ fl fast_line_search(quasi_newton_aux_gpu<GPUCacheInfo>&, sz,
    const conf_gpu&, const change_gpu&, const fl, const change_gpu&, conf_gpu&,
    change_gpu&, fl&);

#ifndef USE_METAL
template __global__ void bfgs_gpu(quasi_newton_aux_gpu<GPUNonCacheInfo>,
    conf_gpu x, conf_gpu x_orig, conf_gpu x_new, change_gpu g,
    change_gpu g_orig, change_gpu g_new, change_gpu p, change_gpu y,
    flmat_gpu h, change_gpu minus_hy, const fl average_required_improvement,
    const minimization_params params, float* out_energy);
template __global__ void bfgs_gpu(quasi_newton_aux_gpu<GPUCacheInfo>,
    conf_gpu x, conf_gpu x_orig, conf_gpu x_new, change_gpu g,
    change_gpu g_orig, change_gpu g_new, change_gpu p, change_gpu y,
    flmat_gpu h, change_gpu minus_hy, const fl average_required_improvement,
    const minimization_params params, float* out_energy);
#endif // USE_METAL

template fl bfgs(quasi_newton_aux_gpu<GPUNonCacheInfo>&, conf_gpu&, change_gpu&,
    const fl, const minimization_params&);
template fl bfgs(quasi_newton_aux_gpu<GPUCacheInfo>&, conf_gpu&, change_gpu&,
    const fl, const minimization_params&);
