#pragma once 

#include <cmath>
#include <iostream>
#include "../lib/matrix.h"

struct flt_int {
  float val;
  int idx;
};

template<typename T> inline T __device__ __host__ zero(void);

template<> inline flt_int zero(void) {
  flt_int tmp = {0., 0};
  return tmp;
}

inline flt_int __device__ __host__  min(void) {
  flt_int tmp = { -100., -1 };
  return tmp;
}

__host__ __device__ inline bool operator==(const flt_int& lhs, const flt_int& rhs){ return lhs.val==rhs.val && lhs.idx == rhs.idx; }
__host__ __device__ inline bool operator!=(const flt_int& lhs, const flt_int& rhs){return !operator==(lhs,rhs);}
__host__ __device__ inline bool operator< (const flt_int& lhs, const flt_int& rhs){ return lhs.val < rhs.val ? true: false; }
__host__ __device__ inline bool operator> (const flt_int& lhs, const flt_int& rhs){return  operator< (rhs,lhs);}
__host__ __device__ inline bool operator<=(const flt_int& lhs, const flt_int& rhs){return !operator> (lhs,rhs);}
__host__ __device__ inline bool operator>=(const flt_int& lhs, const flt_int& rhs){return !operator< (lhs,rhs);}

#ifdef __CUDACC__

__device__ __inline__ flt_int shuffle_down(flt_int pair, int offset) {
  flt_int tmp;
#if __CUDACC_VER_MAJOR__ >= 9
  tmp.val = __shfl_down_sync(0xffffffff, pair.val, offset);
  tmp.idx = __shfl_down_sync(0xffffffff, pair.idx, offset);
#else
  tmp.val = __shfl_down(pair.val,offset);
  tmp.idx = __shfl_down(pair.idx,offset);
#endif
  return tmp;
}

#endif

inline const flt_int& max( const flt_int& a, const flt_int& b) {
      return a.val > b.val ? a : b;
}

#pragma omp declare reduction( maxVal: flt_int: omp_out=max( omp_out, omp_in ) )

inline void cpu_l1(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float diff = fabs(optgrid[k] - screengrid[k]);
    sum += diff;
  }
  *scoregrid = sum;
}

inline void cpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float diff = optgrid[k] - screengrid[k];
    float sqdiff = diff * diff;
    sum += sqdiff;
  }
  *scoregrid = sum;
}

inline void cpu_mult(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float weight = optgrid[k] * screengrid[k];
    sum += weight;
  }
  *scoregrid = sum;
}

inline void cpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid, 
    size_t gsize, float positive_threshold, float negative_threshold) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float threshold = optgrid[k] >=0 ? positive_threshold : negative_threshold;
    float sign = optgrid[k] >= 0 ? 1 : -1;
    float magnitude = fabs(optgrid[k]);
    float weight = ((magnitude > threshold) && screengrid[k]) ? 1 : 0;
    sum += sign * weight;
  }
  *scoregrid = sum;
}

inline void cpu_calcSig(const float* grid, std::vector<float>& sig, unsigned subgrid_dim, 
    unsigned dim, unsigned gsize, unsigned ntypes, unsigned blocks_per_side) {
  // perform reduction to obtain total weight for each cube, normalized...
  // this might not be the "right" signature, but i think we basically want to
  // say, with reasonably good granularity, "this is how much weight is here."
  // TODO: redo this with boost::multi_array_ref to sanity check the indexing
  unsigned n_subcubes = blocks_per_side * blocks_per_side * blocks_per_side;
  unsigned subcube_npts = subgrid_dim * subgrid_dim * subgrid_dim;
  unsigned siglength = n_subcubes * ntypes;
  if (sig.size() != siglength) {
    sig.resize(siglength);
  }

  for (size_t ch=0; ch<ntypes; ++ch) {
#pragma omp parallel for
    for (size_t cube_idx=0; cube_idx<n_subcubes; ++cube_idx) {
      float total = 0;
      unsigned block_x_offset = cube_idx / (blocks_per_side * blocks_per_side);
      unsigned block_y_offset = (cube_idx % (blocks_per_side * blocks_per_side)) / blocks_per_side;
      unsigned block_z_offset = cube_idx % blocks_per_side;
      unsigned cube_offset = block_x_offset * (blocks_per_side * blocks_per_side * 
          subgrid_dim * subgrid_dim * subgrid_dim) + 
          block_y_offset * (blocks_per_side * subgrid_dim * subgrid_dim) + 
          block_z_offset * subgrid_dim;
      // omp simd?
      for (size_t idx=0; idx<subcube_npts; ++idx) {
        unsigned thread_x_offset = idx / (subgrid_dim * subgrid_dim);
        unsigned thread_y_offset = (idx % (subgrid_dim * subgrid_dim)) / subgrid_dim;
        unsigned thread_z_offset = idx % subgrid_dim;
        unsigned tidx = cube_offset + thread_x_offset * (dim*dim) + thread_y_offset *
          dim + thread_z_offset;
        total += grid[ch * gsize + tidx];
      }
      sig[ch * n_subcubes + cube_idx] = total;
    }
  }

  // that's the total amt of density per cube, now let's convert to weights
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t i=0; i<sig.size(); ++i) {
    sum += sig[i];
  }
#pragma omp parallel
  for (size_t i=0; i<sig.size(); ++i) {
    sig[i] /= sum;
  }
}
  
inline void cpu_emd(std::vector<float>& optsig, std::vector<float>& screensig, 
    float* scoregrid, unsigned dim, unsigned subgrid_dim, unsigned ntypes, 
    unsigned gsize, flmat& cost_matrix, std::vector<float>& flow) {
  unsigned maxiter = 10000;
  float reg = 9;
  double tolerance = 1e-9;
  double current_threshval = 1;
  unsigned siglength = optsig.size();
  unsigned sigsq = siglength * siglength;

  // right now require signatures to be the same length, they always should be
  // anyway because they depend on the CNN grid size
  if (screensig.size() != siglength) {
    std::cerr << "EMD signatures need to be the same length, pad if necessary";
    std::exit(1);
  }
  // sanity check cost matrix too
  if (cost_matrix.dim() != siglength) {
    std::cerr << "Cost matrix dimension doesn't match signatures";
    std::exit(1);
  }
  // unlike cost, flow isn't symmetric
  if (flow.size() != sigsq) {
    std::cerr << "Flow matrix is wrong size";
    std::exit(1);
  }
  if (siglength == 0) {
    std::cerr << "Signature is empty";
    std::exit(1);
  }
  flmat K(siglength,0); 
#pragma omp parallel for
  for (size_t i=0; i<siglength; ++i) {
    K(i) = std::exp(cost_matrix(i) / -reg);
  }

  std::vector<float> u(siglength, 1.f/siglength);
  std::vector<float> v(siglength, 1.f/siglength);
#pragma omp parallel for
  for (size_t i=0; i<siglength; ++i) {
    for (size_t j=0; j<siglength; ++j) {
      flow[i*siglength + j] = u[i] * K(i,j) * v[j];
    }
  }

  std::vector<float> viol(siglength, 0);
  std::vector<float> viol_2(siglength, 0);
  // TODO: verify the inner loop is actually being vectorized
  float sum = 0;
#pragma omp parallel for
  for (size_t i=0; i<siglength; ++i) {
#pragma omp simd reduction(+: sum)
    for (size_t j=0; j<siglength; ++j) {
      sum += flow[i*siglength+j];
    }
    viol[i] = sum - optsig[i];
    sum = 0;
  }

#pragma omp parallel for
  for (size_t i=0; i<siglength; ++i) {
#pragma omp simd reduction(+: sum)
    for (size_t j=0; j<siglength; ++j) {
      sum += flow[j*siglength+i];
    }
    viol_2[i] = sum - screensig[i];
    sum = 0;
  }

  for (size_t iter=0; iter<maxiter; ++iter) {
    flt_int argmax_1 = { -100., -1 };
#pragma omp parallel for reduction( maxVal: argmax_1 )
    for (size_t i=0; i<siglength; ++i) {
      float val = std::abs(viol[i]);
      if (val > argmax_1.val) {
        argmax_1.idx = i;
        argmax_1.val = val;
      }
    }
    unsigned i_1 = argmax_1.idx;
    float m_viol_1 = std::abs(viol[i_1]);

    flt_int argmax_2 = { -100., -1 };
#pragma omp parallel for reduction( maxVal: argmax_2 )
    for (size_t i=0; i<siglength; ++i) {
      float val = std::abs(viol_2[i]);
      if (val > argmax_2.val) {
        argmax_2.idx = i;
        argmax_2.val = val;
      }
    }
    unsigned i_2 = argmax_2.idx;
    float m_viol_2 = std::abs(viol_2[i_2]);
    current_threshval = std::max(m_viol_1, m_viol_2);

    if (m_viol_1 > m_viol_2) {
      float old_u = u[i_1];
      float K_dot_v = 0;
#pragma omp parallel for reduction(+: K_dot_v)
      for (size_t i=0; i<siglength; ++i) {
        K_dot_v += K(i_1 * siglength + i) * v[i];
      }
      u[i_1] = optsig[i_1] / K_dot_v;
      float udiff = u[i_1] - old_u;

#pragma omp parallel for
      for (size_t i=0; i<siglength; ++i) {
        flow[i_1 * siglength + i] = u[i_1] * K(i_1 * siglength + i) * v[i];
        viol_2[i] += K(i_1*siglength + i) * v[i] * udiff;
      }
      viol[i_1] = u[i_1] * K_dot_v - optsig[i_1];
    }
    else {
      float old_v = v[i_2];
      float K_dot_u = 0;
#pragma omp parallel for reduction(+: K_dot_u)
      for (size_t i=0; i<siglength; ++i) {
        K_dot_u += K(i*siglength+i_2) * u[i];
      }
      v[i_2] = screensig[i_2] / K_dot_u;
      float vdiff = v[i_2] - old_v;

#pragma omp parallel for
      for (size_t i=0; i<siglength; ++i) {
        float Kval = K(i*siglength + i_2);
        flow[i*siglength + i_2] = u[i] * Kval * v[i_2];
        viol[i] += vdiff * Kval * u[i];
      }
      viol_2[i_2] = v[i_2] * K_dot_u - screensig[i_2];
    }
    
    if (current_threshval <= tolerance)
      break;
    else if (iter == maxiter-1)
      std::cout << "Warning: EMD did not converge"; // everything user-facing calls this EMD
  }

  // now we have the flow and can do \sum F \odot D to get the total cost
  float cost = 0;
#pragma omp parallel for reduction(+: cost)
  for (size_t i=0; i<sigsq; ++i) {
    // gross...
    unsigned idx1 = i / siglength;
    unsigned idx2 = i % siglength;
    cost += cost_matrix(idx1, idx2) * flow[i];
  }
  *scoregrid = cost;
}

inline void do_cpu_emd(const float* optgrid, const float* screengrid, float* scoregrid, 
    unsigned dim, unsigned subgrid_dim, unsigned blocks_per_side, unsigned ntypes, 
    size_t gsize, flmat& cost_matrix) {
  // calculates the Greenkhorn distance, a near-linear time approximation to
  // the entropy-regularized EMD, for the pair of signatures optsig and screensig, 
  // using user-provided cost_matrix. populates the flow matrix for
  // visualization purposes and sets score to the final transportation cost
  // TODO: dump the flow matrix someplace
  std::vector<float> optsig;
  std::vector<float> screensig;
  std::vector<float> flow;

  // get signatures for each grid
  cpu_calcSig(optgrid, optsig, subgrid_dim, dim, gsize, ntypes, blocks_per_side);
  cpu_calcSig(screengrid, screensig, subgrid_dim, dim, gsize, ntypes, blocks_per_side);

  // calculate flow and total cost
  cpu_emd(optsig, screensig, scoregrid, dim, subgrid_dim, ntypes, gsize, cost_matrix, flow);
}

inline float l1(const vec& icoords, const vec& jcoords) {
  float sum = 0.;
  for (size_t k=0; k<3; ++k) {
    sum += icoords[k] - jcoords[k];
  }
  return sum;
}

inline vec get_cube_coords(unsigned cube_idx, float dimension, unsigned blocks_per_side) {
  // wlog have cube(0,0,0) start at (0,0,0) (i.e. ignore grid center) since
  // distance is invariant to translation
  vec coords(0,0,0);
  // map flat cube index to (x,y,z) offsets
  unsigned x_offset = cube_idx / (blocks_per_side * blocks_per_side);
  unsigned y_offset = (cube_idx % (blocks_per_side * blocks_per_side)) / blocks_per_side;
  unsigned z_offset = cube_idx % blocks_per_side;
  float cube_dimension = dimension / blocks_per_side;
  coords[0] = x_offset * cube_dimension;
  coords[1] = y_offset * cube_dimension;
  coords[2] = z_offset * cube_dimension;
  return coords;
}

inline float get_cube_l1(unsigned i_cube, unsigned j_cube, float dimension, unsigned
    blocks_per_side) {
  // get (x,y,z) coords of cubes
  vec icoords = get_cube_coords(i_cube, dimension, blocks_per_side);
  vec jcoords = get_cube_coords(j_cube, dimension, blocks_per_side);
  return l1(icoords, jcoords);
}

inline void populate_cost_matrix(unsigned ncubes, unsigned subgrid_dim,
    unsigned ntypes, unsigned blocks_per_side, float dimension, flmat& cost_matrix) {
  // fill in cost matrix for comparing two spatial signatures
  unsigned n_elements = ntypes * ncubes;
  for (unsigned i=0; i<n_elements; ++i) {
    for (unsigned j=i+1; j<n_elements; ++j) {
      // within the same cube, cost is 0...assumes 0-initialized
      if (i!=j) {
        unsigned i_ch = i / ncubes;
        unsigned j_ch = j / ncubes;
        unsigned i_cube = i % ncubes;
        unsigned j_cube = j % ncubes;
        float cube_l1 = get_cube_l1(i_cube, j_cube, dimension, blocks_per_side);
        // within the same channel but different cube, cost is L1 distance between cubes
        if (i_ch == j_ch)
          cost_matrix(i,j) = cube_l1;
        // across channels, cost is constant, arbitrarily large value + L1
        else
          cost_matrix(i,j) = cube_l1 + 1e6;
      }
    }
  }
}

void do_gpu_l1(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_mult(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize, float positive_threshold, float negative_threshold);

void do_gpu_emd(const float* optgrid, const float* screengrid, float* scoregrid, 
    unsigned dim, unsigned subgrid_dim, unsigned blocks_per_side, unsigned ntypes, 
    size_t gsize, flmat& cost_matrix);
