#include "../lib/gpu_math.h"
#include "../lib/matrix.h"
#include "loss.h"
#define CUDA_NUM_THREADS 256
#define CUDA_NUM_BLOCKS 48
#define warpSize 32
#define KERNEL_CHECK CUDA_CHECK_GNINA(cudaPeekAtLastError())
#define CUDA_KERNEL_LOOP(i, n) \
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; \
       i < (n); \
       i += blockDim.x * gridDim.x)

template<typename T> T __device__ __host__ zero(void);
template<> float3 zero(void) {
  return float3(0, 0, 0);
}

template<> float zero(void) {
  return 0;
}

template<> flt_int zero(void) {
  flt_int tmp = {0, 0};
  return tmp;
}

flt_int __device__ __host__  min(void) {
  flt_int tmp = { -100., -1 };
  return tmp;
}

// device functions for warp-based reduction using shufl operations
// TODO: should probably just be factored out into gpu_math or gpu_util
template<class T>
__device__   __forceinline__ T warp_sum(T mySum) {
  for (int offset = warpSize >> 1; offset > 0; offset >>= 1)
    mySum += shuffle_down(mySum, offset);
  return mySum;
}

template<class T>
__device__   __forceinline__ T warp_max(T myMax) {
  for (int offset = warpSize >> 1; offset > 0; offset >>= 1) {
    T nextMax = shuffle_down(myMax, offset); 
    myMax = nextMax > myMax ? nextMax : myMax;
  }
  return myMax;
}

__device__ __forceinline__
bool isNotDiv32(unsigned int val) {
  return val & 31;
}

/* requires blockDim.x <= 1024, blockDim.y == 1 */
template<class T>
__device__   __forceinline__ T block_sum(T mySum) {
  const unsigned int lane = threadIdx.x & 31;
  const unsigned int wid = threadIdx.x >> 5;

  __shared__ T scratch[32];

  mySum = warp_sum(mySum);
  if (lane == 0) scratch[wid] = mySum;
  __syncthreads();

  if (wid == 0) {
    mySum = (threadIdx.x < blockDim.x >> 5) ? scratch[lane] : zero<T>();
    mySum = warp_sum(mySum);
    if (threadIdx.x == 0 && isNotDiv32(blockDim.x))
      mySum += scratch[blockDim.x >> 5];
  }
  return mySum;
}

/* requires blockDim.x <= 1024, blockDim.y == 1 */
template<class T>
__device__   __forceinline__ T block_max(T myMax) {
  const unsigned int lane = threadIdx.x & 31;
  const unsigned int wid = threadIdx.x >> 5;

  __shared__ T scratch[32];

  myMax = warp_max(myMax);
  if (lane == 0) scratch[wid] = myMax;
  __syncthreads();

  if (wid == 0) {
    myMax = (threadIdx.x < blockDim.x >> 5) ? scratch[lane] : min();
    myMax = warp_max(myMax);
    if (threadIdx.x == 0 && isNotDiv32(blockDim.x)) {
      T nextMax = scratch[blockDim.x >> 5]; 
      myMax = nextMax > myMax ? nextMax : myMax;
    }
  }
  return myMax;
}

__global__
void gpu_l1(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize) {
  unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned nthreads = blockDim.x * gridDim.x;
  // optimized grids
  float sum = 0.;
  for (size_t k=tidx; k<gsize; k+=nthreads) {
    sum += fabsf(optgrid[k] - screengrid[k]);
  }
  float total = block_sum<float>(sum);
  if (threadIdx.x == 0)
    atomicAdd(scoregrid, total);
}

void do_gpu_l1(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize) {
  unsigned block_multiple = gsize / CUDA_NUM_THREADS;
  unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
  gpu_l1<<<nblocks, CUDA_NUM_THREADS>>>(optgrid, screengrid, scoregrid, gsize);
  KERNEL_CHECK;
}

__global__
void gpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize) {
  unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned nthreads = blockDim.x * gridDim.x;
  // optimized grids
  float sum = 0.;
  for (size_t k=tidx; k<gsize; k+=nthreads) {
    float diff = optgrid[k] - screengrid[k];
    float sqdiff = diff * diff;
    sum += sqdiff;
  }
  float total = block_sum<float>(sum);
  if (threadIdx.x == 0)
    atomicAdd(scoregrid, total);
}

void do_gpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize) {
  unsigned block_multiple = gsize / CUDA_NUM_THREADS;
  unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
  gpu_l2sq<<<nblocks, CUDA_NUM_THREADS>>>(optgrid, screengrid, scoregrid, gsize);
  KERNEL_CHECK;
}

__global__
void gpu_mult(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize) {
  unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned nthreads = blockDim.x * gridDim.x;
  // optimized grids
  float sum = 0.;
  for (size_t k=tidx; k<gsize; k+=nthreads) {
    float weight = optgrid[k] * screengrid[k];
    sum += weight;
  }
  float total = block_sum<float>(sum);
  if (threadIdx.x == 0)
    atomicAdd(scoregrid, total);
}

void do_gpu_mult(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize) {
  unsigned block_multiple = gsize / CUDA_NUM_THREADS;
  unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
  gpu_mult<<<nblocks, CUDA_NUM_THREADS>>>(optgrid, screengrid, scoregrid, gsize);
  KERNEL_CHECK;
}

__global__
void gpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize, float positive_threshold, float negative_threshold) {
  unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned nthreads = blockDim.x * gridDim.x;
  // optimized grids
  float sum = 0.;
  for (size_t k=tidx; k<gsize; k+=nthreads) {
    float threshold = optgrid[k] >=0 ? positive_threshold : negative_threshold;
    float sign = optgrid[k] >= 0 ? 1 : -1;
    float magnitude = fabs(optgrid[k]);
    float weight = ((magnitude > threshold) && screengrid[k]) ? 1 : 0;
    sum += sign * weight;
  }
  float total = block_sum<float>(sum);
  if (threadIdx.x == 0)
    atomicAdd(scoregrid, total);
}

void do_gpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid,
    size_t gsize, float positive_threshold, float negative_threshold) {
  unsigned block_multiple = gsize / CUDA_NUM_THREADS;
  unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
  gpu_thresh<<<nblocks, CUDA_NUM_THREADS>>>(optgrid, screengrid, scoregrid,
      gsize, positive_threshold, negative_threshold);
  KERNEL_CHECK;
}

__global__
void gpu_calcSig(const float* grid, float* sig, unsigned subgrid_dim, unsigned dim, 
    unsigned gsize, unsigned ntypes, unsigned blocks_per_side) {
  // perform reduction to obtain total weight for each cube, normalized...
  // this might not be the "right" signature, but i think we basically want to
  // say, with reasonably good granularity, "this is how much weight is here."
  unsigned n_subcubes = blocks_per_side * blocks_per_side * blocks_per_side;
  unsigned ch_idx = blockIdx.x / n_subcubes;
  unsigned cube_idx = blockIdx.x % n_subcubes;
  unsigned block_x_offset = cube_idx / (blocks_per_side * blocks_per_side);
  unsigned block_y_offset = (cube_idx % (blocks_per_side * blocks_per_side)) / blocks_per_side;
  unsigned block_z_offset = cube_idx % blocks_per_side;
  unsigned cube_offset = block_x_offset * (blocks_per_side * blocks_per_side * 
      subgrid_dim * subgrid_dim * subgrid_dim) + 
      block_y_offset * (blocks_per_side * subgrid_dim * subgrid_dim) + 
      block_z_offset * subgrid_dim;
  unsigned thread_x_offset = threadIdx.x / (subgrid_dim * subgrid_dim);
  unsigned thread_y_offset = (threadIdx.x % (subgrid_dim * subgrid_dim)) / subgrid_dim;
  unsigned thread_z_offset = threadIdx.x % subgrid_dim;
  unsigned tidx = cube_offset + thread_x_offset * (dim*dim) + thread_y_offset *
    dim + thread_z_offset;
  unsigned overall_idx = ch_idx * gsize + tidx;
  float val = grid[overall_idx];
  float total = block_sum<float>(val);
  if (threadIdx.x == 0)
    sig[blockIdx.x] = total;
}

__global__
void gpu_normalizeSig(float* sig, unsigned sigsize, unsigned roundsize) {
  __shared__ float sigsum;
  float total;
  // can't just loop on sigsize because block_sum has __syncthreads
  CUDA_KERNEL_LOOP(i,roundsize) {
    float val = 0;
    if (i<sigsize)
      val = sig[i];
    total = block_sum<float>(val);
  }
  if (threadIdx.x == 0)
    sigsum = total;
  __syncthreads();

  CUDA_KERNEL_LOOP(i,sigsize) {
    sig[i] /= sigsum;
  }
}

__global__
void gpu_emd(float* optsig, float* screensig, float* scoregrid, unsigned dim, 
    unsigned subgrid_dim, unsigned ntypes, unsigned gsize, flmat_gpu cost_matrix, 
    unsigned siglength, float* flow, flmat_gpu K, float* u, float* v, 
    float* viol, float* viol_2) {
  unsigned maxiter = 10000;
  double tolerance = 1e-9;
  unsigned sigsq = siglength * siglength;
  __shared__ flt_int argmax_1;
  __shared__ flt_int argmax_2;
  __shared__ float K_dot_v;
  __shared__ float K_dot_u;

  // calculates the Greenkhorn distance, a near-linear time approximation to
  // the entropy-regularized EMD, for the pair of signatures optsig and screensig, 
  // using user-provided cost_matrix. populates the flow matrix for
  // visualization purposes and returns the final transportation cost in
  // location specified by scoregrid
  // unfortunately this isn't a perfect fit for the GPU because there are many
  // substeps that require explicit block-level synchronization

  // K, u and v are already populated; used them to populate flow
  for (size_t i=0; i<siglength; ++i) {
    CUDA_KERNEL_LOOP(j,siglength) {
      flow[i*siglength+j] = u[i] * K(i,j) * v[j];
    }
  }
  __syncthreads();

  // populate viol and viol_2
  float sum = 0;
  for (size_t i=0; i<siglength; ++i) {
    CUDA_KERNEL_LOOP(j,siglength) {
      sum += flow[i*siglength+j];
    }
    float total = block_sum<float>(sum);
    if (threadIdx.x == 0)
      viol[i] = total - optsig[i];
    sum = 0;
  }

  for (size_t i=0; i<siglength; ++i) {
    CUDA_KERNEL_LOOP(j,siglength) {
      sum += flow[j*siglength+i];
    }
    float total = block_sum<float>(sum);
    if (threadIdx.x == 0)
      viol_2[i] = total - screensig[i];
    sum = 0;
  }

  // iterate on flow updates until converged/hits maxiter
  for (size_t iter=0; iter<maxiter; ++iter) {
    flt_int my_argmax_1 = { -100., -1 };
    CUDA_KERNEL_LOOP(i,siglength) {
      float val = fabsf(viol[i]);
      if (val > argmax_1.val) {
        argmax_1.idx = i;
        argmax_1.val = val;
      }
    }
    my_argmax_1 = block_max<flt_int>(my_argmax_1);
    if (threadIdx.x == 0)
      argmax_1 = my_argmax_1;

    __syncthreads();
    unsigned i_1 = argmax_1.idx;
    float m_viol_1 = fabsf(viol[i_1]);

    flt_int my_argmax_2 = { -100., -1 };
    CUDA_KERNEL_LOOP(i,siglength) {
      float val = fabsf(viol_2[i]);
      if (val > my_argmax_2.val) {
        my_argmax_2.idx = i;
        my_argmax_2.val = val;
      }
    }
    my_argmax_2 = block_max<flt_int>(my_argmax_2);
    if (threadIdx.x == 0)
      argmax_2 = my_argmax_2;

    __syncthreads();
    unsigned i_2 = argmax_2.idx;
    float m_viol_2 = fabsf(viol_2[i_2]);
    double current_threshval = m_viol_1 > m_viol_2 ? m_viol_1 : m_viol_2;

    if (m_viol_1 > m_viol_2) {
      float old_u = u[i_1];
      float my_K_dot_v = 0;
      CUDA_KERNEL_LOOP(i,siglength) {
        my_K_dot_v += K(i_1 * siglength + i) * v[i];
      }
      my_K_dot_v = block_sum<float>(my_K_dot_v);
      if (threadIdx.x == 0)
        K_dot_v = my_K_dot_v;
      __syncthreads();

      u[i_1] = optsig[i_1] / K_dot_v;
      float udiff = u[i_1] - old_u;

      CUDA_KERNEL_LOOP(i,siglength) {
        flow[i_1 * siglength + i] = u[i_1] * K(i_1 * siglength + i) * v[i];
        viol_2[i] += K(i_1*siglength + i) * v[i] * udiff;
      }
      viol[i_1] = u[i_1] * K_dot_v - optsig[i_1];
    }
    else {
      float old_v = v[i_2];
      float my_K_dot_u = 0;
      CUDA_KERNEL_LOOP(i,siglength) {
        my_K_dot_u += K(i*siglength+i_2) * u[i];
      }
      my_K_dot_u = block_sum<float>(my_K_dot_u);
      if (threadIdx.x == 0)
        K_dot_u = my_K_dot_u;
      __syncthreads();

      v[i_2] = screensig[i_2] / K_dot_u;
      float vdiff = v[i_2] - old_v;

      CUDA_KERNEL_LOOP(i,siglength) {
        float Kval = K(i*siglength+ i_2);
        flow[i*siglength+i_2] = u[i] * Kval * v[i_2];
        viol[i] += vdiff * Kval * u[i];
      }
      viol_2[i_2] = v[i_2] * K_dot_u - screensig[i_2];
    }

    if (current_threshval <= tolerance)
      break;
    else if (iter == maxiter-1 && threadIdx.x == 0)
      printf("Warning: EMD did not converge\n");
  }

  // use flow and cost matrices to get overall cost
  float cost = 0;
  CUDA_KERNEL_LOOP(i,sigsq) {
    unsigned idx1 = i / siglength;
    unsigned idx2 = i % siglength;
    cost += cost_matrix(idx1, idx2) * flow[i];
  }
  cost = block_sum<float>(cost);
  if (threadIdx.x == 0)
    *scoregrid = cost;
}

void do_gpu_emd(const float* optgrid, const float* screengrid, float* scoregrid, 
    unsigned dim, unsigned subgrid_dim, unsigned blocks_per_side, unsigned ntypes, 
    size_t gsize, flmat& cost_matrix) {
  float reg = 9;
  // allocate memory for signatures, flow matrix, and helper variables
  float* optsig;
  float* screensig;
  float* flow;
  float* u;
  float* v;
  float* viol;
  float* viol_2;

  unsigned ncubes = blocks_per_side * blocks_per_side * blocks_per_side;
  unsigned siglength = ncubes * ntypes;
  unsigned nelems_sq = siglength * siglength;
  CUDA_CHECK_GNINA(cudaMalloc(&optsig, sizeof(float)*siglength));
  CUDA_CHECK_GNINA(cudaMalloc(&screensig, sizeof(float)*siglength));

  CUDA_CHECK_GNINA(cudaMalloc(&flow, sizeof(float)*nelems_sq));
  CUDA_CHECK_GNINA(cudaMalloc(&u, sizeof(float)*siglength));
  CUDA_CHECK_GNINA(cudaMalloc(&v, sizeof(float)*siglength));
  CUDA_CHECK_GNINA(cudaMalloc(&viol, sizeof(float)*siglength));
  CUDA_CHECK_GNINA(cudaMalloc(&viol_2, sizeof(float)*siglength));

  // need a modified version of the cost matrix for iterating
  flmat K_cpu(siglength,0);
#pragma omp parallel for
  for (size_t i=0; i<siglength; ++i) {
    K_cpu(i) = std::exp(cost_matrix(i) / -reg);
  }
  flmat_gpu gpu_cost_matrix(cost_matrix);
  flmat_gpu K(K_cpu);

  // populate and copy u and v
  std::vector<float> u_cpu(siglength, 1.f/siglength);
  std::vector<float> v_cpu(siglength, 1.f/siglength);
  CUDA_CHECK_GNINA(cudaMemcpy(u, &u_cpu[0], sizeof(float)*u_cpu.size(),
        cudaMemcpyHostToDevice));
  CUDA_CHECK_GNINA(cudaMemcpy(v, &v_cpu[0], sizeof(float)*v_cpu.size(),
        cudaMemcpyHostToDevice));

  // accumulate weights for each subcube; this is easiest if we can use 
  // block_sum, which means we need a thread block for each cube with some 
  // funny indexing
  unsigned threads_per_block = subgrid_dim * subgrid_dim * subgrid_dim;
  unsigned remainder = siglength % 1024;
  unsigned roundsize = siglength;
  if (remainder)
    roundsize += 1024 - remainder;
  // TODO: depending on calcSig performance, could store the optgrid signature
  // since it is reused...but it might be pretty fast
  gpu_calcSig<<<siglength, threads_per_block>>>(optgrid, optsig,
      subgrid_dim, dim, gsize, ntypes, blocks_per_side);
  gpu_normalizeSig<<<1, 1024>>>(optsig, siglength, roundsize);

  gpu_calcSig<<<siglength, threads_per_block>>>(screengrid, screensig,
      subgrid_dim, dim, gsize, ntypes, blocks_per_side);
  gpu_normalizeSig<<<1, 1024>>>(screensig, siglength, roundsize);

  gpu_emd<<<1, 1024>>>(optsig, screensig, scoregrid, dim, subgrid_dim, 
      ntypes, gsize, gpu_cost_matrix, siglength, flow, K, u, v, viol,
      viol_2);
  KERNEL_CHECK;
  // TODO: dump the flow matrix someplace
}

__global__
void constant_fill(float* fillgrid, size_t gsize, float fillval) {
  unsigned tidx = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned nthreads = blockDim.x * gridDim.x;
  for (size_t i=tidx; i<gsize; i+=nthreads) {
    if (fillgrid[i] == 0)
      fillgrid[i] = fillval;
  }
}

void do_constant_fill(float* fillgrid, size_t gsize, float fillval) {
  // sets any 0 values to constant fill value
  unsigned block_multiple = gsize / CUDA_NUM_THREADS;
  unsigned nblocks = block_multiple < CUDA_NUM_BLOCKS ? block_multiple : CUDA_NUM_BLOCKS;
  constant_fill<<<nblocks, CUDA_NUM_THREADS>>>(fillgrid, gsize, fillval);
  KERNEL_CHECK;
}
