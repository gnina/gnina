#include "bfgs.h"
#include <cuda_runtime.h>


__global__ void lambdamin_kernel(fl* test, fl* pvalues, fl* xvalues, sz n) {
	for (sz i = 0; i < n; i++)
	{
		fl temp = fabsf(pvalues[i]) / fmaxf(fabsf(xvalues[i]), 1.0f);
		if (temp > *test)
			*test = temp;
	}
}

fl compute_lambdamin(const change_gpu& p, const conf_gpu& x, sz n)
{
	fl* test;
    fl out_test;
    CUDA_CHECK_GNINA(cudaMalloc(&test, sizeof(fl)));
    lambdamin_kernel<<<1,1>>>(test, p.change_values, x.cinfo->values, n);
    CUDA_CHECK_GNINA(cudaMemcpy(&out_test, test, sizeof(fl), cudaMemcpyDeviceToHost));
    CUDA_CHECK_GNINA(cudaFree(test));
	return out_test;
}

__global__ void set_diagonal_kernel(flmat_gpu m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

void set_diagonal(const flmat_gpu& m, fl x)
{
    set_diagonal_kernel<<<1,1>>>(m,x);
}

__global__ void bfgs_update_aux(const sz n, const fl alpha, const fl r, float coef,
        flmat_gpu h, float* pvec, float* minus_hyvec) {
    int i = threadIdx.x;
    int j = blockIdx.x;
    if (j >= i) {
        h(i, j) += alpha * r
        		* (minus_hyvec[i] * pvec[j] + minus_hyvec[j] * pvec[i])
        		+coef * pvec[i]	* pvec[j]; // s * s == alpha * alpha * p * p	}
    }
}

void bfgs_update(const flmat_gpu& h, const change_gpu& p, const
        change_gpu& y, const fl alpha) {
	const fl yp = y.dot(p);
	if (alpha * yp < epsilon_fl)
		return; // FIXME?

	change_gpu minus_hy(y);
	y.minus_mat_vec_product(h, minus_hy);

	const fl yhy = -y.dot(minus_hy);
	const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
	const sz n = p.num_floats();

	float coef = +alpha * alpha * (r * r * yhy + r) ;
    bfgs_update_aux<<<n,n>>>(n, alpha, r, coef, h, p.change_values, minus_hy.change_values);
	return;
}

