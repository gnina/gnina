/*

   GPU optimized versions for conf and change.

*/

#include "conf_gpu.h"

#define GNINA_CUDA_NUM_THREADS (512)
#define WARPSIZE (32)

#define CUDA_KERNEL_LOOP(i, n) \
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; \
       i < (n); \
       i += blockDim.x * gridDim.x)

__global__ void scalar_mult_kernel(float mult, const int n,
		float *vals) {
	CUDA_KERNEL_LOOP(index, n)
	{
		vals[index] *= mult;
	}
}

//compute a -= b (result goes in a
__global__ void vec_sub_kernel(const int n,
		float *a, float *b) {
	CUDA_KERNEL_LOOP(index, n)
	{
		a[index] -= b[index];
	}
}

//compute dot(a,b) using a single warp -> call with exactly WARPSIZE threads
__global__ void warp_dot_kernel(const int n, float *a, float *b, float *out)
{
	int start = blockIdx.x * blockDim.x + threadIdx.x;
	float val = 0.0;
	int warpSize = blockDim.x; //allow for dynamic warp sizes (in practice 16 or 32)
	for(int i = start; i < n; i += warpSize)
	{
		val += a[i]*b[i];
	}
	//now warp reduce with shuffle

	for(uint offset = warpSize/2; offset > 0; offset >>= 1)
		val += __shfl_down(val, offset);

	if(start == 0)
		*out = val;
}

__global__ void minus_mat_vec_product_kernel(const int n, flmat_gpu m, float*
        in, float* out) {
    VINA_FOR(i,n) {
        fl sum = 0;
        VINA_FOR(j,n)
            sum += m(m.index_permissive(i,j)) * in[j];
        out[i] = -sum;
    }
}

change_gpu::change_gpu(const change& src) :
		change_values(NULL), n(0) {
	std::vector<float> data;
	//figure out number of torsions
	assert(src.ligands.size() == 1);
	n = 6; //position + orientation
	const ligand_change& lig = src.ligands[0];

	for (unsigned i = 0; i < 3; i++)
		data.push_back(lig.rigid.position[i]);

	for (unsigned i = 0; i < 3; i++)
		data.push_back(lig.rigid.orientation[i]);

	n += lig.torsions.size();
	for (unsigned i = 0, nn = lig.torsions.size(); i < nn; i++) {
		data.push_back(lig.torsions[i]);
	}

	for (unsigned i = 0, nn = src.flex.size(); i < nn; i++) {
		n += src.flex[i].torsions.size();
		for (unsigned j = 0, m = src.flex[i].torsions.size(); j < m; j++) {
			data.push_back(src.flex[i].torsions[j]);
		}
	}
	//allocate vector
	CUDA_CHECK_GNINA(cudaMalloc(&change_values, sizeof(float) * (n+1))); //leave scratch space for dot
	//and init
	assert(n == data.size());
	CUDA_CHECK_GNINA(
			cudaMemcpy(change_values, &data[0], n * sizeof(float),
					cudaMemcpyHostToDevice));
}

//allocate and copy
change_gpu::change_gpu(const change_gpu& src) :
		n(0), change_values(NULL) {
	*this = src;
}

change_gpu& change_gpu::operator=(const change_gpu& src) {
	if (change_values == NULL || n < src.n) {
		if (change_values) {
			CUDA_CHECK_GNINA(cudaFree(change_values));
		}
		CUDA_CHECK_GNINA(cudaMalloc(&change_values, sizeof(float) * (src.n+1))); //scratch space
	}
	n = src.n;
	CUDA_CHECK_GNINA(
			cudaMemcpy(change_values, src.change_values, sizeof(float) * n,
					cudaMemcpyDeviceToDevice));
	return *this;
}

change_gpu::~change_gpu() {
	//deallocate mem
	CUDA_CHECK_GNINA(cudaFree(change_values));
}

//dkoes - zeros out all differences
void change_gpu::clear() {
	CUDA_CHECK_GNINA(cudaMemset(change_values, 0, sizeof(float) * n));
}

//dkoes - multiply by -1
void change_gpu::invert() {
	scalar_mult_kernel<<<1, min(GNINA_CUDA_NUM_THREADS, n)>>>(-1.0, n,
			change_values);
}

//return dot product
float change_gpu::dot(const change_gpu& rhs) const {
	//since N is small, I think we should do a single warp of threads for this
	warp_dot_kernel<<<1, (n <= WARPSIZE/2 ? WARPSIZE/2 : WARPSIZE)>>>(n, change_values, rhs.change_values, &change_values[n]);
	float gpuval = 0;
	cudaMemcpy(&gpuval, &change_values[n], sizeof(float),cudaMemcpyDeviceToHost);
	return gpuval;
}

//subtract rhs from this
void change_gpu::sub(const change_gpu& rhs) {

	vec_sub_kernel<<<1, min(GNINA_CUDA_NUM_THREADS, n)>>>(n,
				change_values, rhs.change_values);
}

void change_gpu::minus_mat_vec_product(const flmat_gpu& m, change_gpu& out) const {
    minus_mat_vec_product_kernel<<<1,1>>>(n, m, change_values,
            out.change_values);
}

sz change_gpu::num_floats() const {
	return n;
}

//for debugging
void change_gpu::get_data(std::vector<float>& d) const {
	d.resize(n);
	CUDA_CHECK_GNINA(
			cudaMemcpy(&d[0], change_values, n * sizeof(float),
					cudaMemcpyDeviceToHost));

}

void change_gpu::set_data(std::vector<float>& d) const {
	CUDA_CHECK_GNINA(
			cudaMemcpy(change_values, &d[0], n * sizeof(float),
					cudaMemcpyHostToDevice));
}

void change_gpu::print() const {
	std::vector<float> d;
	get_data(d);
	for (unsigned i = 0, n = d.size(); i < n; i++) {
		std::cout << d[i] << " ";
	}
	std::cout << "\n";
}

conf_gpu::conf_gpu(const conf& src) :
		cinfo(NULL), n(0) {
	std::vector<float> data;
	//figure out number of torsions
	assert(src.ligands.size() == 1);
	n = 7; //position + orientation(qt)
	const ligand_conf& lig = src.ligands[0];

	for (unsigned i = 0; i < 3; i++)
		data.push_back(lig.rigid.position[i]);

	data.push_back(lig.rigid.orientation.R_component_1());
	data.push_back(lig.rigid.orientation.R_component_2());
	data.push_back(lig.rigid.orientation.R_component_3());
	data.push_back(lig.rigid.orientation.R_component_4());

	n += lig.torsions.size();
	for (unsigned i = 0, nn = lig.torsions.size(); i < nn; i++) {
		data.push_back(lig.torsions[i]);
	}

	for (unsigned i = 0, nn = src.flex.size(); i < nn; i++) {
		n += src.flex[i].torsions.size();
		for (unsigned j = 0, m = src.flex[i].torsions.size(); j < m; j++) {
			data.push_back(src.flex[i].torsions[j]);
		}
	}

	//allocate vector
	CUDA_CHECK_GNINA(cudaMalloc(&cinfo, sizeof(float) * n));
	//and init
	assert(n == data.size());
	CUDA_CHECK_GNINA(
			cudaMemcpy(cinfo, &data[0], n * sizeof(float),
					cudaMemcpyHostToDevice));
}

//set cpu to gpu values, assumes correctly sized
void conf_gpu::set_cpu(conf& dst) const {
	std::vector<float> d;
	get_data(d);
	assert(dst.ligands.size() == 1);
	unsigned pos = 0;
	if (d.size() >= 7) {
		ligand_conf& lig = dst.ligands[0];
		lig.rigid.position = vec(d[0], d[1], d[2]);
		lig.rigid.orientation = qt(d[3], d[4], d[5], d[6]);
		pos = 7;
		for (unsigned i = 0, nt = lig.torsions.size(); i < nt && pos < n;
				i++) {
			lig.torsions[i] = d[pos];
			pos++;
		}
	}

	for (unsigned r = 0, nr = dst.flex.size(); r < nr; r++) {
		residue_conf& res = dst.flex[r];
		for (unsigned i = 0, nt = res.torsions.size(); i < nt && pos < n;
				i++) {
			res.torsions[i] = d[pos];
			pos++;
		}
	}
}

//allocate and copy
conf_gpu::conf_gpu(const conf_gpu& src) :
		n(0), cinfo(NULL) {
	*this = src;
}

conf_gpu& conf_gpu::operator=(const conf_gpu& src) {
	if (cinfo == NULL || n < src.n) {
		if (cinfo) {
			CUDA_CHECK_GNINA(cudaFree(cinfo));
		}
		CUDA_CHECK_GNINA(cudaMalloc(&cinfo, sizeof(float) * src.n));
	}
	n = src.n;
	CUDA_CHECK_GNINA(
			cudaMemcpy(cinfo, src.cinfo, sizeof(float) * n,
					cudaMemcpyDeviceToDevice));
	return *this;
}

conf_gpu::~conf_gpu() {
	//deallocate mem
	CUDA_CHECK_GNINA(cudaFree(cinfo));
}

__global__ void increment_kernel(float* x, float* c, fl factor, int n) {
	//position
	for(unsigned i = 0; i < 3; i++) {
		x[i] += c[i]*factor;
	}

	//rotation
	qt orientation(x[3],x[4],x[5],x[6]);
	vec rotation(factor * c[3], factor * c[4], factor *
            c[5]);
	quaternion_increment(orientation, rotation);
	x[3] = orientation.R_component_1();
	x[4] = orientation.R_component_2();
	x[5] = orientation.R_component_3();
	x[6] = orientation.R_component_4();

	//torsions
	for(unsigned i = 7; i < n; i++) {
		x[i] += normalized_angle(factor*c[i-1]);
		normalize_angle(x[i]);
	}
}

void conf_gpu::increment(const change_gpu& c, fl factor) {
    increment_kernel<<<1,1>>>(cinfo->values, c.change_values, factor, n);
}

//for debugging (mostly)
void conf_gpu::get_data(std::vector<float>& d) const {
	d.resize(n);
	CUDA_CHECK_GNINA(
			cudaMemcpy(&d[0], cinfo, n * sizeof(float),
					cudaMemcpyDeviceToHost));

}

void conf_gpu::set_data(std::vector<float>& d) const {
	CUDA_CHECK_GNINA(
			cudaMemcpy(cinfo, &d[0], n * sizeof(float),
					cudaMemcpyHostToDevice));
}

void conf_gpu::print() const {
	std::vector<float> d;
	get_data(d);
	for (unsigned i = 0, n = d.size(); i < n; i++) {
		std::cout << d[i] << " ";
	}
	std::cout << "\n";
}
