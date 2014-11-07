#ifdef SMINA_GPU
#ifndef __PRECALCULATE_GPU_H
#define __PRECALCULATE_GPU_H
/* host code for managing spline gpu calculations */

#include "precalculate.h"
#include "gpucode.h" //all gpu device code is broken out into standalone functions declared here
struct FloatSplineData //gpu works best with flaots
{
	float x, a, b, c, d;

	FloatSplineData()
	{
	}
	FloatSplineData(SplineData d) :
			x(d.x), a(d.a), b(d.b), c(d.c), d(d.d)
	{
	}
};

// dkoes - using cubic spline interpolation instead of linear for nice
// smooth gradients; stores spline data on gpu
class precalculate_gpu: public precalculate
{

	//evaluates splines at t1/t2 and r, properly swaping result
	component_pair evaldata(smt t1, smt t2, fl r) const
			{
		std::vector<float> vals, derivs;
		if (t1 <= t2)
		{
			evaluate_splines_host(data(t1, t2), r, vals, derivs);
			return component_pair(result_components(vals),
					result_components(derivs));
		}
		else
		{
			evaluate_splines_host(data(t2, t1), r, vals, derivs);
			component_pair ret = component_pair(result_components(vals),
					result_components(derivs));
			ret.first.swapOrder();
			ret.second.swapOrder();
			return ret;
		}
	}

	//create control points for spline
	//poitns indexed by component first; nonzero indexec by component
	void setup_points(std::vector<std::vector<pr> >& points, smt t1, smt t2,
			double cutoff, unsigned n, const scoring_function& sf) const
			{
		assert(n >= 2);
		fl fraction = cutoff / (fl) n;
		sz numc = sf.num_used_components();

		//clear out arguments
		points.resize(numc);
		for (sz i = 0; i < numc; i++)
		{
			points[i].clear();
			points[i].reserve(numc + 1);
		}

		//compute points
		for (unsigned i = 0; i < n; i++)
		{
			fl xval = i * fraction;
			result_components res = sf.eval_fast(t1, t2, xval);
			for (unsigned c = 0; c < numc; c++)
			{
				points[c].push_back(pr(xval, res[c]));
			}
		}
		//last point at cutoff is zero
		for (unsigned c = 0; c < numc; c++)
		{
			points[c].push_back(pr(cutoff, 0));
		}
	}
public:
	precalculate_gpu(const scoring_function& sf, fl factor_) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
			precalculate(sf),
					data(num_atom_types(), GPUSplineInfo()),
					deviceData(NULL),
					delta(0.000005),
					factor(factor_)
	{
		VINA_CHECK(factor > epsilon_fl);
		unsigned n = factor * m_cutoff;
		unsigned pitch = data.dim();
		unsigned datasize = pitch * (pitch + 1) / 2;
		GPUSplineInfo devicedata[datasize]; //todo, investigate using 2d array

		if(sf.has_slow())
		{
			std::cerr << "GPU acceleration does not support 'slow' scoring terms\n";
			abort();
		}
		VINA_FOR(t1, data.dim())
		{

			VINA_RANGE(t2, t1, data.dim())
			{
				//create the splines and copy the data over to the gpu
				//todo: do spline interpolation on gpu
				Spline spline;
				//create spline
				std::vector<std::vector<pr> > points;
				setup_points(points, (smt) t1, (smt) t2, m_cutoff, n, sf);
				std::vector<void *> deviceMem;
				float fraction = 0, cutoff = 0;
				float **splines = NULL;
				cudaMalloc(&splines, sizeof(float*) * points.size());
				for (sz i = 0, n = points.size(); i < n; i++)
				{
					//the size of points is the number of components
					Spline spline;
					spline.initialize(points[i]);
					fraction = spline.getFraction();
					cutoff = spline.getCutoff();
					//conver to float
					const std::vector<SplineData>& data = spline.getData();
					std::vector<FloatSplineData> fldata;
					fldata.reserve(data.size());
					for (sz j = 0, m = data.size(); j < m; j++)
					{
						fldata.push_back(FloatSplineData(data[j]));
					}
					//create cuda memory
					deviceMem.push_back(NULL);
					cudaMalloc(&deviceMem.back(),
							sizeof(FloatSplineData) * fldata.size());
					cudaMemcpy(deviceMem.back(), &fldata[0],
							sizeof(FloatSplineData) * fldata.size(),
							cudaMemcpyHostToDevice);
				}

				//copy array of pointers over to device
				cudaMemcpy(splines, &deviceMem[0],
						sizeof(float*) * deviceMem.size(),
						cudaMemcpyHostToDevice);
				GPUSplineInfo info;
				info.n = deviceMem.size();
				info.splines = splines;
				info.fraction = fraction;
				info.cutoff = cutoff;
				data(t1, t2) = info;

				unsigned tindex = triangular_matrix_index_permissive(pitch, t1,
						t2);
				devicedata[tindex] = info;
			}
		}
		//copy devicedata to gpu

		cudaMalloc(&deviceData, sizeof(GPUSplineInfo) * datasize);
		cudaMemcpy(deviceData, devicedata, sizeof(GPUSplineInfo) * datasize,
				cudaMemcpyHostToDevice);
	}

	virtual ~precalculate_gpu()
	{
		//need to cuda free splines
		if (deviceData)
			cudaFree(deviceData);

		unsigned n = data.dim();
		for (unsigned i = 0; i < n; i++)
		{
			for (unsigned j = i; j < n; j++)
			{
				//free each spline, which means we have to get the pointers
				const GPUSplineInfo& info = data(i, j);
				if (info.splines)
				{
					void *splines[info.n];
					cudaMemcpy(splines, info.splines, sizeof(float*) * info.n,
							cudaMemcpyDeviceToHost);
					for(unsigned k = 0; k < info.n; k++)
					{
						if(splines[k]) cudaFree(splines[k]);
					}
					cudaFree(info.splines);
				}
			}
		}
	}

	GPUSplineInfo *getDeviceData() const { return deviceData; }

	unsigned num_types() const { return data.dim(); }

	result_components eval_fast(smt t1, smt t2, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		fl r = sqrt(r2);
		return evaldata(t1, t2, r).first;
	}

	pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const
			{
		assert(r2 <= m_cutoff_sqr);
		smt t1 = a.get();
		smt t2 = b.get();
		fl r = sqrt(r2);

		component_pair rets = evaldata(t1, t2, r);
		pr ret(rets.first.eval(a, b), rets.second.eval(a, b));

		if (scoring.has_slow())
		{
			//compute value and numerical derivative directly from function
			fl X = scoring.eval_slow(a, b, r);
			ret.first += X;

			fl rhi = r + delta;
			fl rlo = r - delta;
			if (rlo < 0)
				rlo = 0;
			if (rhi > m_cutoff)
				rhi = m_cutoff;

			fl W = scoring.eval_slow(a, b, rlo);
			fl Y = 0;
			if (rhi < m_cutoff)
				Y = scoring.eval_slow(a, b, rhi);

			fl dx = (Y - W) / (rhi - rlo);
			ret.second += dx;
		}
		ret.second /= r;
		return ret;
	}

private:

	triangular_matrix<GPUSplineInfo> data;
	GPUSplineInfo *deviceData;
	fl delta;
	fl factor;
};

#endif
#endif
