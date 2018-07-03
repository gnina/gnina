#ifndef __PRECALCULATE_GPU_H
#define __PRECALCULATE_GPU_H
/* host code for managing spline gpu calculations */

#include "precalculate.h"
#include "gpucode.h" //all gpu device code is broken out into standalone functions declared here
struct FloatSplineData //gpu works best with flaots
{
    float x, a, b, c, d;

    FloatSplineData() {
    }
    FloatSplineData(SplineData d)
        : x(d.x), a(d.a), b(d.b), c(d.c), d(d.d) {
    }
};

// dkoes - using cubic spline interpolation instead of linear for nice
// smooth gradients; stores spline data on gpu
class precalculate_gpu : public precalculate {

    //evaluates splines at t1/t2 and r, properly swaping result
    component_pair evaldata(smt t1, smt t2, fl r) const;

    void evaluate_splines(const GPUSplineInfo& spInfo, float r,
        std::vector<float>& vals, std::vector<float>& derivs) const;

    //create control points for spline
    //poitns indexed by component first; nonzero indexec by component
    void setup_points(std::vector<std::vector<pr> >& points, smt t1, smt t2,
        double cutoff, unsigned n, const scoring_function& sf) const;
  public:
    precalculate_gpu(const scoring_function& sf, fl factor_);

    virtual ~precalculate_gpu();

    GPUSplineInfo *getDeviceData() const {
      return deviceData;
    }
    unsigned num_types() const {
      return data.dim();
    }
    result_components eval_fast(smt t1, smt t2, fl r2) const;
    pr eval_deriv(const atom_base& a, const atom_base& b, fl r2) const;
    //TODO: reprivate
    GPUSplineInfo *deviceData;

  private:

    triangular_matrix<GPUSplineInfo> data;
    triangular_matrix<spline_cache> cpudata;
    float *device_vals, *device_derivs;
    fl delta;
    fl factor;
};

#endif
