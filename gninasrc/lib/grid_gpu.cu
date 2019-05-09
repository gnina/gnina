#include "grid_gpu.h"
#include "gpucode.h"
#include "curl.h"

/*evaluate using grid; TODO: establish code path for energy-only eval*/
__device__ void grid_gpu::evaluate(const atom_params& a, float slope, float v,
    force_energy_tup& e_and_deriv) const {
  //charge indep
  evaluate_aux(data, a, slope, v, e_and_deriv);
  if (a.charge != 0 && chargedata.dim0() > 0) {
    //charge dependent
    force_energy_tup cderiv(0, 0, 0, 0);
    evaluate_aux(chargedata, a, slope, v, cderiv);
    e_and_deriv.energy += a.charge * cderiv.energy;
    e_and_deriv.minus_force += a.charge * cderiv.minus_force;
  }
}

__device__ void grid_gpu::evaluate_user(const atom_params& a, float slope,
    force_energy_tup& deriv) const {
  evaluate_aux(data, a, slope, (fl) 1000, deriv);
}

__device__ void grid_gpu::evaluate_aux(const array3d_gpu<fl, fl>& data,
    const atom_params& atom, float slope, float v,
    force_energy_tup& e_and_deriv) const {
  gfloat3 s = (atom.coords - m_init) * m_factor;

  gfloat3 miss(0, 0, 0);
  int region[3];
  sz a[3];

  VINA_FOR(i, 3) {
    if (s[i] < 0) {
      miss[i] = -s[i];
      region[i] = -1;
      a[i] = 0;
      s[i] = 0;
    } else
      if (s[i] >= m_dim_fl_minus_1[i]) {
        miss[i] = s[i] - m_dim_fl_minus_1[i];
        region[i] = 1;
        assert(data.dim(i) >= 2);
        a[i] = data.dim(i) - 2;
        s[i] = 1;
      } else {
        region[i] = 0; // now that region is boost::array, it's not initialized
        a[i] = sz(s[i]);
        s[i] -= a[i];
      }
    assert(s[i] >= 0);
    assert(s[i] <= 1);
    assert(
        a[i] + 1 < data.dim(i)
            || (printf("a[i]+1 is %d and data.dim(i) is %d\n", a[i] + 1,
                data.dim(i)) && 0));
  }
  const fl penalty = slope * dot(miss, m_factor_inv); // FIXME check that inv_factor is correctly initialized and serialized
  assert(penalty > -epsilon_fl);

  const sz x0 = a[0];
  const sz y0 = a[1];
  const sz z0 = a[2];

  const sz x1 = x0 + 1;
  const sz y1 = y0 + 1;
  const sz z1 = z0 + 1;

  const fl f000 = data(x0, y0, z0);
  const fl f100 = data(x1, y0, z0);
  const fl f010 = data(x0, y1, z0);
  const fl f110 = data(x1, y1, z0);
  const fl f001 = data(x0, y0, z1);
  const fl f101 = data(x1, y0, z1);
  const fl f011 = data(x0, y1, z1);
  const fl f111 = data(x1, y1, z1);

  const fl x = s[0];
  const fl y = s[1];
  const fl z = s[2];

  const fl mx = 1 - x;
  const fl my = 1 - y;
  const fl mz = 1 - z;

  fl f = f000 * mx * my * mz + f100 * x * my * mz + f010 * mx * y * mz
      + f110 * x * y * mz + f001 * mx * my * z + f101 * x * my * z
      + f011 * mx * y * z + f111 * x * y * z;

  /*TODO: deriv only*/
  const fl x_g = f000 * (-1) * my * mz + f100 * 1 * my * mz
      + f010 * (-1) * y * mz + f110 * 1 * y * mz + f001 * (-1) * my * z
      + f101 * 1 * my * z + f011 * (-1) * y * z + f111 * 1 * y * z;

  const fl y_g = f000 * mx * (-1) * mz + f100 * x * (-1) * mz
      + f010 * mx * 1 * mz + f110 * x * 1 * mz + f001 * mx * (-1) * z
      + f101 * x * (-1) * z + f011 * mx * 1 * z + f111 * x * 1 * z;

  const fl z_g = f000 * mx * my * (-1) + f100 * x * my * (-1)
      + f010 * mx * y * (-1) + f110 * x * y * (-1) + f001 * mx * my * 1
      + f101 * x * my * 1 + f011 * mx * y * 1 + f111 * x * y * 1;

  gfloat3 gradient(x_g, y_g, z_g);
  curl(f, gradient, v);
  gfloat3 gradient_everywhere;

  VINA_FOR(i, 3) {
    gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
    e_and_deriv.minus_force[i] += m_factor[i] * gradient_everywhere[i]
        + slope * region[i];
  }
  e_and_deriv.energy += f + penalty;
}

