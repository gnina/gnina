#include "grid_gpu.h"
#include "gpucode.h"
#include "curl.h"

/*evaluate using grid; TODO: establish code path for energy-only eval*/
__device__ void grid_gpu::evaluate(const atom_params& a, float slope, float v, force_energy_tup&
        e_and_deriv, const grid_gpu* deriv=NULL) const {
	//charge indep
	evaluate_aux(data, a, slope, v, e_and_deriv, deriv);
	if (a.charge != 0 && chargedata.dim0() > 0)
	{
		//charge dependent
		force_energy_tup cderiv(0,0,0,0);
		evaluate_aux(chargedata, a, slope, v, cderiv, NULL);
        e_and_deriv.energy += a.charge * cderiv.energy;
		e_and_deriv.minus_force += a.charge * cderiv.minus_force;
	}
}

__device__ void grid_gpu::evaluate_user(const atom_params& a, float slope, force_energy_tup& deriv) const {
    evaluate_aux(data, a, slope, (fl) 1000, deriv);
}

__device__ void grid_gpu::evaluate_aux(const array3d_gpu<fl, fl>& data, const atom_params& atom, float slope,
        float v, force_energy_tup& e_and_deriv, const grid_gpu* deriv) const {
	float3 s = (atom.coords - m_init) * m_factor;
    float f = tex3D(data.tex, s[0], s[1], s[2]);
	
	float3 miss(0, 0, 0);
	int region[3];
	sz a[3];

	VINA_FOR(i, 3)
	{
		if (s[i] < 0)
		{
			miss[i] = -s[i];
			region[i] = -1;
		}
		else if (s[i] >= m_dim_fl_minus_1[i])
		{
			miss[i] = s[i] - m_dim_fl_minus_1[i];
			region[i] = 1;
			assert(data.dim(i) >= 2);
		}
		else
		{
			region[i] = 0; // now that region is boost::array, it's not initialized
		}
	}
	const fl penalty = slope * dot(miss, m_factor_inv); 
	assert(penalty > -epsilon_fl);

    if (deriv) {
        
    }

	float3 gradient(x_g, y_g, z_g);
	curl(f, gradient, v);
	float3 gradient_everywhere;

	VINA_FOR(i, 3)
	{
		gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
		e_and_deriv.minus_force[i] += m_factor[i] * gradient_everywhere[i]
				+ slope * region[i];
	}
    e_and_deriv.energy += f + penalty;
}

