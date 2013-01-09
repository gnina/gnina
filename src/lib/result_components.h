//dkoes - class for holding charge dependent factors of terms
//in general, this may be used to factor apart components of a term
//calculation that are not exculsively dependent on atom type
//to support better precalculation; any instance of result components
//represents the value at a specific distance and atom type combination
#include "atom_base.h"

//a charge dependent term must be separated into components
//that are dependent on different charges to enable precalculation
class result_components
{
	enum
	{
		TypeDependentOnly,//no need to adjust by charge
		AbsAChargeDependent,//multiply by the absolute value of a's charge
		AbsBChargeDependent,//multiply by abs(b)'s charge
		ABChargeDependent,//multiply by a*b
		Last
	};
	fl components[Last];

public:
	//add in rhs
	result_components& operator+=(const result_components& rhs)
	{
		for(unsigned i = 0; i < Last; i++)
			components[i] += rhs.components[i];

		return *this;
	}

	//add in a non-charge dependent term
	result_components& operator+=(fl val)
	{
		components[TypeDependentOnly] += val;
		return *this;
	}

	result_components()
	{
		memset(components, 0, sizeof(components));
	}

	result_components(const flv& vals)
	{
		assert(vals.size() == Last);
		for(sz i = 0; i < Last; i++)
			components[i] = vals[i];
	}

	fl eval(const atom_base& a, const atom_base& b) const
	{
		return components[TypeDependentOnly] +
					std::abs(a.charge)*components[AbsAChargeDependent] +
					std::abs(b.charge)*components[AbsBChargeDependent] +
					a.charge*b.charge*components[ABChargeDependent];
	}

	static sz size() { return Last; }
	fl operator[](sz i) const
	{
		assert(i < Last);
		return components[i];
	}

	fl& operator[](sz i)
	{
		assert(i < Last);
		return components[i];
	}
};
