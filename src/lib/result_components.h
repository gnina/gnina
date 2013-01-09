//dkoes - class for holding charge dependent factors of terms
//in general, this may be used to factor apart components of a term
//calculation that are not exculsively dependent on atom type
//to support better precalculation; any instance of result components
//represents the value at a specific distance and atom type combination

#ifndef SMINA_SCORING_FUNCTION_H
#define SMINA_SCORING_FUNCTION_H

#include "atom_base.h"

//a charge dependent term must be separated into components
//that are dependent on different charges to enable precalculation
class result_components
{
public:
	enum
	{
		TypeDependentOnly,//no need to adjust by charge
		AbsAChargeDependent,//multiply by the absolute value of a's charge
		AbsBChargeDependent,//multiply by abs(b)'s charge
		ABChargeDependent,//multiply by a*b
		Last
	};

	//add in rhs
	result_components& operator+=(const result_components& rhs)
	{
		for(unsigned i = 0; i < Last; i++)
			components[i] += rhs.components[i];

		return *this;
	}

	//scale ALL components
	result_components& operator*=(fl val)
	{
		for(unsigned i = 0; i < Last; i++)
			components[i] *= val;

		return *this;
	}

	result_components()
	{
		memset(components, 0, sizeof(components));
		std::cout << sizeof(components) << "\n";
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

	friend result_components operator*(const result_components& lhs, fl rhs);
private:
	fl components[Last];
};

inline result_components operator*(const result_components& lhs, fl rhs)
{
	result_components ret(lhs);
	ret *= rhs;
	return ret;
}

#endif
