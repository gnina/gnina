//dkoes - class for holding charge dependent factors of terms
//in general, this may be used to factor apart components of a term
//calculation that are not exculsively dependent on atom type
//to support better precalculation; any instance of result components
//represents the value at a specific distance and atom type combination

#ifndef SMINA_SCORING_FUNCTION_H
#define SMINA_SCORING_FUNCTION_H

#include "atom_base.h"
#include <cstdlib>
#include <cstring>

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
	}

	result_components(const std::vector<double>& vals)
	{
		sz valn = vals.size();
		for(sz i = 0; i < Last; i++) {
			if(i < valn)
				components[i] = vals[i];
			else
				components[i] = 0;
		}
	}

	result_components(const std::vector<float>& vals)
	{
		sz valn = vals.size();
		for(sz i = 0; i < Last; i++) {
			if(i < valn)
				components[i] = vals[i];
			else
				components[i] = 0;
		}
	}

	fl eval(const atom_base& a, const atom_base& b) const
	{
		return components[TypeDependentOnly] +
					std::abs(a.charge)*components[AbsAChargeDependent] +
					std::abs(b.charge)*components[AbsBChargeDependent] +
					a.charge*b.charge*components[ABChargeDependent];
	}

	//if you know the scoring function doesn't have charge dependencies
	//this is faster
	fl eval_charge_independent() const
	{
		return components[TypeDependentOnly];
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

	//the order components may be sensitive to which atom is a and which is b
	//swapOrder to change this to represent the result component if they had
	//been called in the opposite order
	void swapOrder()
	{
		std::swap(components[AbsAChargeDependent],components[AbsBChargeDependent]);
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
