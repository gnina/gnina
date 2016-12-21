/* dkoes - An implementation of cubmic splines.
 * Assume and enforce that x values are evenly spaced from
 * zero to some cutoff (user may specify a larger last step to cutoff to
 * enhance smoothing), derivatives must go to zero at ends, value goes to
 * zero at cutoff.  Eigen is used for matrix calculations.
 *
 * The spline is initialized with a function object.
 *
 * This is loosely inspired from the Spline class of dynamo:
 * https://github.com/toastedcrumpets/DynamO/blob/master/src/magnet/magnet/math/spline.hpp
 */

#include "common.h"
#include "Eigen/Core"
#include "Eigen/Dense"

typedef fl fltype;
struct SplineData
{
	fltype x, a, b, c, d;
};

class Spline
{
	//vector of calculated spline data
	std::vector<SplineData> data;
	fltype cutoff;
	fltype fraction;

public:

	Spline(): cutoff(0), fraction(0){}

	//A spline is initialized by specifying the points to interpolate,
	//the cutoff is inferred from the last point, the fraction between the
	//evenly spaced points is inferred from the first two points
	void initialize(const std::vector<pr>& points)
	{
		if(points.size() == 0) return;
		cutoff = points.back().first;
		assert(points.size() >= 2);
		fraction = points[1].first - points[0].first;
		const unsigned e = points.size() - 1;

		typedef Eigen::Matrix<fltype, Eigen::Dynamic, Eigen::Dynamic> flMatrix;
		flMatrix A = flMatrix::Zero(points.size(), points.size());
		fltype hlast = points[e].first - points[e - 1].first;
		for (unsigned i = 1; i < e; ++i)
		{
			fltype hi = fraction;
			//last point may not have fixed delta due to smoothing
			if (i == e - 1)
				hi = hlast;
			A(i - 1, i) = hi;
			A(i, i) = 2 * (fraction + hi);
			A(i + 1, i) = hi;
		}

		typedef Eigen::Matrix<fltype, 1, Eigen::Dynamic> flVector;
		flVector C = flVector::Zero(points.size());

		for (unsigned i = 1; i < e; ++i)
		{
			fltype hi = fraction;
			if (i == e - 1)
				hi = hlast;

			C(i) = 6
					*
					((points[i + 1].second - points[i].second) / hi
							- (points[i].second - points[i - 1].second)
									/ fraction);
		}

		//Boundary condition: zero first derivative
		C(0) = 6 * ((points[1].second - points[0].second) / fraction);
		A(0, 0) = 2 * fraction;
		A(1, 0) = fraction;

		C(e) = 6 * (-(points[e].second - points[e - 1].second) / hlast);
		A(e, e) = 2 * hlast;
		A(e - 1, e) = hlast;

		flMatrix AInv = A.inverse();

		flVector ddy = C * AInv;

		data.resize(e);
		for (unsigned i = 0; i < e; ++i)
		{
			fltype hi = fraction;
			if (i == e - 1)
				hi = hlast;
			data[i].x = points[i].first;
			data[i].a = (ddy(i + 1) - ddy(i)) / (6 * hi);
			data[i].b = ddy(i) / 2;
			data[i].c = (points[i + 1].second - points[i].second) / hi - ddy(i + 1) * hi / 6
					- ddy(i) * hi / 3;
			data[i].d = points[i].second;
		}
	}

	//return both the value and derivative at xval
	pr eval_deriv(fl xval) const
	{
		//Special cases when we're outside the range of the spline points
		//distances should never be less than zero
		assert(xval >= 0);

		if (xval >= cutoff)
		{
			return pr(0, 0); //so an uninitialized spline returns zero
		}

		unsigned index = xval / fraction; //xval*numpoints/cutoff
		const SplineData& i = data[index];
		const fl lx = xval - i.x;
		fl val = ((i.a * lx + i.b) * lx + i.c) * lx + i.d;
		fl dx = (3*i.a*lx + 2*i.b) * lx + i.c;
		return pr(val,dx);
	}

	//return only the value
	fl eval(fl xval) const
	{
		assert(xval >= 0);

		if (xval >= cutoff)
		{
			return 0; //so an uninitialized spline returns zero
		}

		unsigned index = xval / fraction; //xval*numpoints/cutoff
		const SplineData& i = data[index];
		const fl lx = xval - i.x;
		fl val = ((i.a * lx + i.b) * lx + i.c) * lx + i.d;
		return val;
	}

	const std::vector<SplineData>& getData() const {
		return data;
	}

	fltype getFraction() const {
		return fraction;
	}

	fltype getCutoff() const {
		return cutoff;
	}
};

