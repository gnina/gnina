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

template<class F>
class Spline
{
	struct SplineData
	{
		fl x, a, b, c, d;
	};

	//vector of calculated spline data
	std::vector<SplineData> data;
	fl cutoff;
	fl fraction;
	unsigned numpoints;

	//Function to calculate the value and derivativeof a given spline at a point xval
	inline pr splineCalc(const SplineData& i, double xval) const
	{
		const fl lx = xval - i.x;
		fl val = ((i.a * lx + i.b) * lx + i.c) * lx + i.d;
		fl dx = (3*i.a*lx + 2*i.b) * lx + i.c;
		return pr(val,dx);
	}

public:

	Spline(): cutoff(0), fraction(0), numpoints(0) {}

	//A spline is initialized by specifying a function object to interpolate,
	//a cutoff, the number of control points to generate, and the number of
	//end control points to ignore when smoothing to zero
	void initialize(const F& func, double c, unsigned npoints,
			unsigned smoothed_endpoints = 0)
	{
		cutoff = c;
		numpoints = npoints;
		assert(numpoints >= 2);
		smoothed_endpoints = std::min(smoothed_endpoints, npoints - 2);
		//first generate points
		fraction = cutoff / (fl) numpoints;
		std::vector<pr> points;
		points.reserve(numpoints);
		unsigned numinterp = numpoints - smoothed_endpoints;
		for (unsigned i = 0; i < numpoints - smoothed_endpoints; i++)
		{
			fl xval = i * fraction;
			fl yval = func(xval);
			points.push_back(pr(xval, yval));
		}
		points.push_back(pr(cutoff, 0));

		const unsigned e = points.size() - 1;

		typedef Eigen::Matrix<fl, Eigen::Dynamic, Eigen::Dynamic> flMatrix;
		flMatrix A = flMatrix::Zero(points.size(), points.size());
		fl hlast = points[e].first - points[e - 1].first;
		for (unsigned i = 1; i < e; ++i)
		{
			fl hi = fraction;
			//last point may not have fixed delta due to smoothing
			if (i == e - 1)
				hi = hlast;
			A(i - 1, i) = hi;
			A(i, i) = 2 * (fraction + hi);
			A(i + 1, i) = hi;
		}

		typedef Eigen::Matrix<fl, 1, Eigen::Dynamic> flVector;
		flVector C = flVector::Zero(points.size());

		for (unsigned i = 1; i < e; ++i)
		{
			fl hi = fraction;
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
			fl hi = fraction;
			if (i == e - 1)
				hi = hlast;
			data[i].x = points[i].first;
			data[i].a = (ddy(i + 1) - ddy(i)) / (6 * hi);
			data[i].b = ddy(i) / 2;
			data[i].c = (points[i + 1].second - points[i].second) / hi - ddy(i + 1) * hi / 6
					- ddy(i) * hi / 3;
			data[i].d = points[i].second;
		}

		//if there were smoothed out points, add them here for easy lookup
		for (unsigned i = 0; i < smoothed_endpoints; i++)
		{
			data.push_back(data[e - 1]);
		}

		assert(data.size() == numpoints);

	}

	//return both the value and derivative at xval
	pr operator()(fl xval) const
	{
		//Special cases when we're outside the range of the spline points
		//distances should never be less than zero
		assert(xval >= 0);

		if (xval >= cutoff)
		{
			return pr(0, 0);
		}

		unsigned index = xval / fraction; //xval*numpoints/cutoff
		return splineCalc(data[index], xval);
	}

};

