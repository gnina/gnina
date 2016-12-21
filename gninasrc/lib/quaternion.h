/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_QUATERNION_H
#define VINA_QUATERNION_H

#include <boost/math/quaternion.hpp>
#include <boost/serialization/split_free.hpp>

#include "common.h"
#include "random.h"

typedef boost::math::quaternion<fl> bqt;

struct qt : float4{
    __host__ __device__ qt(){};
    __host__ __device__
    qt(fl _x, fl _y, fl _z, fl _w){
        x = _x;
        y = _y;
        z = _z;
        w = _w;
    };

    /* TODO: get rid of this */
    __host__ __device__
    fl R_component_1() const{return x;}
    __host__ __device__
    fl R_component_2() const{return y;}
    __host__ __device__
    fl R_component_3() const{return z;}
    __host__ __device__
    fl R_component_4() const{return w;}

    __host__ __device__
    qt &operator *=(const fl &r){
        x *= r;
        y *= r;
        z *= r;
        w *= r;
        return *this;
    }
    __host__ __device__
    qt operator /=(const fl &r){
        x /= r;
        y /= r;
        z /= r;
        w /= r;
        return *this;
    }
    __host__ __device__
    qt operator /=(const qt &r);


    operator bqt&(){
        return *(bqt *)this;
    }

    __host__ __device__
    qt operator*(qt r) const;
};

// non-intrusive free function split serialization
namespace boost {
	namespace serialization {
		template<class Archive>
		void save(Archive& ar, const bqt& q, const unsigned version) {
			fl q1 = q.R_component_1();
			fl q2 = q.R_component_2();
			fl q3 = q.R_component_3();
			fl q4 = q.R_component_4();

			ar & q1;
			ar & q2;
			ar & q3;
			ar & q4;
		}
		template<typename Archive>
		void load(Archive& ar, bqt& q, const unsigned version) {
			fl a, b, c, d;
			ar & a;
			ar & b;
			ar & c;
			ar & d;
			q = qt(a, b, c, d);
		}
	}
}
BOOST_SERIALIZATION_SPLIT_FREE(bqt)



bool eq(const qt& a, const qt& b); // elementwise approximate equality - may return false for equivalent rotations
const qt qt_identity(1, 0, 0, 0);
/* qt angle_to_quaternion(const vec& axis, fl angle); // axis is assumed to be a unit vector */
qt angle_to_quaternion(const vec& rotation); // rotation == angle * axis
vec quaternion_to_angle(const qt& q);
qt random_orientation(rng& generator);
void quaternion_increment(qt& q, const vec& rotation);
vec quaternion_difference(const qt& b, const qt& a); // rotation that needs to be applied to convert a to b
void print(const qt& q, std::ostream& out = std::cout); // print as an
                                                        // angle

__host__ __device__
inline fl quaternion_norm_sqr(const qt& q) { // equivalent to sqr(boost::math::abs(const qt&))
	return sqr(q.R_component_1()) + sqr(q.R_component_2()) + sqr(q.R_component_3()) + sqr(q.R_component_4());
}

inline bool quaternion_is_normalized(const qt& q) { // not in the interface, used in assertions
	return eq(quaternion_norm_sqr(q), 1) && eq(boost::math::abs((bqt&) q), 1);
}

inline void quaternion_normalize(qt& q) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs((bqt&)q))));
    const fl a = sqrt(s);
	assert(a > epsilon_fl);
	q *= 1/a;
	assert(quaternion_is_normalized(q));
}

__host__ __device__
inline qt quaternion_normalize_approx(const qt& q, const fl tolerance = 1e-6) {
	const fl s = quaternion_norm_sqr(q);
	/* assert(eq(s, sqr(abs(q)))); */
    if(abs(s - 1) < tolerance)
        return q; // most likely scenario
    else {
        const fl a = sqrt(s);
        assert(a > epsilon_fl);
        qt r = q;
        r *= 1/a;
        return r;
        /* assert(quaternion_is_normalized(q)); */
    }
}

__host__ __device__
inline
void g_normalize_angle(fl& x) { // subtract or add enough 2*pi's to make x be in [-pi, pi]

    if(x >  3*pi) { // very large
		fl n = ( x - pi) / (2*pi); // how many 2*pi's do you want to subtract?
		x -= 2*pi*ceil(n); // ceil can be very slow, but this should not be called often
		g_normalize_angle(x);
	}
	else if(x < -3*pi) { // very small
		fl n = (-x - pi) / (2*pi); // how many 2*pi's do you want to add?
		x += 2*pi*ceil(n); // ceil can be very slow, but this should not be called often
		g_normalize_angle(x);
	}
	else if(x >    pi) { // in (   pi, 3*pi]
		x -= 2*pi;
	}
	else if(x <   -pi) { // in [-3*pi,  -pi)
		x += 2*pi;
	}
	assert(x >= -pi && x <= pi);
	// in [-pi, pi]
}


__host__ __device__
inline
qt angle_to_quaternion(const vec& axis, fl angle) { // axis is assumed to be a unit vector
	//assert(eq(tvmet::norm2(axis), 1));
	/* assert(eq(axis.norm(), 1)); */
	g_normalize_angle(angle); // this is probably only necessary if angles can be very big
	fl c = cos(angle/2);
	fl s = sin(angle/2);
	return qt(c, s*axis[0], s*axis[1], s*axis[2]);
}

__host__ __device__
inline
qt qt::operator*(qt r) const{
    /* TODO: renaming to avoid messing with the actual expression.  */
    const fl &a = x;
    const fl &b = y;
    const fl &c = z;
    const fl &d = w;

    const fl ar = r.x;
    const fl br = r.y;
    const fl cr = r.z;
    const fl dr = r.w;

    return qt(
        +a*ar-b*br-c*cr-d*dr,
        +a*br+b*ar+c*dr-d*cr,
        +a*cr-b*dr+c*ar+d*br,
        +a*dr+b*cr-c*br+d*ar
        );
}

inline
qt qt::operator/=(const qt& r){
    fl &a = x;
    fl &b = y;
    fl &c = z;
    fl &d = w;

    const fl &ar = r.x;
    const fl &br = r.y;
    const fl &cr = r.z;
    const fl &dr = r.w;

    fl denominator = ar*ar+br*br+cr*cr+dr*dr;

    fl at = (+a*ar+b*br+c*cr+d*dr)/denominator;
    fl bt = (-a*br+b*ar-c*dr+d*cr)/denominator;
    fl ct = (-a*cr+b*dr+c*ar-d*br)/denominator;
    fl dt = (-a*dr-b*cr+c*br+d*ar)/denominator;

    a = at;
    b = bt;
    c = ct;
    d = dt;

    return(*this);
}


__host__ __device__
inline
mat quaternion_to_r3(const qt& q) {
	/* assert(quaternion_is_normalized(q)); */

	const fl a = q.R_component_1();
	const fl b = q.R_component_2();
	const fl c = q.R_component_3();
	const fl d = q.R_component_4();

	const fl aa = a*a;
	const fl ab = a*b;
	const fl ac = a*c;
	const fl ad = a*d;
	const fl bb = b*b;
	const fl bc = b*c;
	const fl bd = b*d;
	const fl cc = c*c;
	const fl cd = c*d;
	const fl dd = d*d;

	/* assert(eq(aa+bb+cc+dd, 1)); */

	mat tmp;

	// from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
	tmp(0, 0) = (aa+bb-cc-dd);
	tmp(0, 1) = 2*(-ad+bc);
	tmp(0, 2) = 2*(ac+bd);

	tmp(1, 0) = 2*(ad+bc);
	tmp(1, 1) = (aa-bb+cc-dd);
	tmp(1, 2) = 2*(-ab+cd);

	tmp(2, 0) = 2*(-ac+bd);
	tmp(2, 1) = 2*(ab+cd);
	tmp(2, 2) = (aa-bb-cc+dd);

	return tmp;
}

#endif
