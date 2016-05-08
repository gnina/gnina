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

typedef boost::math::quaternion<fl> qt;

// non-intrusive free function split serialization
namespace boost {
	namespace serialization {
		template<class Archive>
		void save(Archive& ar, const qt& q, const unsigned version) {
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
		void load(Archive& ar, qt& q, const unsigned version) {
			fl a, b, c, d;
			ar & a;
			ar & b;
			ar & c;
			ar & d;
			q = qt(a, b, c, d);
		}
	}
}
BOOST_SERIALIZATION_SPLIT_FREE(qt)

bool eq(const qt& a, const qt& b); // elementwise approximate equality - may return false for equivalent rotations
const qt qt_identity(1, 0, 0, 0);
qt angle_to_quaternion(const vec& axis, fl angle); // axis is assumed to be a unit vector
qt angle_to_quaternion(const vec& rotation); // rotation == angle * axis
vec quaternion_to_angle(const qt& q);
mat quaternion_to_r3(const qt& q);

inline fl quaternion_norm_sqr(const qt& q) { // equivalent to sqr(boost::math::abs(const qt&))
	return sqr(q.R_component_1()) + sqr(q.R_component_2()) + sqr(q.R_component_3()) + sqr(q.R_component_4());
}

inline bool quaternion_is_normalized(const qt& q) { // not in the interface, used in assertions
	return eq(quaternion_norm_sqr(q), 1) && eq(boost::math::abs(q), 1);
}

inline void quaternion_normalize(qt& q) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs(q))));
    const fl a = std::sqrt(s);
	assert(a > epsilon_fl);
	q *= 1/a;
	assert(quaternion_is_normalized(q));
}

inline void quaternion_normalize_approx(qt& q, const fl tolerance = 1e-6) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs(q))));
    if(std::abs(s - 1) < tolerance)
        ; // most likely scenario
    else {
        const fl a = std::sqrt(s);
        assert(a > epsilon_fl);
        q *= 1/a;
        assert(quaternion_is_normalized(q));
    }
}

qt random_orientation(rng& generator);
void quaternion_increment(qt& q, const vec& rotation);
vec quaternion_difference(const qt& b, const qt& a); // rotation that needs to be applied to convert a to b
void print(const qt& q, std::ostream& out = std::cout); // print as an angle

/* TODO: remove */

__host__ __device__
inline
void g_normalize_angle(fl& x) { // subtract or add enough 2*pi's to
     // make x be in [-pi, pi]
    assert(x < 3 * pi);
    assert(x > -3 * pi);
	if(x >    pi) { // in (   pi, 3*pi]
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
qt g_angle_to_quaternion(const vec& axis, fl angle) { // axis is assumed to be a unit vector
	//assert(eq(tvmet::norm2(axis), 1));
	/* assert(eq(axis.norm(), 1)); */
	g_normalize_angle(angle); // this is probably only necessary if angles can be very big
	fl c = cos(angle/2);
	fl s = sin(angle/2);
    float4 ret = make_float4(c, s*axis[0], s*axis[1], s*axis[2]);
	return *(qt *)&ret;
}

__host__ __device__
inline
fl g_quaternion_norm_sqr(const qt& q) { // equivalent to
                                        // sqr(boost::math::abs(const
                                        // qt&))
    float4 punned = *(float4 *)&q;
	return sqr(punned.x) + sqr(punned.y) + sqr(punned.z) + sqr(punned.w);
}

__host__ __device__
inline
qt& g_scale(qt &l, fl r){
    float4 q = *(float4 *)&l;

    q.x *= r;
    q.y *= r;
    q.z *= r;
    q.w *= r;

    *(float4 *)&l = q;

    return(l);
}

__host__ __device__
inline
qt g_quaternion_normalize_approx(qt q, const fl tolerance = 1e-6) {
	const fl s = g_quaternion_norm_sqr(q);
	/* assert(eq(s, sqr(boost::math::abs(q)))); */
    if(abs(s - 1) < tolerance)
        ; // most likely scenario
    else {
        const fl a = sqrt(s);
        /* assert(a > epsilon_fl); */
        g_scale(q, 1/a);
        /* assert(quaternion_is_normalized(q)); */
    }
    return q;
}

__host__ __device__
inline
void g_quaternion_write(qt *to, const qt& q) {
    *(float4 *)to = *(float4 *)&q;
}

__host__ __device__
mat g_quaternion_to_r3(const qt& q) {
	/* assert(quaternion_is_normalized(q)); */

    float4 punned = *(float4 *)&q;

	const fl a = punned.x;
	const fl b = punned.y;
	const fl c = punned.z;
	const fl d = punned.w;

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
