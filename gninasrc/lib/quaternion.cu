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

#include "quaternion.h"

bool eq(const qt& a, const qt& b) { // elementwise approximate equality - may return false for equivalent rotations
  return eq(a.R_component_1(), b.R_component_1())
      && eq(a.R_component_2(), b.R_component_2())
      && eq(a.R_component_3(), b.R_component_3())
      && eq(a.R_component_4(), b.R_component_4());
}

__host__  __device__ qt angle_to_quaternion(const vec& rotation) {
  //fl angle = tvmet::norm2(rotation); 
  fl angle = rotation.norm();
  if (angle > epsilon_fl) {
    //vec axis; 
    //axis = rotation / angle;	
    vec axis = (1 / angle) * rotation;
    return angle_to_quaternion(axis, angle);
  }
  //TODO: should be qt_identity, temporarily changing for device compatibility
  return qt(1, 0, 0, 0);
}

#ifndef __CUDA_ARCH__
vec quaternion_to_angle(const qt& q) {
  assert(quaternion_is_normalized(q));
  const fl c = q.R_component_1();
  if(c > -1 && c < 1) { // c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
    fl angle = 2*std::acos(c);// acos is in [0, pi]
    if(angle > pi)
    angle -= 2*pi;// now angle is in [-pi, pi]
    vec axis(q.R_component_2(), q.R_component_3(), q.R_component_4());
    fl s = std::sin(angle/2);// perhaps not very efficient to calculate sin of acos
    if(std::abs(s) < epsilon_fl)
    return zero_vec;
    axis *= (angle / s);
    return axis;
  }
  else // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
  return zero_vec;
}
#else
__device__ vec quaternion_to_angle(const qt& q) {
  assert(quaternion_is_normalized(q));
  const fl c = q.R_component_1();
  if (c > -1 && c < 1) { // c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
    fl angle = 2 * acosf(c); // acos is in [0, pi]
    if (angle > pi) angle -= 2 * pi; // now angle is in [-pi, pi]
    vec axis(q.R_component_2(), q.R_component_3(), q.R_component_4());
    fl s = sinf(angle / 2); // perhaps not very efficient to calculate sin of acos
    if (fabsf(s) < epsilon_fl) return vec(0, 0, 0);
    axis *= (angle / s);
    return axis;
  } else
    // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
    return vec(0, 0, 0);
}
#endif

qt random_orientation(rng& generator) {
  qt q(random_normal(0, 1, generator), random_normal(0, 1, generator),
      random_normal(0, 1, generator), random_normal(0, 1, generator));
  fl nrm = abs(q);
  if (nrm > epsilon_fl) {
    q /= nrm;
    assert(quaternion_is_normalized(q));
    return q;
  } else
    return random_orientation(generator); // this call should almost never happen
}

__host__ __device__ void quaternion_increment(qt& q, const vec& rotation) {
  assert(quaternion_is_normalized(q));
  q = angle_to_quaternion(rotation) * q;
  q = quaternion_normalize_approx(q); // normalization added in 1.1.2
}

vec quaternion_difference(const qt& b, const qt& a) { // rotation that needs to be applied to convert a to b
  quaternion_is_normalized(a);
  quaternion_is_normalized(b);
  qt tmp = b;
  tmp /= a; // b = tmp * a    =>   b * inv(a) = tmp 
  return quaternion_to_angle(tmp); // already assert normalization
}

void print(const qt& q, std::ostream& out) { // print as an angle
  out << "(" << q.a << "," << q.b << "," << q.c << "," << q.d << ")";
  //print(quaternion_to_angle(q), out);
}

