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

#ifndef VINA_COMMON_H
#define VINA_COMMON_H

#include <cassert>
#include <string>
#include <limits>
#include <utility> // pair#include <algorithm> // too common#include <vector> // used in typedef, and commonly used overall#include <cmath> // commonly used#include <iostream> // various debugging everywhere#include <fstream> // print_coords#include <iomanip> // to_string#include <sstream> // to_string#include <string> // probably included by the above anyway, common anyway#include <boost/serialization/vector.hpp> // can't come before the above two - wart fixed in upcoming Boost versions#include <boost/serialization/base_object.hpp> // movable_atom needs it - (derived from atom)#include <boost/filesystem/path.hpp> // typedef'ed#include "macros.h"
#include "math.h"
#include <cuda_runtime.h>

typedef float fl;

//collection of parameters specifying how minimization should be done
struct minimization_params {
    enum Type {
      BFGSFastLineSearch, BFGSAccurateLineSearch, ConjugateGradient, Simple
    };

    Type type;
    unsigned maxiters; //maximum number of iterations of algorithm
    bool early_term; //terminate early based on different of function values
    bool single_min; //do single full minimization instead of hunt_cap truncated followed by full
    int outputframes;
    minimization_params()
        : type(BFGSFastLineSearch), maxiters(0), early_term(false),
            single_min(false), outputframes(0) {

    }
};

template<typename T>
__host__   __device__
  inline T sqr(T x) {
  return x * x;
}

#define not_a_num 0.0
/* TODO TODO TODO TODO: reactivate. */
/* const fl not_a_num = std::sqrt(fl(-1)); // FIXME? check  */

template<class T1, class T2>
struct gpair {
    T1 first;
    T2 second;
    __host__ __device__ gpair() {
    }
    ;
    __host__ __device__ gpair(T1 f, T2 s)
        : first(f), second(s) {
    }
    ;
};

typedef std::size_t sz;
typedef unsigned short atmidx; //dkoes - to reduce size of smina format
typedef std::pair<fl, fl> pr;

// TODO: remove alignment. Exists so that vec operations in tree_gpu.cu can
// coalesce reads and writes.
#if defined(__CUDACC__)
#define CUDA_ALIGN(n) __align__(n)
#else
#define CUDA_ALIGN(n) alignas(n)
#endif

struct CUDA_ALIGN(4 * sizeof(float)) vec {
    fl data[3];
    /* TODO: remove. Exists so that force_energy_tup * can be
     interpretend as vec *. */
    fl pad[1];
    __host__ __device__ vec() {
#ifndef NDEBUG
      data[0] = data[1] = data[2] = not_a_num;
#endif
    }
    __host__ __device__ vec(fl x, fl y, fl z) {
      data[0] = x;
      data[1] = y;
      data[2] = z;
    }

    __host__ __device__ vec(fl x, fl y, fl z, fl w) {
      data[0] = x;
      data[1] = y;
      data[2] = z;
      pad[0] = w;
    }

    __host__   __device__
	  const fl& operator[](sz i) const {
      assert(i < 3);
      return data[i];
    }
    __host__   __device__ fl& operator[](sz i) {
      assert(i < 3);
      return data[i];
    }
    __host__   __device__ fl norm_sqr() const {
      /* TODO: oleg was using sqr */
      return sqr(data[0]) + sqr(data[1]) + sqr(data[2]);
    }

    __host__ __device__ fl x() const { return data[0]; }
    __host__ __device__ fl y() const { return data[1]; }
    __host__ __device__ fl z() const { return data[2]; }

#ifndef __CUDA_ARCH__
    fl norm() const {
      return std::sqrt(norm_sqr());
    }
#else
    __device__
    float norm() const {
      return sqrtf(norm_sqr());
    }
#endif
    __host__   __device__ fl operator*(const vec& v) const {
      return data[0] * v[0] + data[1] * v[1] + data[2] * v[2];
    }
    __host__   __device__
	  const vec& operator+=(const vec& v) {
      data[0] += v[0];
      data[1] += v[1];
      data[2] += v[2];
      return *this;
    }
    __host__   __device__
	  const vec& operator-=(const vec& v) {
      data[0] -= v[0];
      data[1] -= v[1];
      data[2] -= v[2];
      return *this;
    }
    __host__   __device__
	  const vec& operator+=(fl s) {
      data[0] += s;
      data[1] += s;
      data[2] += s;
      return *this;
    }
    __host__   __device__
	  const vec& operator-=(fl s) {
      data[0] -= s;
      data[1] -= s;
      data[2] -= s;
      return *this;
    }
    __host__   __device__ vec operator+(const vec& v) const {
      return vec(data[0] + v[0], data[1] + v[1], data[2] + v[2]);
    }
    __host__   __device__ vec operator-(const vec& v) const {
      return vec(data[0] - v[0], data[1] - v[1], data[2] - v[2]);
    }
    __host__   __device__
	  const vec& operator*=(fl s) {
      data[0] *= s;
      data[1] *= s;
      data[2] *= s;
      return *this;
    }

    const vec& operator/=(fl s) {
      data[0] /= s;
      data[1] /= s;
      data[2] /= s;
      return *this;
    }

    bool operator==(const vec& rhs) {
      return data[0] == rhs.data[0] && data[1] == rhs.data[1]
          && data[2] == rhs.data[2];
    }
    __host__ __device__
    void assign(fl s) {
      data[0] = data[1] = data[2] = s;
    }
    __host__   __device__ sz size() const {
      return 3;
    }

    void print(std::ostream& out) const {
      out << "<" << data[0] << "," << data[1] << "," << data[2] << ">";
    }
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      for (unsigned i = 0; i < 3; i++) //store as floats
          {
        float x = data[i];
        ar & x;
        data[i] = x;
      }
    }
};

__host__   __device__
  inline vec operator*(fl s, const vec& v) {
  return vec(s * v[0], s * v[1], s * v[2]);
}

__host__   __device__
  inline vec cross_product(const vec& a, const vec& b) {
  return vec(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

__host__   __device__
  inline vec elementwise_product(const vec& a, const vec& b) {
  return vec(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

__host__   __device__
  inline fl vec_distance_sqr(const vec& a, const vec& b) {
  return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
}

__host__   __device__
  inline fl sqr(const vec& v) {
  return sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
}

struct mat {
    fl data[9];
    __host__ __device__ mat() {
#ifndef NDEBUG
      data[0] = data[1] = data[2] = data[3] = data[4] = data[5] = data[6] =
          data[7] = data[8] = not_a_num;
#endif
    }
    // column-major
    __host__   __device__
	  const fl& operator()(sz i, sz j) const {
      assert(i < 3);
      assert(j < 3);
      return data[i + 3 * j];
    }
    __host__   __device__ fl& operator()(sz i, sz j) {
      assert(i < 3);
      assert(j < 3);
      return data[i + 3 * j];
    }

    __host__ __device__ mat(fl xx, fl xy, fl xz, fl yx, fl yy, fl yz, fl zx,
        fl zy, fl zz) {

      data[0] = xx;
      data[3] = xy;
      data[6] = xz;
      data[1] = yx;
      data[4] = yy;
      data[7] = yz;
      data[2] = zx;
      data[5] = zy;
      data[8] = zz;
    }
    __host__   __device__ vec operator*(const vec& v) const {
      return vec(data[0] * v[0] + data[3] * v[1] + data[6] * v[2],
          data[1] * v[0] + data[4] * v[1] + data[7] * v[2],
          data[2] * v[0] + data[5] * v[1] + data[8] * v[2]);
    }
    __host__   __device__
	  const mat& operator*=(fl s) {
      VINA_FOR(i, 9)
        data[i] *= s;
      return *this;
    }
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & data;
    }
};

typedef std::vector<vec> vecv;
typedef gpair<vec, vec> vecp;
typedef std::vector<fl> flv;
typedef std::vector<pr> prv;
typedef std::vector<sz> szv;
typedef boost::filesystem::path path;

//template instantiation, mostly for cdt indexer
template class std::vector<vec>;
template class gpair<vec, vec> ;
template class std::vector<fl>;
template class std::vector<pr>;
template class std::vector<sz>;

struct internal_error {
    std::string file;
    unsigned line;
    internal_error(const std::string& file_, unsigned line_)
        : file(file_), line(line_) {
    }
};

struct numerical_error : public std::runtime_error {
    numerical_error(const std::string& message)
        : std::runtime_error(message) {
    }
};

#ifdef NDEBUG
#define VINA_CHECK(P) do { if(!(P)) throw internal_error(__FILE__, __LINE__); } while(false)
#else
#define VINA_CHECK(P) assert(P)
#endif

const fl pi = fl(3.1415926535897931);

inline sz fl_to_sz(fl x, sz max_sz) { // return a value in [0, max_sz]
  if (x <= 0) return 0;
  if (x >= max_sz) return max_sz;
  sz tmp = static_cast<sz>(x);
  return (std::min)(tmp, max_sz); // sz -> fl cast loses precision. 'min' is to guard against returning values out of range
}

const fl fl_tolerance = fl(0.001);

#ifndef __CUDA_ARCH__
__host__ inline bool eq(fl a, fl b) {
  return std::abs(a - b) < fl_tolerance;
}
#else
__device__ inline bool eq(float a, float b) {
  return fabsf(a - b) < fl_tolerance;
}

__device__ inline bool eq(double a, double b) {
  return fabs(a - b) < fl_tolerance;
}
#endif

__host__ __device__ inline bool eq(const vec& a, const vec& b) {
  return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);
}

template<typename T>
bool eq(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) return false;
  VINA_FOR_IN(i, a)
    if (!eq(a[i], b[i])) return false;
  return true;
}

const fl max_fl = (std::numeric_limits<fl>::max)();
const sz max_sz = (std::numeric_limits<sz>::max)();
const unsigned max_unsigned = (std::numeric_limits<unsigned>::max)();
const fl epsilon_fl = std::numeric_limits<fl>::epsilon();

const vec zero_vec(0, 0, 0);
const vec max_vec(max_fl, max_fl, max_fl);

inline bool not_max(fl x) {
  return (x < 0.1 * max_fl);
}

template<typename T, typename A>
sz vector_append(std::vector<T, A>& x, const std::vector<T, A>& y) { // return old size
  sz old_size = x.size();
  x.insert(x.end(), y.begin(), y.end());
  return old_size;
}

template<typename T>
sz find_min(const std::vector<T>& v) { // returns v.size() i.e. 0 for empty vectors; the index of the smallest elt otherwise
  sz tmp = v.size();
  VINA_FOR_IN(i, v)
    if (i == 0 || v[i] < v[tmp]) tmp = i;
  return tmp;
}

#ifndef __CUDA_ARCH__
__host__ inline void normalize_angle(fl& x) { // subtract or add enough 2*pi's to make x be in [-pi, pi]
  if (x > 3*pi) { // very large
    fl n = ( x - pi) / (2*pi);// how many 2*pi's do you want to subtract?
    x -= 2*pi*std::ceil(n);// ceil can be very slow, but this should not be called often
    normalize_angle(x);
  }
  else if(x < -3*pi) { // very small
    fl n = (-x - pi) / (2*pi);// how many 2*pi's do you want to add?
    x += 2*pi*std::ceil(n);// ceil can be very slow, but this should not be called often
    normalize_angle(x);
  }
  else if(x > pi) { // in (   pi, 3*pi]
    x -= 2*pi;
  }
  else if(x < -pi) { // in [-3*pi,  -pi)
    x += 2*pi;
  }
  if(!(x >= -pi && x <= pi)) {
    throw numerical_error("Numerical degeneracy encountered. Check for non-physical inputs.");
  }
  // in [-pi, pi]
}
#else
__device__ inline void normalize_angle(fl& x) {
  if (x > 3 * pi) {
    fl n = (x - pi) / (2 * pi);
    x -= 2 * pi * ceilf(n);
    normalize_angle(x);
  } else
    if (x < -3 * pi) {
      fl n = (-x - pi) / (2 * pi);
      x += 2 * pi * ceilf(n);
      normalize_angle(x);
    } else
      if (x > pi) {
        x -= 2 * pi;
      } else
        if (x < -pi) {
          x += 2 * pi;
        }
  //TODO: put the assert back?
}
#endif

__host__   __device__   inline fl normalized_angle(fl x) {
  normalize_angle(x);
  return x;
}

template<typename T>
std::string to_string(const T& x, std::streamsize width = 0, char fill = ' ') { // default 0 width means no restrictions on width
  std::ostringstream out;
  out.fill(fill);
  if (width > 0) out << std::setw(width);
  out << x;
  return out.str();
}

template<typename T>
T sum(const std::vector<T>& v) {
  T acc = 0;
  VINA_FOR_IN(i, v)
    acc += v[i];
  return acc;
}

// multiply pK by this to get free energy in kcal/mol:
// K = exp(E/RT)  -- lower K and lower E == better binder
// pK = -log10(K)   => K = 10^(-pK)
// E = RT ln(K) = RT ln (10^(-pK)) = - RT * ln(10) * pK
const fl pK_to_energy_factor = -8.31 /* RT in J/K/mol */* 0.001 /* kilo */* 300 /* K */
/ 4.184 /* J/cal */* std::log(10.0); //  -0.6 kcal/mol * log(10) = -1.38

inline fl pK_to_energy(fl pK) {
  return pK_to_energy_factor * pK;
}

inline void print(fl x, std::ostream& out = std::cout) {
  out << x;
}

inline void print(sz x, std::ostream& out = std::cout) {
  out << x;
}

inline void print(const vec& v, std::ostream& out = std::cout) {
  out << "(";
  VINA_FOR_IN(i, v) {
    if (i != 0) out << ", ";
    print(v[i], out);
  }
  out << ")";
}

template<typename T>
void print(const std::vector<T>& v, std::ostream& out = std::cout) {
  out << "[";
  VINA_FOR_IN(i, v) {
    if (i != 0) out << " ";
    print(v[i], out);
  }
  out << "]";
}

template<typename T>
void printnl(const T& x, std::ostream& out = std::cout) {
  print(x, out);
  out << '\n';
}

inline bool starts_with(const std::string& str, const std::string& start) {
  return str.size() >= start.size() && str.substr(0, start.size()) == start;
}

template<typename T>
bool has(const std::vector<T>& v, const T& element) {
  return std::find(v.begin(), v.end(), element) != v.end();
}

struct usage_error : public std::runtime_error {
    usage_error(const std::string& message)
        : std::runtime_error(message) {
    }
};

#endif
