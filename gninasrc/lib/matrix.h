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

#ifndef VINA_MATRIX_H
#define VINA_MATRIX_H

#include <vector>
#include "triangular_matrix_index.h"
#include <cuda_runtime.h>
#include "gpu_util.h"
#include "device_buffer.h"

// these 4 lines are used 3 times verbatim - defining a temp macro to ease the pain
#define VINA_MATRIX_DEFINE_OPERATORS \
	const T& operator()(sz i) const { return m_data[i]; } \
	      T& operator()(sz i)       { return m_data[i]; } \
	const T& operator()(sz i, sz j) const { return m_data[index(i, j)]; } \
	      T& operator()(sz i, sz j)       { return m_data[index(i, j)]; } 

template<typename T>
class matrix {
    std::vector<T> m_data;
    sz m_i, m_j;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & m_data;
      ar & m_i;
      ar & m_j;
    }
  public:
    sz index(sz i, sz j) const {
      assert(j < m_j);
      assert(i < m_i);
      return i + m_i * j; // column-major
    }
    matrix()
        : m_i(0), m_j(0) {
    }
    matrix(sz i, sz j, const T& filler_val)
        : m_data(i * j, filler_val), m_i(i), m_j(j) {
    }
    void resize(sz m, sz n, const T& filler_val) { // new sizes should be the same or greater than the old
      if (m == dim_1() && n == dim_2()) return; // no-op
      VINA_CHECK(m >= dim_1());
      VINA_CHECK(n >= dim_2());
      std::vector<T> tmp(m * n, filler_val);
      VINA_FOR(i, m_i)
        VINA_FOR(j, m_j)
          tmp[i + m * j] = (*this)(i, j);
      m_data = tmp;
      m_i = m;
      m_j = n;
    }
    void append(const matrix<T>& x, const T& filler_val) {
      sz m = dim_1();
      sz n = dim_2();
      resize(m + x.dim_1(), n + x.dim_2(), filler_val);
      VINA_FOR(i, x.dim_1())
        VINA_FOR(j, x.dim_2())
          (*this)(i + m, j + n) = x(i, j);
    }
    VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above
    sz dim_1() const {
      return m_i;
    }
    sz dim_2() const {
      return m_j;
    }
};

template<typename T>
class triangular_matrix_gpu;

template<typename T>
class triangular_matrix {
    std::vector<T> m_data;
    sz m_dim;
  public:
    sz index(sz i, sz j) const {
      return triangular_matrix_index(m_dim, i, j);
    }
    sz index_permissive(sz i, sz j) const {
      return (i < j) ? index(i, j) : index(j, i);
    }
    triangular_matrix()
        : m_dim(0) {
    }
    triangular_matrix(sz n, const T& filler_val)
        : m_data(n * (n + 1) / 2, filler_val), m_dim(n) {
    }
    VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above
    sz dim() const {
      return m_dim;
    }
    friend triangular_matrix_gpu<T> ;
};

typedef triangular_matrix<fl> flmat;

template<typename T>
struct triangular_matrix_gpu {
    T* m_data;
    sz m_dim;

    __device__ sz index(sz i, sz j) const {
      return triangular_matrix_index(m_dim, i, j);
    }
    __device__ sz index_permissive(sz i, sz j) const {
      return (i < j) ? index(i, j) : index(j, i);
    }
    triangular_matrix_gpu(sz n)
        : m_data(NULL), m_dim(n) {
      triangular_matrix<T> cpu_matrix(n, 0);
      set_diagonal(cpu_matrix, 1);
      CUDA_CHECK_GNINA(
          device_malloc(&m_data, sizeof(T) * cpu_matrix.m_data.size()));
      CUDA_CHECK_GNINA(
          definitelyPinnedMemcpy(m_data, &cpu_matrix.m_data[0],
              sizeof(T) * cpu_matrix.m_data.size(), cudaMemcpyHostToDevice));
    }

    ~triangular_matrix_gpu() {
      device_free(m_data);
    }

    __device__  const T& operator()(sz i) const {
      return m_data[i];
    }
    __device__ T& operator()(sz i) {
      return m_data[i];
    }
    __device__  const T& operator()(sz i, sz j) const {
      return m_data[index(i, j)];
    }
    __device__ T& operator()(sz i, sz j) {
      return m_data[index(i, j)];
    }

    __device__ sz dim() const {
      return m_dim;
    }
};

typedef triangular_matrix_gpu<fl> flmat_gpu;

template<typename T>
class strictly_triangular_matrix {
    std::vector<T> m_data;
    sz m_dim;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & m_data;
      ar & m_dim;
    }
  public:
    sz index(sz i, sz j) const {
      assert(j < m_dim);
      assert(i < j);
      assert(j >= 1); // by implication, really
      return i + j * (j - 1) / 2;
    }
    sz index_permissive(sz i, sz j) const {
      return (i < j) ? index(i, j) : index(j, i);
    }
    strictly_triangular_matrix()
        : m_dim(0) {
    }
    strictly_triangular_matrix(sz n, const T& filler_val)
        : m_data(n * (n - 1) / 2, filler_val), m_dim(n) {
    }
    void resize(sz n, const T& filler_val) {
      if (n == m_dim) return; // no-op
      VINA_CHECK(n > m_dim);
      m_dim = n;
      m_data.resize(n * (n - 1) / 2, filler_val); // preserves original data
    }
    void append(const strictly_triangular_matrix<T>& m, const T& filler_val) {
      sz n = dim();
      resize(n + m.dim(), filler_val);
      VINA_FOR(i, m.dim())
        VINA_RANGE(j, i+1, m.dim())
          (*this)(i + n, j + n) = m(i, j);
    }
    void append(const matrix<T>& rectangular,
        const strictly_triangular_matrix<T>& triangular) {
      VINA_CHECK(dim() == rectangular.dim_1());
      VINA_CHECK(rectangular.dim_2() == triangular.dim());

      // a filler value is needed by append or resize
      // we will use a value from rectangular as the filler value
      // but it can not be obtained if dim_1 or dim_2 is 0
      // these cases have to be considered separately
      if (rectangular.dim_2() == 0) return;
      if (rectangular.dim_1() == 0) {
        (*this) = triangular;
        return;
      }
      const T& filler_val = rectangular(0, 0); // needed by 'append below'

      sz n = dim();
      append(triangular, filler_val);
      VINA_FOR(i, rectangular.dim_1())
        VINA_FOR(j, rectangular.dim_2())
          (*this)(i, n + j) = rectangular(i, j);
    }
    VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above
    sz dim() const {
      return m_dim;
    }
};

#undef VINA_MATRIX_DEFINE_OPERATORS

#endif
