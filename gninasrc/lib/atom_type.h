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

#ifndef VINA_ATOM_TYPE_H
#define VINA_ATOM_TYPE_H

#include "atom_constants.h"
#include "triangular_matrix_index.h"

struct atom_type {
    smt sm; //single smina atom type
    atom_type()
        : sm(smina_atom_type::NumTypes) {
    }
    inline smt get() const {
      return sm;
    }

    bool is_hydrogen() const {
      return ::is_hydrogen(sm);
    }
    bool is_heteroatom() const {
      return ::is_heteroatom(sm);
    }
    bool acceptable_type() const {
      return sm < smina_atom_type::NumTypes;
    }

    bool same_element(const atom_type& a) const { // does not distinguish metals or unassigned types
      return smina_atom_type::data[sm].anum == smina_atom_type::data[a.sm].anum;
    }
    fl covalent_radius() const {
      return ::covalent_radius(sm);
    }
    fl optimal_covalent_bond_length(const atom_type& x) const {
      return covalent_radius() + x.covalent_radius();
    }
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      unsigned char c = (int) sm;
      assert((int)sm < UCHAR_MAX);
      ar & c;
      sm = (smt) c;
    }
};

inline sz num_atom_types() {
  return smina_atom_type::NumTypes;
}

inline sz get_type_pair_index(const atom_type& a, const atom_type& b) { // throws error if any arg is unassigned in the given typing scheme
  sz n = num_atom_types();

  smt i = a.get();
  VINA_CHECK(i < n);
  smt j = b.get();
  VINA_CHECK(j < n);

  if (i <= j)
    return triangular_matrix_index(n, i, j);
  else
    return triangular_matrix_index(n, j, i);
}

#endif
