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

#ifndef VINA_SZV_GRID_H
#define VINA_SZV_GRID_H

#include "model.h"
#include "grid_dim.h"
#include "array3d.h"
#include "brick.h"

#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>

namespace boost {
//overload for using array3 as hash key
inline bool operator==(const array<int, 3> &a, const array<int, 3> &b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

// hash value
inline std::size_t hash_value(const array<int, 3> &e) {
  std::size_t seed = 0;
  boost::hash_combine(seed, e[0]);
  boost::hash_combine(seed, e[1]);
  boost::hash_combine(seed, e[2]);
  return seed;
}
}

//dkoes - this is a 'global' cache of receptor atoms that are within a cutoff
//distance from global grid points; the atom lists are calculated on demand
//and stored in a hash
class szv_grid_cache {
    typedef boost::array<int, 3> ijk;
    typedef boost::unordered_map<ijk, szv*> cache_type;
    mutable cache_type cache;
    const model& m;
    fl cutoff_sqr;
    static constexpr fl granularity = 3.0; //good balance of cache locality and avoiding redundant computation
  public:
    szv_grid_cache(const model& m_, fl cut)
        : m(m_), cutoff_sqr(cut) {

    }

    ~szv_grid_cache() {
      //clear out szv vectors
      for (cache_type::iterator itr = cache.begin(), end = cache.end();
          itr != end; ++itr) {
        if (itr->second != NULL) {
          delete itr->second;
          itr->second = NULL;
        }
      }
    }

    const model& getModel() const {
      return m;
    }

    //compute all the receptor atoms that may be reachable by passed grid dims
    void compute_relevant(const grid_dims& gd, szv& relevant_indices) const {
      vec start, end;
      for (sz i = 0; i < 3; i++) {
        start[i] = gd[i].begin;
        end[i] = gd[i].end;
      }

      VINA_FOR_IN(i, m.grid_atoms) {
        const atom& a = m.grid_atoms[i];
        if (a.acceptable_type() && !a.is_hydrogen()
            && brick_distance_sqr(start, end, a.coords) < cutoff_sqr)
          relevant_indices.push_back(i);
      }
    }

    //given grid dimensions, fill an offset to adjust indicies to local
    //grid (gs) and the range of these values (range is end point, not last value)
    static void get_local_dims(const grid_dims& gd, ijk& offset, ijk& dims) {
      for (sz i = 0; i < 3; i++) {
        offset[i] = std::floor(gd[i].begin / granularity);
        dims[i] = std::ceil(gd[i].end / granularity) - offset[i] + 1;
      }
    }

    //return index of coord in local reference (offset comes from get_local_dims)
    static ijk local_index(const vec& coord, const ijk& offset) {
      ijk ret;
      for (sz i = 0; i < 3; i++) {
        ret[i] = std::floor(coord[i] / granularity) - offset[i];
      }
      return ret;
    }

    //return pointer to possibilities vector from cache
    //the value is generated on-demand looking just at the receptor
    //atoms in relvant_indices if necessary
    const szv* get(const vec& coord, const szv& relevant_indices) const {
      //get unique global index for coord
      ijk index;
      for (sz i = 0; i < 3; i++) {
        index[i] = std::floor(coord[i] / granularity);
      }

      if (cache.count(index) == 0) {
        //fill out the list of close enough receptor atoms
        szv *atoms = new szv();
        //compute lower and upper coordinates of this grid point
        vec lower, upper;
        for (sz i = 0; i < 3; i++) {
          lower[i] = std::floor(coord[i] / granularity) * granularity;
          upper[i] = std::ceil(coord[i] / granularity) * granularity;
        }
        VINA_FOR_IN(ri, relevant_indices) {
          const sz i = relevant_indices[ri];
          const atom& a = m.grid_atoms[i];
          if (!a.is_hydrogen() && a.acceptable_type()) {
            if (brick_distance_sqr(lower, upper, a.coords) < cutoff_sqr)
              atoms->push_back(i);
          }
        }
        cache[index] = atoms;
      }
      return cache[index];
    }

    //return the dimension of the grid for given dimensions
    static grid_dims szv_grid_dims(const grid_dims& gd) {
      ijk off, range;
      get_local_dims(gd, off, range);
      grid_dims tmp;
      VINA_FOR_IN(i, tmp) {
        tmp[i].begin = gd[i].begin;
        tmp[i].end = gd[i].end;
        tmp[i].n = range[i];
      }
      return tmp;
    }

};

//dkoes - this keeps track of what receptor atoms are possibly close enough
//to grid points to matter and caches their indices
struct szv_grid {
    szv_grid(szv_grid_cache& c, const grid_dims& gd)
        : cache(c) {
      cache.get_local_dims(gd, offset, range);
      m_data.resize(range[0], range[1], range[2]);

      cache.compute_relevant(gd, relevant_indexes);
      //don't precompute - this is particularly inefficient for minimization
    }

    const szv& possibilities(const vec& coords) const {
      boost::array<int, 3> index = cache.local_index(coords, offset);
      assert(index[0] < m_data.dim0());
      assert(index[1] < m_data.dim1());
      assert(index[2] < m_data.dim2());
      const szv* ret = m_data(index[0], index[1], index[2]);
      if (ret == NULL) {
        //fetch from cache
        ret = cache.get(coords, relevant_indexes);
        m_data(index[0], index[1], index[2]) = ret;
      }
      return *ret;
    }
  private:
    szv_grid_cache& cache;
    szv relevant_indexes; //rec atoms within distance of docking grid
    mutable array3d<const szv*> m_data; //this is updated as needed, does NOT own memory
    boost::array<int, 3> offset;
    boost::array<int, 3> range;

};

#endif
