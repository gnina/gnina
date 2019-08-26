#pragma once 

#include <cmath>
#include <iostream>
#include "../lib/matrix.h"

inline void cpu_l1(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float diff = fabs(optgrid[k] - screengrid[k]);
    sum += diff;
  }
  *scoregrid = sum;
}

inline void cpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float diff = optgrid[k] - screengrid[k];
    float sqdiff = diff * diff;
    sum += sqdiff;
  }
  *scoregrid = sum;
}

inline void cpu_mult(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float weight = optgrid[k] * screengrid[k];
    sum += weight;
  }
  *scoregrid = sum;
}

inline void cpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid, 
    size_t gsize, float positive_threshold, float negative_threshold) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float threshold = optgrid[k] >=0 ? positive_threshold : negative_threshold;
    float sign = optgrid[k] >= 0 ? 1 : -1;
    float magnitude = fabs(optgrid[k]);
    float weight = ((magnitude > threshold) && screengrid[k]) ? 1 : 0;
    sum += sign * weight;
  }
  *scoregrid = sum;
}

float l1(const vec& icoords, const vec& jcoords) {
  float sum = 0.;
  for (size_t k=0; k<3; ++k) {
    sum += icoords[k] - jcoords[k];
  }
  return sum;
}

inline vec get_cube_coords(unsigned cube_idx, float dimension, unsigned blocks_per_side) {
  // wlog have cube(0,0,0) start at (0,0,0) (i.e. ignore grid center) since
  // distance is invariant to translation
  vec coords(0,0,0);
  // map flat cube index to (x,y,z) offsets
  unsigned x_offset = cube_idx / (blocks_per_side * blocks_per_side);
  unsigned y_offset = (cube_idx % (blocks_per_side * blocks_per_side)) / blocks_per_side;
  unsigned z_offset = cube_idx % blocks_per_side;
  float cube_dimension = dimension / blocks_per_side;
  coords[0] = x_offset * cube_dimension;
  coords[1] = y_offset * cube_dimension;
  coords[2] = z_offset * cube_dimension;
  return coords;
}

inline float get_cube_l1(unsigned i_cube, unsigned j_cube, float dimension, unsigned
    blocks_per_side) {
  // get (x,y,z) coords of cubes
  vec icoords = get_cube_coords(i_cube, dimension, blocks_per_side);
  vec jcoords = get_cube_coords(j_cube, dimension, blocks_per_side);
  return l1(icoords, jcoords);
}

inline void populate_cost_matrix(unsigned ncubes, unsigned subgrid_dim,
    unsigned ntypes, unsigned blocks_per_side, float dimension, flmat& cost_matrix) {
  // fill in cost matrix for comparing two spatial signatures
  unsigned n_elements = ntypes * ncubes;
  for (unsigned i=0; i<n_elements; ++i) {
    for (unsigned j=i+1; j<n_elements; ++j) {
      // within the same cube, cost is 0...assumes 0-initialized
      if (i!=j) {
        unsigned i_ch = i / ncubes;
        unsigned j_ch = j / ncubes;
        unsigned i_cube = i % ncubes;
        unsigned j_cube = j % ncubes;
        float cube_l1 = get_cube_l1(i_cube, j_cube, dimension, blocks_per_side);
        // within the same channel but different cube, cost is L1 distance between cubes
        if (i_ch == j_ch)
          cost_matrix(i,j) = cube_l1;
        // across channels, cost is constant, arbitrarily large value + L1
        else
          cost_matrix(i,j) = cube_l1 + 1e6;
      }
    }
  }
}

void do_gpu_l1(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_mult(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize, float positive_threshold, float negative_threshold);

void do_gpu_emd(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);
