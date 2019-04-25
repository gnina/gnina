#pragma once 

#include <cmath>
#include <iostream>

inline void cpu_l2(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float diff = optgrid[k] - screengrid[k];
    float sqdiff = diff * diff;
    sum += sqdiff;
  }
  *scoregrid = std::sqrt(sum);
}

