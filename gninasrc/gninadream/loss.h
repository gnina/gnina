#pragma once 

#include <cmath>
#include <iostream>

inline void cpu_l1(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize) {
  float sum = 0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t k=0; k<gsize; ++k) {
    float diff = optgrid[k] - screengrid[k];
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

void do_gpu_l1(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_l2sq(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_mult(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize);

void do_gpu_thresh(const float* optgrid, const float* screengrid, float* scoregrid, size_t gsize, float positive_threshold, float negative_threshold);
