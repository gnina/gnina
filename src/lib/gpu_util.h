#ifndef GPU_UTIL_H
#define GPU_UTIL_H

#define shared __shared__
#define global __global__
#define device __device__
#define host __host__

#include "gpu_math.h"
#include <stdio.h>

static inline void abort_on_gpu_err(void){
    cudaError err = cudaGetLastError();
    if (cudaSuccess != err)
    {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, cudaGetErrorString(err));
        exit(-1);
    }
}

#endif
