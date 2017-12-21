/*
 * Interface for reusing a single memory buffer 
 * to avoid repeated calls to cudaMalloc
 */

#ifndef __DEVICE_BUFFER_H
#define __DEVICE_BUFFER_H
#include <cuda_runtime.h>

template<typename T>
struct device_buffer
{
    T* begin; //pointer to the beginning of the buffer
    size_t capacity; //buffer size in bytes is capacity*sizeof(type),
    T* current;  //pointer to the beginning of the unused region

    device_buffer(size_t capacity=10000);

    bool has_space(size_t n_requested);
    void resize(size_t new_requested);
    T* copy(T* cpu_object, size_t n_requested, cudaMemcpyKind kind); //returns ptr to segment
    void reinitialize() {current = begin;}

    ~device_buffer();
};

typedef device_buffer<float> float_buffer;
#endif
