/*
 * Interface for reusing a single memory buffer 
 * to avoid repeated calls to cudaMalloc
 */
#include "device_buffer.h"
#include <cmath>
#include "gpu_util.h"
#include <cassert>

template<typename T>
device_buffer<T>::device_buffer(size_t capacity) : capacity(capacity){
    CUDA_CHECK_GNINA(cudaMalloc(&begin, capacity*sizeof(T)));
    current = begin;
}

template<typename T>
bool device_buffer<T>::has_space(size_t n_requested) {
    size_t space_used = current - begin; //in elements
    return capacity - space_used >= n_requested;
}

template<typename T>
void device_buffer<T>::resize(size_t n_requested) {
    if (begin == current && capacity >= n_requested)
        return;
    //This doesn't memcpy the old contents because of the baked-in assumption
    //that you are resizing before you start copying. 
    else if (!has_space(n_requested)) {
        int factor = n_requested % capacity ? ((n_requested/capacity)+1)*capacity: n_requested;
        size_t new_capacity = capacity * factor;
        CUDA_CHECK_GNINA(cudaFree(begin));
        CUDA_CHECK_GNINA(cudaMalloc(&begin, sizeof(T)*new_capacity));
        capacity = new_capacity;
    }
}

template<typename T>
T* device_buffer<T>::copy(T* cpu_object, size_t n_requested, cudaMemcpyKind kind) {
    //N.B. you need to resize appropriately before starting to copy into the
    //buffer. This function avoids resizing so data structures in the buffer can use
    //pointers, and the only internal protection is the following assert.
    assert(has_space(n_requested));
    CUDA_CHECK_GNINA(cudaMemcpy(current, cpu_object, n_requested*sizeof(T), kind));
    T* old_current = current;
    current += n_requested;
    return old_current;
}

template<typename T>
device_buffer<T>::~device_buffer() {
    CUDA_CHECK_GNINA(cudaFree(begin));
}

template struct device_buffer<float>;
