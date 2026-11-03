/*
 * device_buffer.cpp
 *
 * Arena allocator backed by:
 *   USE_METAL=ON  → Metal MTLBuffer (unified memory, Apple Silicon).
 *   USE_METAL=OFF → CUDA cudaMalloc (original implementation).
 *
 * On Apple Silicon, MTLStorageModeShared buffers occupy a single physical
 * memory region accessible from both the CPU and the GPU.  There is no DMA
 * transfer cost — memcpy() between a host pointer and the buffer's contents()
 * pointer is the fastest (and only) copy needed.
 */

#include "device_buffer.h"
#include "gpu_util.h"
#include <cassert>
#include <cmath>
#include <iostream>

// ─────────────────────────────────────────────────────────────────────────────
// Common alignment helpers (shared by both backends)
// ─────────────────────────────────────────────────────────────────────────────
#define _align_down_pow2(n, size) \
    ((decltype(n))((uintptr_t)(n) & ~((size) - 1)))
#define _align_up_pow2(n, size) \
    ((decltype(n))_align_down_pow2((uintptr_t)(n) + (size) - 1, size))

// ─────────────────────────────────────────────────────────────────────────────
// Common: has_space, alloc_bytes, copy_bytes (independent of backend)
// ─────────────────────────────────────────────────────────────────────────────
bool device_buffer::has_space(size_t n_bytes) {
    return (size_t)(next_alloc - begin) + n_bytes <= capacity;
}

cudaError_t device_buffer::alloc_bytes(void** alloc, size_t n_bytes) {
    if (!has_space(n_bytes)) {
        std::cerr << "[gnina] device_buffer: alloc of " << n_bytes
                  << " bytes failed — "
                  << (capacity - (size_t)(next_alloc - begin))
                  << "/" << capacity << " bytes free\n";
        std::abort();
    }
    *alloc     = next_alloc;
    next_alloc = _align_up_pow2(next_alloc + n_bytes, 128);
    return cudaSuccess;
}

void* device_buffer::copy_bytes(void* cpu_object, size_t n_bytes,
                                 cudaMemcpyKind kind) {
    assert(has_space(n_bytes));
    void* r;
    alloc_bytes(&r, n_bytes);

#ifdef USE_METAL
    // Unified memory: the slab is already CPU-accessible, so memcpy suffices
    // for all transfer directions (H2D, D2H, D2D are all the same physical RAM).
    (void)kind;
    memcpy(r, cpu_object, n_bytes);
#else
    CUDA_CHECK_GNINA(definitelyPinnedMemcpy(r, cpu_object, n_bytes, kind));
#endif
    return r;
}

// ─────────────────────────────────────────────────────────────────────────────
#ifdef USE_METAL
// ─────────────────────────────────────────────────────────────────────────────
// Metal backend
// ─────────────────────────────────────────────────────────────────────────────
#import <Metal/Metal.h>
#include "metal_context.h"

size_t available_mem(size_t num_cpu_threads) {
    // MTLDevice.recommendedMaxWorkingSetSize is Apple's guidance for safe GPU
    // memory use (typically ~65% of total GPU memory).
    MTLDeviceHandle devH = MetalContext::instance().device();
    id<MTLDevice> dev = (__bridge id<MTLDevice>)devH;
    size_t budget = (size_t)dev.recommendedMaxWorkingSetSize;
    return (size_t)((double)budget / (double)num_cpu_threads * 0.8);
}

thread_local device_buffer thread_buffer;

device_buffer::device_buffer()
    : begin(nullptr), capacity(0), next_alloc(nullptr),
      _metal_buffer_handle(nullptr) {}

void device_buffer::init(size_t cap) {
    capacity = cap;

    MTLDeviceHandle devH = MetalContext::instance().device();
    id<MTLDevice> dev = (__bridge id<MTLDevice>)devH;

    id<MTLBuffer> buf = [dev newBufferWithLength:cap
                                        options:MTLResourceStorageModeShared];
    if (!buf) {
        std::cerr << "[gnina] Metal: failed to allocate buffer of "
                  << cap << " bytes\n";
        std::abort();
    }
    // CFBridgingRetain transfers ownership to our void* handle.
    _metal_buffer_handle = (void*)CFBridgingRetain(buf);
    begin      = (char*)[buf contents];
    next_alloc = begin;
}

void device_buffer::resize(size_t n_bytes) {
    assert(begin == next_alloc &&
           "device_buffer: resize only allowed when buffer is empty");
    if (n_bytes > capacity) {
        // Release old buffer.
        if (_metal_buffer_handle) {
            CFBridgingRelease(_metal_buffer_handle);
            _metal_buffer_handle = nullptr;
        }
        capacity   = 0;
        begin      = nullptr;
        next_alloc = nullptr;
        init(n_bytes);
    }
}

device_buffer::~device_buffer() {
    if (_metal_buffer_handle) {
        CFBridgingRelease(_metal_buffer_handle);
        _metal_buffer_handle = nullptr;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
#else
// ─────────────────────────────────────────────────────────────────────────────
// CUDA backend (original implementation, unchanged)
// ─────────────────────────────────────────────────────────────────────────────
#include <cuda.h>
#include <boost/thread/thread.hpp>

size_t available_mem(size_t num_cpu_threads) {
    size_t free_mem, total;
    cudaError_t res = cudaMemGetInfo(&free_mem, &total);
    if (res != cudaSuccess) {
        std::cerr << "cudaMemGetInfo returned status " << res << "\n";
        return 1;
    }
    return (size_t)((double)free_mem / (double)num_cpu_threads * 0.8);
}

thread_local device_buffer thread_buffer;

device_buffer::device_buffer()
    : begin(nullptr), capacity(0), next_alloc(nullptr) {}

void device_buffer::init(size_t cap) {
    capacity = cap;
    CUDA_CHECK_GNINA(cudaMalloc(&begin, cap));
    next_alloc = begin;
}

void device_buffer::resize(size_t n_bytes) {
    assert(begin == next_alloc &&
           "device_buffer: resize only allowed when buffer is empty");
    if (n_bytes > capacity) {
        CUDA_CHECK_GNINA(cudaFree(begin));
        CUDA_CHECK_GNINA(cudaMalloc(&begin, n_bytes));
        capacity   = n_bytes;
        next_alloc = begin;
    }
}

device_buffer::~device_buffer() {
    CUDA_CHECK_GNINA(cudaFree(begin));
}

#endif // USE_METAL

// ─────────────────────────────────────────────────────────────────────────────
// Global helpers (both backends)
// ─────────────────────────────────────────────────────────────────────────────

cudaError_t device_alloc_bytes(void** alloc, size_t n_bytes) {
    return thread_buffer.alloc((char**)alloc, n_bytes);
}

cudaError_t device_free(void* /*buf*/) {
    // Deallocation is deferred: the arena is reset en-masse by reinitialize().
    return cudaSuccess;
}
