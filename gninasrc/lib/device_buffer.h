/*
 * device_buffer.h
 *
 * Arena allocator for GPU memory. Avoids repeated allocation calls by
 * managing one large slab per CPU thread.
 *
 * Backends:
 *   USE_METAL=ON  → Metal MTLBuffer (MTLStorageModeShared, unified memory).
 *                   Contents pointer is valid for both CPU and GPU on Apple
 *                   Silicon — no explicit copies needed.
 *   USE_METAL=OFF → CUDA cudaMalloc (original behaviour, unchanged).
 */

#ifndef __DEVICE_BUFFER_H
#define __DEVICE_BUFFER_H

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstddef>

#ifdef USE_METAL
  #include "cuda_metal_compat.h"   // provides cudaMemcpyKind, cudaError_t, cudaSuccess
#else
  #include <cuda_runtime.h>
#endif

// Returns an estimate of safe allocatable GPU memory per CPU thread.
size_t available_mem(size_t num_cpu_threads);

// ─────────────────────────────────────────────────────────────────────────────
// device_buffer — bump-pointer arena for GPU memory
// ─────────────────────────────────────────────────────────────────────────────
class device_buffer {
    char*   begin;       // start of the slab (host-visible on Metal/AppleSilicon)
    size_t  capacity;    // total slab size in bytes
    char*   next_alloc;  // bump pointer (next free byte)

#ifdef USE_METAL
    void*   _metal_buffer_handle;  // CFBridgingRetained id<MTLBuffer>
#endif

    void*        copy_bytes(void* from, size_t n_bytes, cudaMemcpyKind kind);
    cudaError_t  alloc_bytes(void** alloc, size_t n_bytes);
    bool         has_space(size_t n_bytes);

public:
#ifdef USE_METAL
    // Returns the underlying CFBridgingRetained id<MTLBuffer> as a void*.
    // Cast to MTLBufferHandle (from metal_context.h) before passing to dispatch1D.
    void* metalBufferHandle() const { return _metal_buffer_handle; }

    // Byte offset of a sub-allocation from the start of the Metal buffer's slab.
    // Use to compute the 'offset' parameter for MetalArg::fromBuffer().
    size_t offsetOf(const void* ptr) const {
        return static_cast<size_t>(static_cast<const char*>(ptr) - begin);
    }
#endif

public:
    device_buffer();
    ~device_buffer();

    // Allocate the slab (call once after constructing, before any alloc/copy).
    void init(size_t capacity);

    // Grow the slab. Only valid when the buffer is empty (next_alloc == begin).
    void resize(size_t n_bytes);

    // Reset bump pointer — reuse the slab without freeing it.
    void reinitialize() { next_alloc = begin; }

    // Allocate n_bytes * sizeof(T) from the slab. Returns device pointer.
    template<typename T>
    cudaError_t alloc(T** alloc, size_t n_bytes) {
        return alloc_bytes((void**)alloc, n_bytes);
    }

    // Allocate and copy n_requested elements of type T.
    template<typename T>
    T* copy(T* cpu_object, size_t n_requested, cudaMemcpyKind kind) {
        return (T*)copy_bytes(cpu_object, n_requested * sizeof(T), kind);
    }

    // No-op "free" — deallocation happens en-masse via reinitialize().
    template<typename T>
    cudaError_t dealloc(T* /*alloc*/) {
        return cudaSuccess;
    }
};

// One arena per CPU thread.
extern thread_local device_buffer thread_buffer;

// ─────────────────────────────────────────────────────────────────────────────
// Global helpers — thin wrappers that forward to thread_buffer.
// ─────────────────────────────────────────────────────────────────────────────
cudaError_t device_alloc_bytes(void** alloc, size_t n_bytes);
cudaError_t device_free(void* buf);

template<typename T>
cudaError_t device_malloc(T** alloc, size_t n_bytes) {
    return device_alloc_bytes((void**)alloc, n_bytes);
}

#endif // __DEVICE_BUFFER_H
