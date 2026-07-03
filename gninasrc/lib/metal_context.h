/*
 * metal_context.h
 *
 * Singleton that owns the MTLDevice and provides per-thread MTLCommandQueues.
 *
 * Design principles for Apple Silicon (M1/M2/M3):
 *  - One MTLDevice per process (created once via MTLCreateSystemDefaultDevice).
 *  - One MTLCommandQueue per CPU thread (thread_local, lazy init).
 *  - MTLStorageModeShared buffers: CPU and GPU share the same physical memory,
 *    so there are NO explicit host↔device copies — a plain memcpy suffices.
 *
 * Usage from C++ code:
 *   #include "metal_context.h"
 *   MetalContext& ctx = MetalContext::instance();
 *   MTL::CommandQueue* q = ctx.commandQueue();  // thread-local
 */

#pragma once

#ifdef USE_METAL
#ifdef __APPLE__

// Use the Metal C++ headers (metal-cpp) when available, otherwise fall back
// to Objective-C forward declarations wrapped in an opaque handle.
// We use the opaque-handle approach here to keep this header pure C++
// (implementation is in the .mm file).

#include <cstddef>

// Opaque handle types — the .mm file casts these to the real ObjC types.
struct _MTLDevice_opaque;
struct _MTLCommandQueue_opaque;
struct _MTLBuffer_opaque;
struct _MTLCommandBuffer_opaque;
struct _MTLComputePipelineState_opaque;

using MTLDeviceHandle                = _MTLDevice_opaque*;
using MTLCommandQueueHandle          = _MTLCommandQueue_opaque*;
using MTLBufferHandle                = _MTLBuffer_opaque*;
using MTLComputePipelineStateHandle  = _MTLComputePipelineState_opaque*;

// ─────────────────────────────────────────────────────────────────────────────
// MetalArg — one argument to a compute kernel dispatch.
//
// Either a slice of an MTLBuffer (large arrays) or inline bytes (small
// constants such as a single float or uint).  Corresponds to buffer(N) in MSL.
// ─────────────────────────────────────────────────────────────────────────────
struct MetalArg {
    bool             is_buffer;
    MTLBufferHandle  buf;        // valid when is_buffer == true
    size_t           offset;     // byte offset into buf
    const void*      bytes;      // valid when is_buffer == false
    size_t           length;     // byte length of inline data

    // Construct from a buffer slice.
    static MetalArg fromBuffer(MTLBufferHandle b, size_t off = 0) {
        MetalArg a;
        a.is_buffer = true;
        a.buf       = b;
        a.offset    = off;
        a.bytes     = nullptr;
        a.length    = 0;
        return a;
    }

    // Construct from a small value held by the caller (copied inline).
    static MetalArg fromBytes(const void* ptr, size_t len) {
        MetalArg a;
        a.is_buffer = false;
        a.buf       = nullptr;
        a.offset    = 0;
        a.bytes     = ptr;
        a.length    = len;
        return a;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// MetalContext
// ─────────────────────────────────────────────────────────────────────────────
class MetalContext {
public:
    // Singleton accessor — thread-safe (C++11 static local).
    static MetalContext& instance();

    // The system default MTLDevice (the M1 GPU).
    MTLDeviceHandle device() const { return _device; }

    // Per-thread command queue (lazy init on first call from each thread).
    MTLCommandQueueHandle commandQueue();

    // Allocate a shared-mode Metal buffer of `bytes` bytes.
    // On Apple Silicon, contents() returns a CPU pointer that is also valid
    // on the GPU — no explicit copy is needed.
    MTLBufferHandle newBuffer(size_t bytes);

    // Release a Metal buffer previously obtained from newBuffer().
    static void releaseBuffer(MTLBufferHandle buf);

    // CPU-accessible pointer into a Metal buffer (just buf->contents()).
    static void* bufferContents(MTLBufferHandle buf);

    // Commit a command buffer and wait for GPU completion (used for
    // synchronisation points between host-dispatch sub-passes).
    void synchronize();

    // True if Metal is available on this machine.
    static bool available();

    // ── Compute dispatch ─────────────────────────────────────────────────────

    // Create a compute pipeline for a named MSL function from the embedded
    // gnina_kernels.metallib.  The returned handle is owned by the caller;
    // release with releasePipeline() when done (or hold for the process life).
    MTLComputePipelineStateHandle makePipeline(const char* functionName);

    // Release a pipeline state returned by makePipeline().
    static void releasePipeline(MTLComputePipelineStateHandle pso);

    // Dispatch a 1-D compute kernel and wait for completion.
    //
    //   pso         — pipeline state (from makePipeline)
    //   args        — array of nArgs MetalArg descriptors (one per buffer(N))
    //   nArgs       — length of args[]
    //   threadCount — total number of GPU threads to launch
    //
    // The function commits and waits synchronously so the caller can read
    // results immediately after returning.
    void dispatch1D(MTLComputePipelineStateHandle pso,
                    const MetalArg* args,
                    int nArgs,
                    size_t threadCount);

    // Deleted copy/move — singleton only.
    MetalContext(const MetalContext&)            = delete;
    MetalContext& operator=(const MetalContext&) = delete;

private:
    MetalContext();
    ~MetalContext();

    // Load the embedded gnina_kernels.metallib (called once from constructor).
    void _loadLibrary();

    MTLDeviceHandle  _device;
    void*            _library;   // CFBridgingRetained id<MTLLibrary> or nullptr
};

// ─────────────────────────────────────────────────────────────────────────────
// Free helper: returns the thread-local command queue (shorthand).
// ─────────────────────────────────────────────────────────────────────────────
MTLCommandQueueHandle metalThreadCommandQueue();

#endif // __APPLE__
#endif // USE_METAL
