/*
 * gnina_kernels.metal
 *
 * Metal compute kernels for gnina (Apple Silicon port).
 *
 * Phase 2 kernels:
 *   scalar_mult             — element-wise array * scalar
 *
 * Phase 4 kernels (non-cached pairwise interaction energy):
 *   noncache_pair_energy    — accumulate (lig,rec) pair forces/energies
 *   noncache_postprocess    — curl + OOB penalty + total energy via atomic add
 *
 * Phase 5 kernels (tree derivative):
 *   tree_accum_forces       — per-atom force/torque accumulation onto owning node
 *   tree_propagate_layer    — bottom-up propagation for one tree layer
 *   tree_write_change       — project force_torques to change_gpu DOFs
 */

#include <metal_stdlib>
using namespace metal;

// ── Phase 2: scalar multiply ──────────────────────────────────────────────────
//
// Buffer layout:
//   0  device float*    vals — array to scale in-place
//   1  constant float&  mult — scale factor
//   2  constant uint&   n    — element count
kernel void scalar_mult(
    device   float*  vals [[ buffer(0) ]],
    constant float&  mult [[ buffer(1) ]],
    constant uint&   n    [[ buffer(2) ]],
    uint tid              [[ thread_position_in_grid ]])
{
    if (tid < n) vals[tid] *= mult;
}

// ── Phase 4: non-cached pairwise interaction energy ───────────────────────────

// Matches GPUSplineInfo_MSL in gpucode.h (16 bytes, no pointers).
struct SplineInfoMSL {
    uint  n;         // number of active components (1..4)
    float fraction;  // knot bin width
    float cutoff;    // cutoff distance
    uint  _pad;
};

// Scalar dispatch arguments (matches NoncacheArgs_MSL in gpucode.cu).
struct NoncacheArgs {
    uint  num_movable_atoms;
    uint  nrec_atoms;
    float cutoff_sq;
    float slope;
    float gbx, gby, gbz;  // gridbegins
    uint  spline_numc;     // components per type-pair
    float gex, gey, gez;  // gridends
    uint  _pad;
};

// Evaluate one cubic spline component.
// flat_knots: packed floats, 5 per knot (x, a, b, c, d).
// knot_start: index of the first knot (units of "1 knot = 5 floats").
static inline float eval_spline_msl(
    device const float* flat_knots,
    uint                knot_start,
    float               r,
    float               fraction,
    float               cutoff,
    thread float&       deriv)
{
    deriv = 0.0f;
    if (r >= cutoff || r < 0.0f) return 0.0f;
    uint  idx  = uint(r / fraction);
    uint  base = (knot_start + idx) * 5u;
    float x    = flat_knots[base];
    float a    = flat_knots[base + 1u];
    float b    = flat_knots[base + 2u];
    float c    = flat_knots[base + 3u];
    float d    = flat_knots[base + 4u];
    float lx   = r - x;
    deriv = (3.0f * a * lx + 2.0f * b) * lx + c;
    return ((a * lx + b) * lx + c) * lx + d;
}

// Compute energy and dor (derivative / r) for one atom pair.
// spline_offsets[tindex * numc + component] = knot_start for that spline.
static float eval_deriv_msl(
    device const SplineInfoMSL* splineInfo,
    device const float*         flat_knots,
    device const uint*          spline_offsets,
    uint                        numc,
    uint t,  float charge,
    uint rt, float rcharge,
    float r2, thread float& dor)
{
    float r = sqrt(r2);
    uint  t1, t2;
    float charge1, charge2;
    if (t < rt) {
        t1 = t;  t2 = rt; charge1 = fabs(charge);  charge2 = fabs(rcharge);
    } else {
        t1 = rt; t2 = t;  charge1 = fabs(rcharge); charge2 = fabs(charge);
    }
    uint tindex = t1 + t2 * (t2 + 1u) / 2u;
    SplineInfoMSL info = splineInfo[tindex];
    uint  n        = info.n;
    float fraction = info.fraction;
    float cutoff   = info.cutoff;
    uint  base_off = tindex * numc;

    float ret = 0.0f, d = 0.0f, val, deriv;
    if (n > 0u) {
        val = eval_spline_msl(flat_knots, spline_offsets[base_off],        r, fraction, cutoff, deriv);
        ret += val;           d += deriv;
        if (n > 1u) {
            val = eval_spline_msl(flat_knots, spline_offsets[base_off+1u], r, fraction, cutoff, deriv);
            ret += val * charge1; d += deriv * charge1;
            if (n > 2u) {
                val = eval_spline_msl(flat_knots, spline_offsets[base_off+2u], r, fraction, cutoff, deriv);
                ret += val * charge2; d += deriv * charge2;
                if (n > 3u) {
                    val = eval_spline_msl(flat_knots, spline_offsets[base_off+3u], r, fraction, cutoff, deriv);
                    ret += val * charge1 * charge2; d += deriv * charge1 * charge2;
                }
            }
        }
    }
    dor = d / r;
    return ret;
}

// Curl correction — matches the GPU curl() in curl.h.
static inline void curl_msl(thread float& e, thread float3& f, float v) {
    if (e > 0.0f) {
        float tmp = v / (v + e);
        e *= tmp;
        tmp *= tmp;
        f *= tmp;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// noncache_pair_energy  (Phase 4, pass 1)
//
// One threadgroup per lig atom (NC_TGROUP threads per group).
// Each group iterates over ALL receptor atoms with stride NC_TGROUP, then
// reduces locally and writes force+energy into out[l] via atomic float add.
// out[] must be pre-zeroed by the caller.
//
// Buffer layout:
//   0  constant NoncacheArgs     — scalar args
//   1  device  float4[]          — lig atom_params (coords.xyz + charge)
//   2  device  float4[]          — out: force_energy_tup (force.xyz + energy)
//   3  device  uint[]            — lig atom types
//   4  device  float4[]          — rec atom_params (coords.xyz + charge)
//   5  device  uint[]            — rec atom types
//   6  device  SplineInfoMSL[]   — one per type-pair (tindex order)
//   7  device  float[]           — flat knot data (5 floats per knot: x,a,b,c,d)
//   8  device  uint[]            — offset table [tindex * numc + component]
// ─────────────────────────────────────────────────────────────────────────────
#define NC_TGROUP 256u

kernel void noncache_pair_energy(
    constant NoncacheArgs&      args           [[ buffer(0) ]],
    device const float4*        lig_atoms      [[ buffer(1) ]],
    device float4*              out            [[ buffer(2) ]],
    device const uint*          lig_types      [[ buffer(3) ]],
    device const float4*        rec_atoms      [[ buffer(4) ]],
    device const uint*          rec_types      [[ buffer(5) ]],
    device const SplineInfoMSL* splineInfo     [[ buffer(6) ]],
    device const float*         flat_knots     [[ buffer(7) ]],
    device const uint*          spline_offsets [[ buffer(8) ]],
    uint lid [[ thread_index_in_threadgroup ]],   // 0 .. NC_TGROUP-1
    uint gid [[ threadgroup_position_in_grid ]])  // = lig atom index
{
    uint l = gid;
    if (l >= args.num_movable_atoms) return;
    uint t = lig_types[l];
    if (t <= 1u) return;  // skip hydrogens

    float4 la         = lig_atoms[l];
    float3 lig_xyz    = la.xyz;
    float  lig_charge = la.w;

    float3 local_force  = float3(0.0f);
    float  local_energy = 0.0f;

    // Each thread strides over receptor atoms.
    for (uint r = lid; r < args.nrec_atoms; r += NC_TGROUP) {
        float4 ra   = rec_atoms[r];
        float3 diff = lig_xyz - ra.xyz;
        float  rSq  = dot(diff, diff);
        if (rSq < args.cutoff_sq) {
            float dor;
            float e = eval_deriv_msl(splineInfo, flat_knots, spline_offsets,
                                     args.spline_numc,
                                     t, lig_charge, rec_types[r], ra.w,
                                     rSq, dor);
            local_force  += diff * dor;
            local_energy += e;
        }
    }

    // Two-level threadgroup reduction: simd_sum → shared scratch → simd_sum.
    threadgroup float3 tg_force [NC_TGROUP / 32u];
    threadgroup float  tg_energy[NC_TGROUP / 32u];

    local_force  = simd_sum(local_force);
    local_energy = simd_sum(local_energy);

    uint warp = lid >> 5u;
    uint lane = lid & 31u;
    if (lane == 0u) {
        tg_force [warp] = local_force;
        tg_energy[warp] = local_energy;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (warp == 0u) {
        local_force  = (lane < (NC_TGROUP / 32u)) ? tg_force [lane] : float3(0.0f);
        local_energy = (lane < (NC_TGROUP / 32u)) ? tg_energy[lane] : 0.0f;
        local_force  = simd_sum(local_force);
        local_energy = simd_sum(local_energy);

        if (lane == 0u) {
            // Atomic add: force_energy_tup layout = float4 (force.xyz + energy).
            device atomic_float* af = (device atomic_float*)(&out[l]);
            atomic_fetch_add_explicit(af,     local_force.x,  memory_order_relaxed);
            atomic_fetch_add_explicit(af + 1, local_force.y,  memory_order_relaxed);
            atomic_fetch_add_explicit(af + 2, local_force.z,  memory_order_relaxed);
            atomic_fetch_add_explicit(af + 3, local_energy,   memory_order_relaxed);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// noncache_postprocess  (Phase 4, pass 2)
//
// One thread per lig atom. Applies curl + out-of-bounds penalty, then
// atomically accumulates each atom's energy into total_energy[0].
// The C++ caller writes total_energy[0] into out[0].energy after the dispatch.
//
// Buffer layout:
//   0  constant NoncacheArgs  — scalar args
//   1  device  float4[]       — lig atom_params (for OOB coord check)
//   2  device  float4[]       — out: force_energy_tup (in/out)
//   3  constant float&        — v (curl strength)
//   4  device  atomic_float*  — total_energy accumulator (pre-zeroed, 1 float)
// ─────────────────────────────────────────────────────────────────────────────
kernel void noncache_postprocess(
    constant NoncacheArgs&   args         [[ buffer(0) ]],
    device const float4*     lig_atoms    [[ buffer(1) ]],
    device float4*           out          [[ buffer(2) ]],
    constant float&          v            [[ buffer(3) ]],
    device atomic_float*     total_energy [[ buffer(4) ]],
    uint gid [[ thread_position_in_grid ]])
{
    if (gid >= args.num_movable_atoms) return;

    float4 val    = out[gid];
    float3 force  = val.xyz;
    float  energy = val.w;

    // Curl correction.
    curl_msl(energy, force, v);

    // Out-of-bounds penalty.
    float3 xyz   = lig_atoms[gid].xyz;
    float3 gb    = float3(args.gbx, args.gby, args.gbz);
    float3 ge    = float3(args.gex, args.gey, args.gez);
    float3 oob_d = float3(0.0f);
    float  oob_p = 0.0f;
    for (int i = 0; i < 3; i++) {
        if (xyz[i] < gb[i]) {
            oob_d[i] = -1.0f; oob_p += gb[i] - xyz[i];
        } else if (xyz[i] > ge[i]) {
            oob_d[i] =  1.0f; oob_p += xyz[i] - ge[i];
        }
        oob_d[i] *= args.slope;
    }
    oob_p  *= args.slope;
    energy += oob_p;
    force  += oob_d;

    out[gid] = float4(force, energy);
    atomic_fetch_add_explicit(total_energy, energy, memory_order_relaxed);
}

// ── Phase 5: tree derivative ──────────────────────────────────────────────────

// Compact node data for derivative pass.
// Matches segment_node_deriv in tree_gpu.h (48 bytes, no doubles, no size_t).
struct NodeDeriv {
    float origin[4];  // world-space origin (xyz + 0)
    float axis[4];    // world-space rotation axis (xyz + 0)
    int   parent;     // BFS index of parent node (-1 for roots)
    int   layer;      // BFS layer index (0 = root)
    int   _pad[2];
};

// Force/torque pair. Matches gfloat4p = gpair<gfloat4,gfloat4> (32 bytes).
struct ForceTorquePair {
    float4 force;
    float4 torque;
};

// Scalar arguments shared by all tree derivative kernels.
// Matches TreeDerivArgs_C in tree_gpu.cu.
struct TreeDerivArgs {
    uint num_atoms;
    uint num_nodes;
    uint num_layers;
    uint nlig_roots;
};

// ─────────────────────────────────────────────────────────────────────────────
// tree_accum_forces  (Phase 5, step 1)
//
// One thread per atom. Atomically accumulates force and torque from each atom
// into the force_torques entry for that atom's owning node.
// force_torques[] must be pre-zeroed by the caller.
//
// Buffer layout:
//   0  constant TreeDerivArgs    — scalar args
//   1  device  float4[]          — coords (world-space atom positions)
//   2  device  float4[]          — forces (per-atom force, i.e. minus_forces)
//   3  device  uint[]            — owners[atom] = owning node BFS index
//   4  device  NodeDeriv[]       — compact node data
//   5  device  ForceTorquePair[] — force_torques (accumulated in-place)
// ─────────────────────────────────────────────────────────────────────────────
kernel void tree_accum_forces(
    constant TreeDerivArgs&    args   [[ buffer(0) ]],
    device const float4*       coords [[ buffer(1) ]],
    device const float4*       forces [[ buffer(2) ]],
    device const uint*         owners [[ buffer(3) ]],
    device const NodeDeriv*    nodes  [[ buffer(4) ]],
    device ForceTorquePair*    ft     [[ buffer(5) ]],
    uint gid [[ thread_position_in_grid ]])
{
    if (gid >= args.num_atoms) return;

    uint   nid    = owners[gid];
    float3 origin = float3(nodes[nid].origin[0],
                           nodes[nid].origin[1],
                           nodes[nid].origin[2]);
    float3 f      = forces[gid].xyz;
    float3 r      = coords[gid].xyz - origin;
    float3 torque = cross(r, f);

    device atomic_float* af = (device atomic_float*)(&ft[nid].force);
    atomic_fetch_add_explicit(af,     f.x,      memory_order_relaxed);
    atomic_fetch_add_explicit(af + 1, f.y,      memory_order_relaxed);
    atomic_fetch_add_explicit(af + 2, f.z,      memory_order_relaxed);

    device atomic_float* at = (device atomic_float*)(&ft[nid].torque);
    atomic_fetch_add_explicit(at,     torque.x, memory_order_relaxed);
    atomic_fetch_add_explicit(at + 1, torque.y, memory_order_relaxed);
    atomic_fetch_add_explicit(at + 2, torque.z, memory_order_relaxed);
}

// ─────────────────────────────────────────────────────────────────────────────
// tree_propagate_layer  (Phase 5, step 2 — one dispatch per layer)
//
// One thread per node. For nodes at the target layer, adds their force_torque
// contribution to their parent's entry (with atomic float add).
// Dispatched from deepest layer down to layer 1.
//
// Buffer layout:
//   0  constant TreeDerivArgs    — scalar args
//   1  device  NodeDeriv[]       — compact node data
//   2  device  ForceTorquePair[] — force_torques (updated in-place)
//   3  constant uint&            — target_layer
// ─────────────────────────────────────────────────────────────────────────────
kernel void tree_propagate_layer(
    constant TreeDerivArgs&   args         [[ buffer(0) ]],
    device const NodeDeriv*   nodes        [[ buffer(1) ]],
    device ForceTorquePair*   ft           [[ buffer(2) ]],
    constant uint&            target_layer [[ buffer(3) ]],
    uint gid [[ thread_position_in_grid ]])
{
    if (gid >= args.num_nodes) return;
    if ((uint)nodes[gid].layer != target_layer) return;

    int    pid      = nodes[gid].parent;
    float3 my_orig  = float3(nodes[gid].origin[0],
                             nodes[gid].origin[1],
                             nodes[gid].origin[2]);
    float3 par_orig = float3(nodes[pid].origin[0],
                             nodes[pid].origin[1],
                             nodes[pid].origin[2]);

    float3 f     = ft[gid].force.xyz;
    float3 t     = ft[gid].torque.xyz;
    float3 r     = my_orig - par_orig;
    float3 add_t = cross(r, f) + t;

    device atomic_float* pf = (device atomic_float*)(&ft[pid].force);
    atomic_fetch_add_explicit(pf,     f.x,     memory_order_relaxed);
    atomic_fetch_add_explicit(pf + 1, f.y,     memory_order_relaxed);
    atomic_fetch_add_explicit(pf + 2, f.z,     memory_order_relaxed);

    device atomic_float* pt = (device atomic_float*)(&ft[pid].torque);
    atomic_fetch_add_explicit(pt,     add_t.x, memory_order_relaxed);
    atomic_fetch_add_explicit(pt + 1, add_t.y, memory_order_relaxed);
    atomic_fetch_add_explicit(pt + 2, add_t.z, memory_order_relaxed);
}

// ─────────────────────────────────────────────────────────────────────────────
// tree_write_change  (Phase 5, step 3)
//
// One thread per node. Projects the accumulated force_torques onto the
// change_gpu DOF layout:
//   Lig-root nodes (gid < nlig_roots): write full 6-DOF (force + torque).
//   All other nodes:                   write scalar torque projected onto axis.
//
// Buffer layout:
//   0  constant TreeDerivArgs      — scalar args
//   1  device  NodeDeriv[]         — compact node data (for axis vectors)
//   2  device  ForceTorquePair[]   — force_torques (read-only)
//   3  device  float[]             — change_gpu.values (written here)
// ─────────────────────────────────────────────────────────────────────────────
kernel void tree_write_change(
    constant TreeDerivArgs&       args  [[ buffer(0) ]],
    device const NodeDeriv*       nodes [[ buffer(1) ]],
    device const ForceTorquePair* ft    [[ buffer(2) ]],
    device float*                 cvals [[ buffer(3) ]],
    uint gid [[ thread_position_in_grid ]])
{
    if (gid >= args.num_nodes) return;

    ForceTorquePair ftp = ft[gid];

    if (gid < args.nlig_roots) {
        // Ligand rigid-body root: full 6-DOF force + torque.
        cvals[6u * gid + 0u] = ftp.force.x;
        cvals[6u * gid + 1u] = ftp.force.y;
        cvals[6u * gid + 2u] = ftp.force.z;
        cvals[6u * gid + 3u] = ftp.torque.x;
        cvals[6u * gid + 4u] = ftp.torque.y;
        cvals[6u * gid + 5u] = ftp.torque.z;
    } else {
        // Residue root or torsion node: scalar torque projected onto axis.
        float3 axis = float3(nodes[gid].axis[0],
                             nodes[gid].axis[1],
                             nodes[gid].axis[2]);
        cvals[gid + 5u * args.nlig_roots] = dot(ftp.torque.xyz, axis);
    }
}
