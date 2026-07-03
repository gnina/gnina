#ifndef VINA_TREE_GPU_H
#define VINA_TREE_GPU_H

#include "tree.h"
#include "atom.h"
#include <vector>
#include <algorithm>
#include "gpu_util.h"
#include "conf_gpu.h"
#include "gpu_math.h"

struct gfloat4 : float4 {
    gfloat4() = default;
    __host__ __device__ gfloat4(float x, float y, float z, float w)
        : float4 { x, y, z, w } {
    }
    __host__ __device__ gfloat4(const vec &v)
        : float4 { v[0], v[1], v[2], 0 } {
    }
    __host__ __device__ gfloat4(const gfloat3 &v)
        : float4 { v[0], v[1], v[2], 0 } {
    }

    __host__ __device__ operator vec() const {
      return vec(x, y, z);
    }

    __host__ __device__
    float& operator[](int b) {
      return b == 0 ? x : b == 1 ? y : b == 2 ? z : w;
    }
    ;

    __host__ __device__
    const float& operator[](int b) const {
      return b == 0 ? x : b == 1 ? y : b == 2 ? z : w;
    }
    ;
};

typedef gpair<gfloat4, gfloat4> gfloat4p;

struct gpu_mat {
    gfloat4 vecs[3];

    gpu_mat() = default;
    __host__ __device__
    gpu_mat(const mat &m);
};

struct __align__(sizeof(float4)) marked_coord {
    gfloat3 coords;
    uint owner_idx;
};

struct ligand;
struct residue;

struct segment_node {
    //a segment is a rigid collection of atoms with an orientation
    //from an axis in a torsion tree
    gfloat4 relative_axis;
    gfloat4 relative_origin;
    gfloat4 axis;
    gfloat4 origin;
    qt orientation_q;
    gpu_mat orientation_m;
    sz begin;
    sz end;
    int parent; //location of parent in node array, -1 if root 
    int layer; //layer in BFS tree

    segment_node()
        : parent(-1), layer(0), begin(0), end(0) {
    }

    segment_node(const segment& node, int p, int layer);

    segment_node(const rigid_body& node);

    segment_node(const first_segment& node);

    __device__ gfloat4p sum_force_and_torque(const gfloat4 *coords,
        const gfloat4 *forces) const;

    __device__ gfloat4 local_to_lab_direction(
        const gfloat4& local_direction) const;

    template<typename VecT>
    __device__ VecT local_to_lab(const VecT& local_coords) const;

    __device__
    void set_coords(const gfloat4 *atom_coords, gfloat4 *coords) const;

    __device__
    void set_orientation(const qt& q);

    __device__
    void set_orientation(float x, float y, float z, float w);

};

//just a branch with parent info convenient for creating a bfs-ordered nodes
//array for the tree_gpu struct
struct pinfo_branch {
    branch* bptr;
    int parent;
    int layer;

    pinfo_branch(branch* b, const int pidx, const int lidx)
        : bptr(b), parent(pidx), layer(lidx) {
    }
};

struct __align__(sizeof(uint2)) atom_node_indices {
    uint atom_idx;
    uint node_idx;

    __host__ __device__ atom_node_indices(void)
        : atom_idx(0), node_idx(0) {
    }
    __host__ __device__ atom_node_indices(uint ai, uint ni)
        : atom_idx(ai), node_idx(ni) {
    }
    ;
};

#ifdef USE_METAL
// Compact, Metal-shader-compatible node data for the derivative pass.
// segment_node contains doubles (qt) and size_t — Metal doesn't support these.
// We extract only the fields needed by the three tree derivative shaders.
// Must match NodeDeriv in gnina_kernels.metal (48 bytes).
struct segment_node_deriv {
    float origin[4];  // world-space origin (xyz + padding)
    float axis[4];    // world-space rotation axis (xyz + padding)
    int   parent;     // BFS index of parent node (-1 for roots)
    int   layer;      // BFS layer index (0 = root layer)
    int   _pad[2];
};
#endif // USE_METAL

struct tree_gpu {

    segment_node *device_nodes;
    gfloat4p *force_torques;
    unsigned num_nodes;
    unsigned num_layers;
    unsigned num_atoms;
    unsigned nlig_roots;
    unsigned nres_roots;

    uint *owners;

    tree_gpu() = default;
    tree_gpu(vector_mutable<ligand> &ligands, vector_mutable<residue> &residues,
        vec *atom_coords, size_t*& dfs_order_bfs_indices,
        size_t*& bfs_order_dfs_indices);

    static void deallocate(tree_gpu *t);

    __device__
    void derivative(const vec *coords, const vec* forces, change_gpu *c);
    __device__
    void set_conf(const vec *atom_coords, vec *coords, const conf_gpu *c);

#ifdef USE_METAL
    // CPU serial implementations for Apple Silicon Metal port (Phase 5).
    // All device_nodes / force_torques / owners arrays live in unified memory
    // (MTLStorageModeShared) and are directly readable/writable from the CPU.
    void set_conf_cpu(const vec *atom_coords, vec *coords, const conf_gpu *c);
    void derivative_cpu(const vec *coords, const vec *forces, change_gpu *c);

    // Metal GPU derivative (replaces derivative_cpu when GPU kernels are available).
    void derivative_metal(const vec *coords, const vec *forces, change_gpu *c);
#endif

  private:
    void do_dfs(branch& branch, size_t& dfs_idx);

    __device__
    void _derivative(const gfloat4 *coords, const gfloat4* forces,
        change_gpu *c);
    __device__
    void _set_conf(const marked_coord *atom_coords, gfloat4 *coords,
        const conf_gpu *c);
};

#endif
