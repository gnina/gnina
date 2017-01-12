#ifndef VINA_TREE_GPU_H
#define VINA_TREE_GPU_H

#include "tree.h"
#include "atom.h"
#include <vector>
#include <queue>
#include "gpu_util.h"
#include "conf_gpu.h"
#include "gpu_math.h"

struct gfloat4 : float4{
    __host__ __device__
    gfloat4() = default;
    __host__ __device__
    gfloat4(float x, float y, float z, float w) : float4{x,y,z,w}{}
    __host__ __device__
    gfloat4(const vec &v) : float4{v[0], v[1], v[2], 0}{}
    __host__ __device__
    gfloat4(const float3 &v) : float4{v[0], v[1], v[2], 0}{}

    __host__ __device__
    operator vec() const{
        return vec(x, y, z);
    }
    
    __host__ __device__
    float& operator[](int b){
        return b == 0 ? x :
               b == 1 ? y :
               b == 2 ? z :
               w;
    };
    
    __host__ __device__
    const float& operator[](int b) const{
        return b == 0 ? x :
               b == 1 ? y :
               b == 2 ? z :
               w;
    };
};

typedef gpair<gfloat4, gfloat4> gfloat4p;

struct gpu_mat{
    gfloat4 vecs[3];

    __host__ __device__
    gpu_mat() = default;
    __host__ __device__
    gpu_mat(const mat &m);
};

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

	segment_node() :
			parent(-1), layer(0), begin(0), end(0){
	}

	segment_node(const segment& node,int p,segment_node* pnode);
	segment_node(const rigid_body& node);

	__device__
	gfloat4p sum_force_and_torque(const gfloat4 *coords, const gfloat4 *forces) const;

	__device__
	gfloat4 local_to_lab_direction(const gfloat4& local_direction) const;

    template <typename VecT>
    __device__
    VecT local_to_lab(const VecT& local_coords) const;

	__device__
	void set_coords(const gfloat4 *atom_coords, gfloat4 *coords) const;

	__device__
	void set_orientation(const qt& q);

	__device__
	void set_orientation(float x, float y, float z, float w); 

};

struct __align__(sizeof(uint2)) atom_node_indices{
    uint atom_idx;
    uint node_idx;

    __host__ __device__ atom_node_indices(void): atom_idx(0), node_idx(0) {}
    __host__ __device__ atom_node_indices(uint ai, uint ni) : atom_idx(ai), node_idx(ni) {};
};

struct tree_gpu {

	segment_node *device_nodes;
	gfloat4p *force_torques;
	unsigned num_nodes;
    unsigned num_layers;

    uint *owners;

    tree_gpu() = default;
	tree_gpu(const heterotree<rigid_body> &ligand, vec *atom_coords);

	static void deallocate(tree_gpu *t);

    __device__
	void derivative(const vec *coords,const vec* forces, float *c);
    __device__
	void set_conf(const vec *atom_coords, vec *coords,
                  const conf_info *c);

    struct __align__(sizeof(float4)) marked_coord{
        float3 coords;
        uint owner_idx;
    };

private:
    void do_dfs(int parent, const branch& branch, 
                std::vector<segment_node> &nodes);

	__device__
	void _derivative(const gfloat4 *coords,const gfloat4* forces, float *c);
	__device__
	void _set_conf(const marked_coord *atom_coords, gfloat4 *coords,
                   const conf_info *c);
};

#endif
