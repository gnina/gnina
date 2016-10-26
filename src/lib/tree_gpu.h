#ifndef VINA_TREE_GPU_H
#define VINA_TREE_GPU_H

#include "tree.h"
#include "atom.h"
#include <vector>
#include <queue>
#include "gpu_util.h"
#include "conf_gpu.h"
#include "gpu_math.h"

struct segment_node {
	//a segment is a rigid collection of atoms with an orientation
	//from an axis in a torsion tree
	vec relative_axis;
	vec relative_origin;
	vec axis;
	vec origin;
	qt orientation_q;
	mat orientation_m;
	sz begin;
	sz end;
	int parent; //location of parent in node array, -1 if root
	int layer; //layer in BFS tree

	segment_node() :
			parent(-1), layer(0), begin(0), end(0){
	}

	segment_node(const segment& node,int p,segment_node* pnode) :
			relative_axis(node.relative_axis), relative_origin(node.relative_origin), axis(
					node.axis), origin(node.origin), orientation_q(node.orientation_q), orientation_m(
					node.orientation_m), begin(node.begin), end(node.end), parent(p){
		layer = pnode->layer + 1;

	}

	segment_node(const rigid_body& node) : //root node
			relative_axis(0, 0, 0), relative_origin(0, 0, 0), axis(0, 0, 0), origin(
					node.origin), orientation_q(node.orientation_q), orientation_m(
					node.orientation_m), begin(node.begin), end(node.end), parent(-1), layer(0){

	}

	__device__
	vecp sum_force_and_torque(const vec *coords, const vec *forces) const;

	__device__
	vec local_to_lab_direction(const vec& local_direction) const;

	__device__
	vec local_to_lab(const vec& local_coords) const;

	__device__
	void set_coords(const vec *atom_coords, vec *coords) const;

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
	vecp *force_torques;
	unsigned num_nodes;
    unsigned num_layers;
    atom_node_indices *atom_node_list; // atom and corresponding node indices in bfs order
    unsigned *atoms_per_layer;
    unsigned max_atoms_per_layer; 

	tree_gpu() :
			device_nodes(NULL), force_torques(NULL), num_nodes(0), num_layers(0), atom_node_list(NULL), atoms_per_layer(NULL), max_atoms_per_layer(0) {
	}

	void do_dfs(int parent, const branch& branch, std::vector<segment_node>& nodes, std::vector<unsigned>& atoms_per_layer_host, std::map<int,std::vector<atom_node_indices>>& atom_node_map);

	tree_gpu(const heterotree<rigid_body> &ligand);

	static void deallocate(tree_gpu *t);

	__device__
	void derivative(const vec *coords,const vec* forces, float *c);

	__device__
	void set_conf(const vec *atom_coords, vec *coords, const conf_info *c, unsigned nlig_atoms);

};

#endif
