#ifndef VINA_TREE_GPU_H
#define VINA_TREE_GPU_H

#include "tree.h"
#include "atom.h"
#include <vector>
#include <queue>

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

	segment_node() :
			parent(-1), begin(0), end(0){
	}

	segment_node(const segment& node,int p) :
			relative_axis(node.relative_axis), relative_origin(node.relative_origin), axis(
					node.axis), origin(node.origin), orientation_q(node.orientation_q), orientation_m(
					node.orientation_m), begin(node.begin), end(node.end), parent(p){

	}

	segment_node(const rigid_body& node) : //root node
			relative_axis(0, 0, 0), relative_origin(0, 0, 0), axis(0, 0, 0), origin(
					node.origin), orientation_q(node.orientation_q), orientation_m(
					node.orientation_m), begin(node.begin), end(node.end), parent(-1){

	}

	__device__
	vecp sum_force_and_torque(const gvecv& coords, const gvecv& forces) const {
		vecp tmp(vec(0,0,0), vec(0,0,0));
		VINA_RANGE(i, begin, end) {
			tmp.first  += forces[i];
			tmp.second += cross_product(coords[i] - origin, forces[i]);
		}
		return tmp;
	}

	__device__
	vec local_to_lab_direction(const vec& local_direction) const{
		vec tmp;
		tmp = orientation_m * local_direction;
		return tmp;
	}

	__device__
	vec local_to_lab(const vec& local_coords) const{
		vec tmp;
		tmp = origin + orientation_m * local_coords;
		return tmp;
	}

	__device__
	void set_coords(const gatomv& atoms, gvecv& coords) const{
		VINA_RANGE(i, begin, end)
			coords[i] = local_to_lab(atoms[i].coords);
	}

	__device__
	void set_orientation(const qt& q) { // does not normalize the orientation
		orientation_q = q;
		orientation_m = quaternion_to_r3(orientation_q);
	}
};

struct tree_gpu {

	segment_node *device_nodes;
	vecp *force_torques;
	unsigned num_nodes;

	tree_gpu() :
			device_nodes(NULL), force_torques(NULL), num_nodes(0) {
	}

	void do_dfs(int parent, const branch& branch, std::vector<segment_node>& nodes) {
		segment_node node(branch.node, parent);
		unsigned index = nodes.size();
		nodes.push_back(node);

		VINA_FOR_IN(i, branch.children) {
			do_dfs(index, branch.children[i], nodes);
		}
	}

	tree_gpu(const heterotree<rigid_body> &ligand){
		//populate nodes in DFS order from ligand, where node zero is the root
		std::vector<segment_node> nodes;
		segment_node root(ligand.node);
		nodes.push_back(root);

		VINA_FOR_IN(i, ligand.children) {
			do_dfs(0,ligand.children[i], nodes);
		}

		num_nodes = nodes.size();
		//allocate device memory and copy
		//nodes
		cudaMalloc(&device_nodes, sizeof(segment_node)*nodes.size());
		cudaMemcpy(device_nodes, &nodes[0], sizeof(segment_node)*nodes.size(), cudaMemcpyHostToDevice);

		//forcetorques
		cudaMalloc(&force_torques, sizeof(vecp)*nodes.size());
		cudaMemset(force_torques, 0, sizeof(vecp)*nodes.size());

	}

	__device__
	void derivative(const gvecv& coords,const gvecv& forces,ligand_change& c){

		assert(c.torsions.size() == num_nodes-1);
		//calculate each segments individual force/torque
		for(unsigned i = 0; i < num_nodes; i++) {
			force_torques[i] = device_nodes[i].sum_force_and_torque(coords, forces);
		}

		//have each child add its contribution to its parents force_torque
		for(unsigned i = num_nodes-1; i > 0; i--) {
			unsigned parent = device_nodes[i].parent;
			const vecp& ft = force_torques[i];
			force_torques[parent].first += ft.first;

			const segment_node& pnode = device_nodes[parent];
			const segment_node& cnode = device_nodes[i];

			vec r = cnode.origin - pnode.origin;
			force_torques[parent].second += cross_product(r, ft.first)+ft.second;

			//set torsions
			c.torsions[i-1] = ft.second * cnode.axis;
		}

		c.rigid.position = force_torques[0].first;
		c.rigid.orientation = force_torques[0].second;
	}

	__device__
	void set_conf(const gatomv& atoms, gvecv& coords, const ligand_conf& c){
		assert(c.torsions.size() == num_nodes-1);

		segment_node& root = device_nodes[0];
		root.origin = c.rigid.position;
		root.set_orientation(c.rigid.orientation);
		root.set_coords(atoms, coords);

		for(unsigned i = 1; i < num_nodes; i++) {
			segment_node& node = device_nodes[i];
			segment_node& parent = device_nodes[node.parent];
			fl torsion = c.torsions[i-1];
			node.origin = parent.local_to_lab(node.relative_origin);
			node.axis = parent.local_to_lab_direction(node.relative_axis);
			node.set_orientation(
					quaternion_normalize_approx(
							angle_to_quaternion(node.axis, torsion) * parent.orientation_q));
			node.set_coords(atoms, coords);
		}
	}
};

static __global__
void derivatives_kernel(tree_gpu &t, const gvecv& coords,
		const gvecv& forces, ligand_change& c){

	t.derivative(coords, forces, c);
}

static __global__
void set_conf_kernel(tree_gpu &t, const gatomv& atoms,
		gvecv& coords, const ligand_conf& c){
	t.set_conf(atoms, coords, c);
}

#endif
