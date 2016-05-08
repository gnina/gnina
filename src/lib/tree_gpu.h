#ifndef VINA_TREE_GPU_H
#define VINA_TREE_GPU_H

#include "tree.h"
#include "atom.h"
#include <vector>
#include <queue>

struct segment_node {
	segment s;
	int child_range[2];
	segment_node() {}
	segment_node(const tree<segment>& subtree, int& index) {
		s = subtree.node;
		child_range[0] = index;
		child_range[1] = index + subtree.children.size()-1;
		index = child_range[1] + 1;
	}
};

struct rigid_body_node {
	rigid_body rb;
	int child_range[2];
	rigid_body_node() {}
	rigid_body_node(const heterotree<rigid_body>& root) {
		rb = root.node;
		child_range[0] = 0;
		child_range[1] = root.children.size()-1;
	}
};

struct tree_gpu {
	rigid_body_node root;
	gvector<segment_node> nodes;
    vecp force_torques[32];

	tree_gpu() {}	
	tree_gpu(const heterotree<rigid_body> &ligand) : root(ligand) {
		// Does BFS starting from the root to construct a 
		// flat tree representation in which each node contains
		// the relevant information for its segment as well as
		// the starting and ending array indices associated
		// with its children
		int index = root.child_range[1] + 1;
		std::queue<const tree<segment> *> tree_traversal;
		for (int i=0; i<index; i++) {
			tree_traversal.push(&ligand.children[i]);
		}
		while (!tree_traversal.empty()) {
			const tree<segment> *next_node = tree_traversal.front();
			tree_traversal.pop();
			nodes.push_back(segment_node(*next_node,index));
			int num_children = next_node->children.size();
			for (int i=0; i<num_children; i++) {
				tree_traversal.push(&next_node->children[i]);
			}
		}
	}	

    __device__ __host__    
	void do_derivatives(const segment_node& node, fl *&p,
                        vecp *force_torques, int index, const gvecv& coords,
                        const gvecv& forces) {

		vecp tmp = node.s.sum_force_and_torque(coords, forces);

		fl& d = *p;
		++p;
		for (int i=node.child_range[0]; i<=node.child_range[1]; i++) {
			do_derivatives(nodes[i], p, force_torques, i, coords, forces);
			vecp child_result = force_torques[i];
			tmp.first += child_result.first;
			vec r = nodes[i].s.get_origin() - node.s.get_origin();
			tmp.second += cross_product(r, child_result.first) + 
				child_result.second;
		}
		force_torques[index].first = tmp.first;
		force_torques[index].second = tmp.second;

		d = force_torques[index].second * node.s.axis;	
	}

    __device__ __host__
	void derivative(const gvecv& coords, const gvecv& forces, ligand_change& c) {
        /* TODO: NODES SIZE. Gotta configure with kern launch */
        
		// Don't put the root's values in the force_torque array
		vecp force_torque = root.rb.sum_force_and_torque(coords, forces);	
		fl *p = &c.torsions[0];
		for (int i=0; i<=root.child_range[1]; i++) {
			// Recursively calculate derivatives for the children using
			// DFS in order to properly update the torsions array
			do_derivatives(nodes[i], p, force_torques, i, coords, forces);
			vec r = nodes[i].s.get_origin() - root.rb.get_origin();
			vecp child_result = force_torques[i];
			force_torque.first += child_result.first;
			force_torque.second += cross_product(r, child_result.first) + 
				child_result.second;
		}
		c.rigid.position = force_torque.first;
		c.rigid.orientation = force_torque.second;
	}

    __device__ __host__    
	void do_confs(segment_node& node, int my_index,
                  const frame& parent, const gatomv& atoms, 
                  gvecv& coords, const fl *&c) {
		node.s.set_conf(parent, atoms, coords, c);
		for (int i=node.child_range[0]; i<=node.child_range[1]; i++) {
			do_confs(nodes[i], i, nodes[my_index].s, atoms, coords, c);
		}
	}

    __device__ __host__
	void set_conf(const gatomv& atoms, gvecv& coords, const ligand_conf& c) {
		root.rb.set_conf(atoms, coords, c.rigid);
		const fl* p = &c.torsions[0];
		for (int i=0; i<=root.child_range[1]; i++) {
			do_confs(nodes[i], i, root.rb, atoms, coords, p);
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
                     gvecv& coords, const ligand_conf& c) {
    t.set_conf(atoms, coords, c);
}

#endif
