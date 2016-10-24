#include "tree_gpu.h"
#include <algorithm>

__device__
vecp segment_node::sum_force_and_torque(const vec *coords, const vec *forces) const {
	vecp tmp(vec(0,0,0), vec(0,0,0));
	VINA_RANGE(i, begin, end) {
		tmp.first  += forces[i];
		tmp.second += cross_product(coords[i] - origin, forces[i]);
	}
	return tmp;
}

__device__
vec segment_node::local_to_lab_direction(const vec& local_direction) const{
	vec tmp;
	tmp = orientation_m * local_direction;
	return tmp;
}

__device__
vec segment_node::local_to_lab(const vec& local_coords) const{
	vec tmp;
	tmp = origin + orientation_m * local_coords;
	return tmp;
}

__device__
void segment_node::set_coords(const vec *atom_coords, vec *coords) const{
	VINA_RANGE(i, begin, end)
		coords[i] = local_to_lab(atom_coords[i]);
}

__device__
void segment_node::set_orientation(const qt& q) { // does not normalize the orientation
	orientation_q = q;
	orientation_m = quaternion_to_r3(orientation_q);
}

__device__
void segment_node::set_orientation(float x, float y, float z, float w) { // does not normalize the orientation
	orientation_q = qt(x,y,z,w);
	orientation_m = quaternion_to_r3(orientation_q);
}

void tree_gpu::do_dfs(int parent, const branch& branch,
        std::vector<segment_node>& nodes, std::vector<unsigned>&
        atoms_per_layer_host, std::vector<atom_node_indices>&
        atom_node_prelist) {
	segment_node node(branch.node, parent, &nodes[parent]);
	unsigned index = nodes.size();
	nodes.push_back(node);

    if (atoms_per_layer_host.size() - 1 < node.layer) {
        atoms_per_layer_host.push_back(node.end - node.begin);
    }
    else {
        atoms_per_layer_host[node.layer] += node.end - node.begin;
    }

    if (node.layer > num_layers)
        num_layers = node.layer;

    for (unsigned i=node.begin; i<node.end; i++) {
        atom_node_prelist.push_back(atom_node_indices(i,nodes.size()-1));
    }

	VINA_FOR_IN(i, branch.children) {
		do_dfs(index, branch.children[i], nodes, atoms_per_layer_host,
                atom_node_prelist);
	}
}

tree_gpu::tree_gpu(const heterotree<rigid_body> &ligand){
	//populate nodes in DFS order from ligand, where node zero is the root
	std::vector<segment_node> nodes;
    std::vector<unsigned> atoms_per_layer_host;
    std::vector<atom_node_indices> atom_node_prelist;
    num_layers = 0;

	segment_node root(ligand.node);
	nodes.push_back(root);
    atoms_per_layer_host.push_back(root.end - root.begin);
    for (unsigned i=root.begin; i<root.end; i++) {
        atom_node_prelist.push_back(atom_node_indices(i,nodes.size()-1));
    }

	VINA_FOR_IN(i, ligand.children) {
		do_dfs(0,ligand.children[i], nodes, atoms_per_layer_host, atom_node_prelist);
	}

    max_atoms_per_layer = *std::max_element(std::begin(atoms_per_layer_host),
            std::end(atoms_per_layer_host));
	num_nodes = nodes.size();
	//allocate device memory and copy
	//nodes
	cudaMalloc(&device_nodes, sizeof(segment_node)*nodes.size());
	cudaMemcpy(device_nodes, &nodes[0], sizeof(segment_node)*nodes.size(), cudaMemcpyHostToDevice);

	//forcetorques
	cudaMalloc(&force_torques, sizeof(vecp)*nodes.size());
	cudaMemset(force_torques, 0, sizeof(vecp)*nodes.size());

    //atom and node values
    cudaMalloc(&atom_node_list,
            sizeof(atom_node_indices)*atom_node_prelist.size());
    cudaMemcpy(atom_node_list, &atom_node_prelist[0],
            sizeof(atom_node_indices)*atom_node_prelist.size(),
            cudaMemcpyHostToDevice);

    //atoms per layer
    cudaMalloc(&atoms_per_layer, sizeof(unsigned)*atoms_per_layer_host.size());
    cudaMemcpy(atoms_per_layer, &atoms_per_layer_host[0],
            sizeof(unsigned)*atoms_per_layer_host.size(),
            cudaMemcpyHostToDevice);
}

//given a gpu point, deallocate all the memory
void tree_gpu::deallocate(tree_gpu *t) {
	tree_gpu cpu;
	cudaMemcpy(&cpu, t, sizeof(tree_gpu), cudaMemcpyDeviceToHost);
	cudaFree(cpu.device_nodes);
	cudaFree(cpu.force_torques);
	cudaFree(t);
}

__device__
void tree_gpu::derivative(const vec *coords,const vec* forces, float *c){

	//calculate each segments individual force/torque
    int index = threadIdx.x;
    force_torques[index] = device_nodes[index].sum_force_and_torque(coords,
            forces);

    unsigned current_layer = num_layers - 1;
	//have each child add its contribution to its parents force_torque
    const segment_node& cnode = device_nodes[index];
    __syncthreads();
    while (current_layer > 0) {
        if (cnode.layer == current_layer) {
            unsigned parent = device_nodes[index].parent;
            const segment_node& pnode = device_nodes[parent];
            const vecp& ft = force_torques[index];
            pseudoAtomicAdd(&force_torques[parent].first, ft.first);
            
            vec r = cnode.origin - pnode.origin;
            pseudoAtomicAdd(&force_torques[parent].second, cross_product(r,
                        ft.first)+ft.second);
            
            //set torsions
            c[6+index-1] = ft.second * cnode.axis;
        }
        current_layer--;
        __syncthreads();
    }

    if (index == 0) {
	    c[0] = force_torques[0].first[0];
	    c[1] = force_torques[0].first[1];
	    c[2] = force_torques[0].first[2];

	    c[3] = force_torques[0].second[0];
	    c[4] = force_torques[0].second[1];
	    c[5] = force_torques[0].second[2];
    }
}

__device__
void tree_gpu::set_conf(const vec *atom_coords, vec *coords, const conf_info
		*c, unsigned nlig_atoms){
	// assert(c.torsions.size() == num_nodes-1);
	int index = threadIdx.x;
    segment_node* node = NULL;
	uint natoms;
	uint current_layer = 0;

	if (index < atoms_per_layer[current_layer]) {
        node = &device_nodes[current_layer];
		for(unsigned i = 0; i < 3; i++)
			node->origin[i] = c->position[i];
		node->set_orientation(c->orientation[0],c->orientation[1],c->orientation[2],c->orientation[3]);
        coords[index] = node->local_to_lab(atom_coords[index]);
	}

    natoms = atoms_per_layer[current_layer];
	__syncthreads();
	while (current_layer < num_layers) {
	    current_layer++;

		if (index < atoms_per_layer[current_layer]) {
            atom_node_indices idx_pair = atom_node_list[natoms + index];
            node = &device_nodes[idx_pair.node_idx];
			segment_node& parent = device_nodes[node->parent];
			fl torsion = c->torsions[idx_pair.node_idx-1];
			node->origin = parent.local_to_lab(node->relative_origin);
			node->axis = parent.local_to_lab_direction(node->relative_axis);
			node->set_orientation(
					quaternion_normalize_approx(
							angle_to_quaternion(node->axis, torsion) * parent.orientation_q));
            coords[idx_pair.atom_idx] =
                node->local_to_lab(atom_coords[idx_pair.atom_idx]);
		}

        natoms += atoms_per_layer[current_layer];
        __syncthreads();
	}
}
