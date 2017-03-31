#include "tree_gpu.h"
#include "model.h"
#include <algorithm>
#include <map>

__device__
gfloat4 pseudoAtomicAdd(gfloat4* address, gfloat4 value) {
    return gfloat4(atomicAdd(&((*address)[0]), value[0]),
                   atomicAdd(&((*address)[1]), value[1]),
                   atomicAdd(&((*address)[2]), value[2]),
                   0);
}

__host__ __device__
gpu_mat::gpu_mat(const mat &m) :
    vecs{gfloat4(m.data[0], m.data[1], m.data[2], 0),
         gfloat4(m.data[3], m.data[4], m.data[5], 0),
         gfloat4(m.data[6], m.data[7], m.data[8], 0)}{}

__host__ __device__
gfloat4 operator*(const gpu_mat &m, const gfloat4&_v){
    gfloat4 v = _v;
    return
    gfloat4(m.vecs[0][0] * v[0] + m.vecs[1][0] * v[1] + m.vecs[2][0] * v[2], 
            m.vecs[0][1] * v[0] + m.vecs[1][1] * v[1] + m.vecs[2][1] * v[2],
            m.vecs[0][2] * v[0] + m.vecs[1][2] * v[1] + m.vecs[2][2] * v[2],
            0);
}

__host__ __device__
gfloat4 operator*(const gpu_mat &m, const float3&v){
    return
    gfloat4(m.vecs[0][0] * v[0] + m.vecs[1][0] * v[1] + m.vecs[2][0] * v[2], 
            m.vecs[0][1] * v[0] + m.vecs[1][1] * v[1] + m.vecs[2][1] * v[2],
            m.vecs[0][2] * v[0] + m.vecs[1][2] * v[1] + m.vecs[2][2] * v[2],
            0);
}

inline __host__ __device__ gfloat4 operator-(gfloat4 a, gfloat4 b)
{
    return gfloat4(a.x - b.x, a.y - b.y, a.z - b.z,  0);
}

inline __host__ __device__ gfloat4 operator+(gfloat4 a, gfloat4 b)
{
    return gfloat4(a.x + b.x, a.y + b.y, a.z + b.z,  0);
}

inline __host__ __device__ gfloat4 &operator+=(gfloat4 &a, gfloat4 b)
{
    a = a + b;
    return a;
}

inline __host__ __device__ gfloat4 operator*(gfloat4 a, float s)
{
    return gfloat4(a.x * s, a.y * s, a.z * s,  0);
}

inline __host__ __device__ float operator*(gfloat4 a, gfloat4 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


/* inline __host__ __device__ gfloat4 operator-(gfloat4 a, gfloat4 b) */
/* { */
/*     return gfloat4(a.x - b.x, a.y - b.y, a.z - b.z,  a.w - b.w); */
/* } */

/* inline __host__ __device__ gfloat4 operator+(gfloat4 a, gfloat4 b) */
/* { */
/*     return gfloat4(a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w); */
/* } */

/* inline __host__ __device__ gfloat4 &operator+=(gfloat4 &a, gfloat4 b) */
/* { */
/*     a = a + b; */
/*     return a; */
/* } */

/* inline __host__ __device__ gfloat4 operator*(gfloat4 a, float s) */
/* { */
/*     return gfloat4(a.x * s, a.y * s, a.z * s,  a.w * s); */
/* } */

/* inline __host__ __device__ float operator*(gfloat4 a, gfloat4 b) */
/* { */
/*     return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; */
/* } */

segment_node::segment_node(const segment& node,int p,segment_node* pnode) :
    relative_axis(node.relative_axis),
    relative_origin(node.relative_origin),
    axis(node.axis),
    origin(node.origin),
    orientation_q(node.orientation_q),
    orientation_m(node.orientation_m),
    begin(node.begin),
    end(node.end),
    parent(p),
    layer(pnode->layer + 1),
    root(pnode->root){}

segment_node::segment_node(const first_segment& node, int root) : //root node
    relative_axis(0,0,0,0), relative_origin(0,0,0,0),
    axis(node.axis), origin(node.origin),
    orientation_q(node.orientation_q), orientation_m(node.orientation_m),
    begin(node.begin), end(node.end), parent(-1), layer(0), root(root){

}
segment_node::segment_node(const rigid_body& node, int root) : //root node
    relative_axis(0,0,0,0), relative_origin(0,0,0,0),
    axis(0,0,0,0), origin(node.origin),
    orientation_q(node.orientation_q), orientation_m(node.orientation_m),
    begin(node.begin), end(node.end), parent(-1), layer(0), root(root){

}

__device__
gfloat4 segment_node::local_to_lab_direction(const gfloat4& local_direction) const{
	gfloat4 tmp;
	tmp = orientation_m * local_direction;
	return tmp;
}

template <typename VecT>
__device__
VecT segment_node::local_to_lab(const VecT& local_coords) const{
	VecT tmp;
	tmp = origin + orientation_m * local_coords;
	return tmp;
}

__device__
void segment_node::set_coords(const gfloat4 *atom_coords, gfloat4 *coords) const{
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

template<typename cpu_root>
void tree_gpu<cpu_root>::do_dfs(int parent, const branch& branch,
                      std::vector<segment_node> &nodes){
	segment_node node(branch.node, parent, &nodes[parent]);
	unsigned index = nodes.size();
	nodes.push_back(node);
	VINA_FOR_IN(i, branch.children) {
		do_dfs(index, branch.children[i], nodes);
	}
}
__device__ 
void ligand_root::set_conf(const float* c) {
	for(unsigned i = 0; i < 3; i++)
	    origin[i] = c[i];
	set_orientation(c[3], c[4],
                      c[5], c[6]);
}
__device__ 
void residue_root::set_conf(const float* c) {
	set_orientation(angle_to_quaternion(axis, *c));
}

template<typename cpu_root>
tree_gpu<cpu_root>::tree_gpu(const vector_mutable<cpu_root> &source, vec *atom_coords){
	//populate nodes in DFS order
	std::vector<segment_node> nodes;
    uint cpu_subtree_sizes[source.size()+1];
   
    VINA_FOR_IN(i, source)
	    nodes.push_back(segment_node(source[i].node,i));

    cpu_subtree_sizes[0] = nodes.size();
    for (unsigned i=0; i < source.size(); i++) {
        for (unsigned j=0; j < source[i].children.size(); j++) 
	    	do_dfs(i, source[i].children[j], nodes);
        cpu_subtree_sizes[i+1] = nodes.size() - cpu_subtree_sizes[0];
    }

    uint natoms = 0;
    uint max_layer = 0;
    for(auto n : nodes){
        natoms += n.end - n.begin;
        //TODO: just BFS-order the nodes
        max_layer = std::max((uint) n.layer, max_layer);
    }
    num_atoms = natoms;
    
    marked_coord *atoms = (marked_coord *) atom_coords;
    uint cpu_owners[natoms];
    for(int i = 0; i < nodes.size(); i++)
        for(int ai = nodes[i].begin; ai < nodes[i].end; ai++)
            cpu_owners[ai] = atoms[ai].owner_idx = i;

	num_nodes = nodes.size();
    num_layers = max_layer + 1;
   
    //TODO: except for subtree_sizes, these should only be malloc'ed if
    //nodes.size()
	cudaMalloc(&device_nodes, sizeof(segment_node)*nodes.size());
	cudaMemcpy(device_nodes, &nodes[0], sizeof(segment_node)*nodes.size(), cudaMemcpyHostToDevice);

	cudaMalloc(&force_torques, sizeof(gfloat4p)*nodes.size());
	cudaMemset(force_torques, 0, sizeof(gfloat4p)*nodes.size());

    cudaMalloc(&owners, sizeof(uint) * natoms);
    cudaMemcpy(owners, cpu_owners, sizeof(uint) * natoms, cudaMemcpyHostToDevice);

    cudaMalloc(&subtree_sizes, sizeof(uint) * (source.size()+1));
    cudaMemcpy(subtree_sizes, cpu_subtree_sizes, sizeof(uint) *
            (source.size()+1), cudaMemcpyHostToDevice);

    /* assert(num_nodes == source.count_torsions() + 1); */
    
    for(uint &o : cpu_owners){
        segment_node &n = nodes[o];
        assert(&o >= &cpu_owners[n.begin] &&
               &o < &cpu_owners[n.end]);
        assert(o < nodes.size());
        
    }
}

//given a gpu point, deallocate all the memory
template<typename cpu_root>
void tree_gpu<cpu_root>::deallocate(tree_gpu *t) {
	tree_gpu cpu;
	cudaMemcpy(&cpu, t, sizeof(tree_gpu), cudaMemcpyDeviceToHost);
	cudaFree(cpu.device_nodes);
	cudaFree(cpu.force_torques);
	cudaFree(t);
}

template<typename cpu_root>
__device__
void tree_gpu<cpu_root>::derivative(const vec *coords,const vec* forces, change_gpu *c){
    _derivative((gfloat4 *) coords, (gfloat4 *) forces, c);
}

template<typename cpu_root>
__device__
void tree_gpu<cpu_root>::_derivative(const gfloat4 *coords,const gfloat4* forces, change_gpu *c){
	// assert(c.torsions.size() == num_nodes-1);
	
	//calculate each segments individual force/torque
    uint tid = threadIdx.x;
    uint nid = owners[tid];
    segment_node &owner = device_nodes[nid];

    if(tid < num_nodes)
        force_torques[tid] = gfloat4p(gfloat4(0,0,0,0), gfloat4(0,0,0,0));
    __syncthreads();
    
    pseudoAtomicAdd(&force_torques[nid].first, forces[tid]);
    pseudoAtomicAdd(&force_torques[nid].second,
                    cross_product(coords[tid] - owner.origin, forces[tid]));


	//have each child add its contribution to its parents force_torque
    uint layer = 0;
    uint parent_id;
    gfloat4 axis;
    gfloat4 origin;
    gfloat4 parent_origin;
    gfloat4p ft;
    uint root = 0;
    
    if(tid > (subtree_sizes[0]-1) && tid < num_nodes){
        segment_node &n = device_nodes[tid];
        parent_id = n.parent;
        parent_origin = device_nodes[parent_id].origin;
        layer = n.layer;
        axis = n.axis;
        origin = n.origin;
        root = n.root;
    }
    
    for(uint l = num_layers - 1; l > 0; l--){
        __syncthreads();
        if(l != layer)
            continue;
        ft = force_torques[tid];
        gfloat4 r = origin - parent_origin;
        
        pseudoAtomicAdd(&force_torques[parent_id].first, ft.first);
        pseudoAtomicAdd(&force_torques[parent_id].second,
                        ::operator+(cross_product(r, ft.first), ft.second));
    }
   
    if (gpu_root::in_unison()) {
        unsigned relative_offset = tid < subtree_sizes[0] ? tid +
            subtree_sizes[tid] : tid + root - subtree_sizes[0];
        if (tid == 0)
            relative_offset = 0;
        c->values[c->flex_offset + relative_offset] = ft.second * axis;
    }
    else if(tid > (subtree_sizes[0]-1) && tid < num_nodes)
        c->values[tid - subtree_sizes[0] + (root+1)*6] = ft.second * axis;
    else if (tid < num_nodes) {
        ft = force_torques[tid];
        uint offset = tid == 0 ? 0 : subtree_sizes[tid];
	    c->values[6*tid + offset] = ft.first[0];
	    c->values[6*tid + offset + 1] = ft.first[1];
	    c->values[6*tid + offset + 2] = ft.first[2];

	    c->values[6*tid + offset + 3] = ft.second[0];
	    c->values[6*tid + offset + 4] = ft.second[1];
	    c->values[6*tid + offset + 5] = ft.second[2];
    }
}

__device__
gfloat4p segment_node::sum_force_and_torque(const gfloat4 *coords, const gfloat4 *forces) const {
	gfloat4p tmp(gfloat4(0, 0, 0, 0), gfloat4(0, 0, 0, 0));
	VINA_RANGE(i, begin, end) {
		tmp.first  += forces[i];
		tmp.second += cross_product(coords[i] - origin, forces[i]);
	}
	return tmp;
}

template<typename cpu_root>
__device__
void tree_gpu<cpu_root>::set_conf(const vec *atom_coords, vec *coords,
                        const conf_gpu* c){
    _set_conf((marked_coord *) atom_coords, (gfloat4 *) coords, c);
}


template<typename cpu_root>
__device__
void tree_gpu<cpu_root>::_set_conf(const marked_coord *atom_coords, gfloat4 *coords,
                         const conf_gpu *c){
    
	uint tid = threadIdx.x;
    segment_node &n = device_nodes[tid];
    marked_coord a = atom_coords[tid];

    uint layer = 0;
    segment_node *parent = NULL;
    fl torsion = 0;
    uint rootid = 0;
    
    if(tid < subtree_sizes[0]){
        gpu_root* root = static_cast<gpu_root*>(&n);
        unsigned offset = tid == 0 ? 0 : subtree_sizes[tid+1];
        root->set_conf(&c->values[7*tid + offset]);
    }else if(tid < num_nodes){
        layer = n.layer;
        parent = &device_nodes[n.parent];
        rootid = n.root;
        torsion = c->values[(rootid+1) * gpu_root::get_offset() + tid -
            subtree_sizes[0]];
    }

    for(uint i = 1; i < num_layers; i++){
        __syncthreads();
        if(layer != i)
            continue;
        n.origin = parent->local_to_lab(n.relative_origin);
        n.axis = parent->local_to_lab_direction(n.relative_axis);
        n.set_orientation(
            quaternion_normalize_approx(
                angle_to_quaternion(n.axis, torsion) * parent->orientation_q));
    }
    __syncthreads();

    assert(owners[tid] == a.owner_idx);
    assert(tid >= device_nodes[a.owner_idx].begin &&
           tid < device_nodes[a.owner_idx].end);

    coords[tid] = device_nodes[a.owner_idx].local_to_lab(a.coords);
}

template struct tree_gpu<ligand>;
template struct tree_gpu<residue>;
