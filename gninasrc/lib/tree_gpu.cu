#include "tree_gpu.h"
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
    layer(pnode->layer + 1){}

segment_node::segment_node(const rigid_body& node) : //root node
    relative_axis(0,0,0,0), relative_origin(0,0,0,0),
    axis(0,0,0,0), origin(node.origin),
    orientation_q(node.orientation_q), orientation_m(node.orientation_m),
    begin(node.begin), end(node.end), parent(-1), layer(0){

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

void tree_gpu::do_dfs(int parent, const branch& branch,
                      std::vector<segment_node> &nodes){
	segment_node node(branch.node, parent, &nodes[parent]);
	unsigned index = nodes.size();
	nodes.push_back(node);
	VINA_FOR_IN(i, branch.children) {
		do_dfs(index, branch.children[i], nodes);
	}
}

struct __align__(sizeof(float4)) marked_coord{
    float3 coords;
    uint owner_idx;
};

tree_gpu::tree_gpu(const heterotree<rigid_body> &ligand, vec *atom_coords){
	//populate nodes in DFS order from ligand, where node zero is the root
	std::vector<segment_node> nodes;
    
	nodes.push_back(segment_node(ligand.node));

	VINA_FOR_IN(i, ligand.children) {
		do_dfs(0,ligand.children[i], nodes);
	}

    uint natoms = 0;
    uint max_layer = 0;
    for(auto n : nodes){
        natoms += n.end - n.begin;
        //TODO: just BFS-order the nodes
        max_layer = std::max((uint) n.layer, max_layer);
    }
    
    marked_coord *atoms = (marked_coord *) atom_coords;
    uint cpu_owners[natoms];
    for(int i = 0; i < nodes.size(); i++)
        for(int ai = nodes[i].begin; ai < nodes[i].end; ai++)
            cpu_owners[ai] = atoms[ai].owner_idx = i;

	num_nodes = nodes.size();
    num_layers = max_layer + 1;
    
	cudaMalloc(&device_nodes, sizeof(segment_node)*nodes.size());
	cudaMemcpy(device_nodes, &nodes[0], sizeof(segment_node)*nodes.size(), cudaMemcpyHostToDevice);

	cudaMalloc(&force_torques, sizeof(gfloat4p)*nodes.size());
	cudaMemset(force_torques, 0, sizeof(gfloat4p)*nodes.size());

    cudaMalloc(&owners, sizeof(uint) * natoms);
    cudaMemcpy(owners, cpu_owners, sizeof(uint) * natoms, cudaMemcpyHostToDevice);


    /* assert(num_nodes == ligand.count_torsions() + 1); */
    
    for(uint &o : cpu_owners){
        segment_node &n = nodes[o];
        assert(&o >= &cpu_owners[n.begin] &&
               &o < &cpu_owners[n.end]);
        assert(o < nodes.size());
        
    }
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
    _derivative((gfloat4 *) coords, (gfloat4 *) forces, c);
}

__device__
void tree_gpu::_derivative(const gfloat4 *coords,const gfloat4* forces, float *c){
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
    
    if(tid && tid < num_nodes){
        segment_node &n = device_nodes[tid];
        parent_id = n.parent;
        parent_origin = device_nodes[parent_id].origin;
        layer = n.layer;
        axis = n.axis;
        origin = n.origin;
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
    
    if(tid && tid < num_nodes)
        c[6 + tid - 1] = ft.second * axis;
    else if(tid == 0) {
        ft = force_torques[0];
	    c[0] = ft.first[0];
	    c[1] = ft.first[1];
	    c[2] = ft.first[2];

	    c[3] = ft.second[0];
	    c[4] = ft.second[1];
	    c[5] = ft.second[2];
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

__device__
void tree_gpu::set_conf(const vec *atom_coords, vec *coords,
                        const conf_info *c){
    _set_conf((marked_coord *) atom_coords, (gfloat4 *) coords, c);
}


__device__
void tree_gpu::_set_conf(const marked_coord *atom_coords, gfloat4 *coords,
                         const conf_info *c){
    
	uint tid = threadIdx.x;
    segment_node &n = device_nodes[tid];
    marked_coord a = atom_coords[tid];

    uint layer = 0;
    segment_node *parent = NULL;
    fl torsion = 0;

    
    if(tid == 0){
		for(unsigned i = 0; i < 3; i++)
			n.origin[i] = c->position[i];
		n.set_orientation(c->orientation[0], c->orientation[1],
                          c->orientation[2], c->orientation[3]);
    }else if(tid < num_nodes){
        layer = n.layer;
        parent = &device_nodes[n.parent];
        torsion = c->torsions[tid - 1];
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
