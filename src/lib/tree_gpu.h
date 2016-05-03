/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_TREE_H
#define VINA_TREE_H

#include "conf_gpu.h"
#include "atom.h"

struct frame_gpu {
	frame_gpu(const vec& origin_) : origin(origin_), orientation_q(qt_identity), orientation_m(quaternion_to_r3(qt_identity)) {}
	vec local_to_lab(const vec& local_coords) const {
		vec tmp;
		tmp = origin + orientation_m*local_coords; 
		return tmp;
	}
	vec local_to_lab_direction(const vec& local_direction) const {
		vec tmp;
		tmp = orientation_m * local_direction;
		return tmp;
	}
	const qt& orientation() const { return orientation_q; }
	const vec& get_origin() const { return origin; }
protected:
	vec origin;
	void set_orientation(const qt& q) { // does not normalize the orientation
		orientation_q = q;
		orientation_m = quaternion_to_r3(orientation_q);
	}
	qt  orientation_q;
	mat orientation_m;

	frame_gpu() {}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & origin;
		ar & orientation_m;
		ar & orientation_q;
	}
};

struct atom_range_gpu {
    sz begin;
    sz end;
	atom_range_gpu(sz begin_, sz end_) : begin(begin_), end(end_) {}
	template<typename F>
	void transform(const F& f) {
		sz diff = end - begin;
		begin = f(begin);
		end   = begin + diff;
	}

	atom_range_gpu(): begin(0), end(0) {}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & begin;
		ar & end;
	}
};

struct atom_frame_gpu : public frame_gpu, public atom_range_gpu {
	atom_frame_gpu(const vec& origin_, sz begin_, sz end_) : frame_gpu(origin_), atom_range_gpu(begin_, end_) {}
	void set_coords(const atomv& atoms, vecv& coords) const {
		VINA_RANGE(i, begin, end)
			coords[i] = local_to_lab(atoms[i].coords);
	}
	vecp sum_force_and_torque(const vecv& coords, const vecv& forces) const {
		vecp tmp;
		tmp.first.assign(0);
		tmp.second.assign(0);
		VINA_RANGE(i, begin, end) {
			tmp.first  += forces[i];
			tmp.second += cross_product(coords[i] - origin, forces[i]);
		}
		return tmp;
	}

	atom_frame_gpu() {}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<frame_gpu>(*this);
		ar & boost::serialization::base_object<atom_range_gpu>(*this);
	}
};

struct rigid_body_gpu : public atom_frame_gpu {
	rigid_body_gpu() {}
	rigid_body_gpu(const vec& origin_, sz begin_, sz end_) : atom_frame_gpu(origin_, begin_, end_) {}
	void set_conf(const atomv& atoms, vecv& coords, const rigid_conf_gpu& c) {
		origin = c.position;
		set_orientation(c.orientation);
		set_coords(atoms, coords);
	}
	void count_torsions(sz& s) const {} // do nothing
	void set_derivative(const vecp& force_torque, rigid_change_gpu& c) const {
		c.position     = force_torque.first;
		c.orientation  = force_torque.second;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<atom_frame_gpu>(*this);
	}
};

struct axis_frame_gpu : public atom_frame_gpu {
	axis_frame_gpu(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : atom_frame_gpu(origin_, begin_, end_) {
		vec diff; diff = origin - axis_root;
		fl nrm = diff.norm();
		VINA_CHECK(nrm >= epsilon_fl);
		axis = (1/nrm) * diff;
	}
	void set_derivative(const vecp& force_torque, fl& c) const {
		c = force_torque.second * axis;
	}
protected:
	vec axis;

	axis_frame_gpu() {}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & axis;
		ar & boost::serialization::base_object<atom_frame_gpu>(*this);;
	}
};

struct segment_gpu : public axis_frame_gpu {
	segment_gpu() {} //for serialization
	segment_gpu(const vec& origin_, sz begin_, sz end_, const vec& axis_root, const frame_gpu& parent) : axis_frame_gpu(origin_, begin_, end_, axis_root) {
		VINA_CHECK(eq(parent.orientation(), qt_identity)); // the only initial parent orientation this c'tor supports
		relative_axis = axis;
		relative_origin = origin - parent.get_origin();
	}
	void set_conf(const frame_gpu& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		const fl torsion = *c;
		++c;
		origin = parent.local_to_lab(relative_origin);
		axis = parent.local_to_lab_direction(relative_axis);
		qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();
		quaternion_normalize_approx(tmp); // normalization added in 1.1.2
		//quaternion_normalize(tmp); // normalization added in 1.1.2
		set_orientation(tmp);
		set_coords(atoms, coords);
	}
	void count_torsions(sz& s) const {
		++s;
	}
private:
	vec relative_axis;
	vec relative_origin;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & relative_axis;
		ar & relative_origin;
		ar & boost::serialization::base_object<axis_frame_gpu>(*this);;
	}
};

struct first_segment_gpu : public axis_frame_gpu {
	first_segment_gpu() {}
	first_segment_gpu(const segment_gpu& s) : axis_frame_gpu(s) {}
	first_segment_gpu(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : axis_frame_gpu(origin_, begin_, end_, axis_root) {}
	void set_conf(const atomv& atoms, vecv& coords, fl torsion) {
		set_orientation(angle_to_quaternion(axis, torsion));
		set_coords(atoms, coords);
	}
	void count_torsions(sz& s) const {
		++s;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<axis_frame_gpu>(*this);;
	}
};

template<typename T> // T == branch
void branches_set_conf_gpu(std::vector<T>& b, const frame_gpu& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
	VINA_FOR_IN(i, b)
		b[i].set_conf(parent, atoms, coords, c);
}

template<typename T> // T == branch
void branches_derivative_gpu(const std::vector<T>& b, const vec& origin, const vecv& coords, const vecv& forces, vecp& out, flv::iterator& d) { // adds to out
	VINA_FOR_IN(i, b) {
		vecp force_torque = b[i].derivative(coords, forces, d);
		out.first  += force_torque.first;
		vec r; r = b[i].node.get_origin() - origin;
		out.second += cross_product(r, force_torque.first) + force_torque.second;
	}
}

template<typename T> // T == segment
struct tree_gpu {
	T node;
	std::vector< tree_gpu<T> > children;

	tree_gpu() {} //for serialization
	tree_gpu(const T& node_) : node(node_) {}
	void set_conf(const frame_gpu& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		node.set_conf(parent, atoms, coords, c);
		branches_set_conf_gpu(children, node, atoms, coords, c);
	}
	vecp derivative(const vecv& coords, const vecv& forces, flv::iterator& p) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		fl& d = *p; // reference
		++p;
		branches_derivative_gpu(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		return force_torque;
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & node;
		ar & children;
	}
};

typedef tree_gpu<segment_gpu> branch_gpu;
typedef std::vector<branch_gpu> branches_gpu;

template<typename Node> // Node == first_segment || rigid_body
struct heterotree_gpu {
	Node node;
	branches_gpu children;

	heterotree_gpu() {} //for serialization
	heterotree_gpu(const Node& node_) : node(node_) {}
	void set_conf(const atomv& atoms, vecv& coords, const ligand_conf_gpu& c) {
		node.set_conf(atoms, coords, c.rigid);
		flv::const_iterator p = c.torsions.begin();
		branches_set_conf_gpu(children, node, atoms, coords, p);
		assert(p == c.torsions.end());
	}
	void set_conf(const atomv& atoms, vecv& coords, const residue_conf_gpu& c) {
		flv::const_iterator p = c.torsions.begin();
		node.set_conf(atoms, coords, *p);
		++p;
		branches_set_conf_gpu(children, node, atoms, coords, p);
		assert(p == c.torsions.end());
	}
	void derivative(const vecv& coords, const vecv& forces, ligand_change_gpu& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();
		branches_derivative_gpu(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, c.rigid);
		assert(p == c.torsions.end());
	}
	void derivative(const vecv& coords, const vecv& forces, residue_change_gpu& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();
		fl& d = *p; // reference
		++p;
		branches_derivative_gpu(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		assert(p == c.torsions.end());
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & node;
		ar & children;
	}
};

template<typename T> // T = main_branch, branch, flexible_body
void count_torsions(const T& t, sz& s) {
	t.node.count_torsions(s);
	VINA_FOR_IN(i, t.children)
		count_torsions(t.children[i], s);
}

typedef heterotree_gpu<rigid_body_gpu> flexible_body_gpu;
typedef heterotree_gpu<first_segment_gpu> main_branch_gpu;

template<typename T> // T == flexible_body || main_branch
struct vector_mutable_gpu : public std::vector<T> {
	template<typename C>
	void set_conf(const atomv& atoms, vecv& coords, const std::vector<C>& c) { // C == ligand_conf || residue_conf
		VINA_FOR_IN(i, (*this))
			(*this)[i].set_conf(atoms, coords, c[i]);
	}
	szv count_torsions() const {
		szv tmp(this->size(), 0);
		VINA_FOR_IN(i, (*this))
			::count_torsions((*this)[i], tmp[i]);
		return tmp;
	}
	template<typename C>
	void derivative(const vecv& coords, const vecv& forces, std::vector<C>& c) const { // C == ligand_change || residue_change
		VINA_FOR_IN(i, (*this))
			(*this)[i].derivative(coords, forces, c[i]);
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & boost::serialization::base_object<std::vector<T> >(*this);
	}
};

template<typename T, typename F> // tree or heterotree - like structure
void transform_ranges_gpu(T& t, const F& f) {
	t.node.transform(f);
	VINA_FOR_IN(i, t.children)
		transform_ranges_gpu(t.children[i], f);
}

//these implement a way to restore I what I consider erroneous behavior
//on the part of Vina concerning OH groups
void set_fixed_rotable_hydrogens_gpu(bool set);
bool get_fixed_rotable_hydrogens_gpu();
#endif
