#ifndef CAFFE_MOLGRID_DATA_LAYER_HPP_
#define CAFFE_MOLGRID_DATA_LAYER_HPP_

#include <string>
#include <utility>
#include <vector>

#include <boost/array.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include "caffe/blob.hpp"
#include "caffe/data_transformer.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/base_data_layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "gninasrc/lib/atom_constants.h"
#include "gninasrc/lib/gridmaker.h"

namespace caffe {

/**
 * @brief Provides data to the Net from n-dimension  files of raw floating point data.
 *
 * TODO(dox): thorough documentation for Forward and proto params.
 */
template <typename Dtype>
class MolGridDataLayer : public BaseDataLayer<Dtype> {
 public:
  explicit MolGridDataLayer(const LayerParameter& param)
      : BaseDataLayer<Dtype>(param), actives_pos_(0),
        decoys_pos_(0), all_pos_(0), num_rotations(0), current_rotation(0),
        example_size(0),balanced(false),inmem(false),
				resolution(0.5), dimension(23.5), radiusmultiple(1.5), randtranslate(0),
				binary(false), randrotate(false), dim(0), numgridpoints(0),
				numReceptorTypes(0),numLigandTypes(0), gpu_alloc_size(0),
				gpu_gridatoms(NULL), gpu_gridwhich(NULL) {}
  virtual ~MolGridDataLayer();
  virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "MolGridData"; }
  virtual inline int ExactNumBottomBlobs() const { return 0; }
  virtual inline int ExactNumTopBlobs() const { return 2; }

  virtual inline void resetRotation() { current_rotation = 0; }

  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  vector<float4> getReceptorAtoms(int batch_idx)
  {
    vector<float4> atoms;
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
        atoms.push_back(mol.atoms[i]);
    return atoms;
  }

  vector<float4> getLigandAtoms(int batch_idx)
  {
    vector<float4> atoms;
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
        atoms.push_back(mol.atoms[i]);
    return atoms;
  }

  vector<short> getReceptorChannels(int batch_idx)
  {
    vector<short> whichGrid;
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
        whichGrid.push_back(mol.whichGrid[i]);
    return whichGrid;
  }

  vector<short> getLigandChannels(int batch_idx)
  {
    vector<short> whichGrid;
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
        whichGrid.push_back(mol.whichGrid[i]);
    return whichGrid;
  }

  vector<float3> getReceptorGradient(int batch_idx)
  {
    vector<float3> gradient;
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
        gradient.push_back(mol.gradient[i]);
    return gradient;
  }

  vector<float3> getLigandGradient(int batch_idx)
  {
    vector<float3> gradient;
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
        gradient.push_back(mol.gradient[i]);
    return gradient;
  }

  //set in memory buffer
  template<typename Atom>
  void setReceptor(const vector<Atom>& receptor)
  {
    //make this a template mostly so I don't have to pull in gnina atom class
    mem_rec.atoms.clear();
    mem_rec.whichGrid.clear();
    mem_rec.gradient.clear();

    //receptor atoms
    for(unsigned i = 0, n = receptor.size(); i < n; i++)
    {
      const Atom& a = receptor[i];
      smt t = a.sm;
      if (rmap[t] >= 0)
      {
        float4 ainfo;
        ainfo.x = a.coords[0];
        ainfo.y = a.coords[1];
        ainfo.z = a.coords[2];
        ainfo.w = xs_radius(t);
        float3 gradient;
        gradient.x = 0.0;
        gradient.y = 0.0;
        gradient.z = 0.0;
        mem_rec.atoms.push_back(ainfo);
        mem_rec.whichGrid.push_back(rmap[t]);
        mem_rec.gradient.push_back(gradient);
      }
    }
  }

  //set in memory buffer
  template<typename Atom, typename Vec3>
  void setLigand(const vector<Atom>& ligand, const vector<Vec3>& coords)
  {
    mem_lig.atoms.clear();
    mem_lig.whichGrid.clear();
    mem_lig.gradient.clear();

    //ligand atoms, grid positions offset and coordinates are specified separately
    vec center(0,0,0);
    unsigned acnt = 0;
    for(unsigned i = 0, n = ligand.size(); i < n; i++)
    {
      smt t = ligand[i].sm;
      if(lmap[t] >= 0)
      {
        const Vec3& coord = coords[i];
        float4 ainfo;
        ainfo.x = coord[0];
        ainfo.y = coord[1];
        ainfo.z = coord[2];
        ainfo.w = xs_radius(t);
        float3 gradient;
        gradient.x = 0.0;
        gradient.y = 0.0;
        gradient.z = 0.0;
        mem_lig.atoms.push_back(ainfo);
        mem_lig.whichGrid.push_back(lmap[t]+numReceptorTypes);
        mem_lig.gradient.push_back(gradient);
        center += coord;
        acnt++;
      }
    }
    center /= acnt; //not ligand.size() because of hydrogens

    mem_lig.center = center;
  }


 protected:

  typedef GridMaker::quaternion quaternion;
  typedef typename boost::multi_array_ref<Dtype, 4>  Grids;

  struct example
	{
  	string receptor;
  	string ligand;
  	Dtype label;

  	example(): label(0) {}
  	example(Dtype l, const string& r, const string& lig): receptor(r), ligand(lig), label(l) {}
	};

  virtual void Shuffle();

  vector<example> actives_;
  vector<example> decoys_;
  vector<example> all_;
  string root_folder;
  int actives_pos_, decoys_pos_, all_pos_;
  unsigned num_rotations;
  unsigned current_rotation;
  unsigned example_size; //channels*numgridpoints
  vector<int> top_shape;
  bool balanced;
  bool inmem;
  vector<Dtype> labels;

  //grid stuff
  GridMaker gmaker;
  double resolution;
  double dimension;
  double radiusmultiple; //extra to consider past vdw radius
  double randtranslate;
  bool binary; //produce binary occupancies
  bool randrotate;

  unsigned dim; //grid points on one side
  unsigned numgridpoints; //dim*dim*dim

  vector<int> rmap; //map atom types to position in grid vectors
  vector<int> lmap;
  unsigned numReceptorTypes;
  unsigned numLigandTypes;


  unsigned gpu_alloc_size;
  float4 *gpu_gridatoms;
  short *gpu_gridwhich;

  void allocateGPUMem(unsigned sz);

  struct mol_info {
    vector<float4> atoms;
    vector<short> whichGrid; //separate for better memory layout on gpu
    vector<float3> gradient;
    vec center; //precalculate centroid, includes any random translation
    boost::array< pair<float, float>, 3> dims;

    mol_info() { center[0] = center[1] = center[2] = 0;}

    void append(const mol_info& a)
    {
      atoms.insert(atoms.end(), a.atoms.begin(), a.atoms.end());
      whichGrid.insert(whichGrid.end(), a.whichGrid.begin(), a.whichGrid.end());
      gradient.insert(gradient.end(), a.gradient.begin(), a.gradient.end());
    }
  };

  struct mol_transform {
    mol_info mol;
    quaternion Q;  // rotation
    vec center; // translation

    mol_transform() {
      mol = mol_info();
      Q = quaternion(0,0,0,0);
      center[0] = center[1] = center[2] = 0;
    }
  };

  //need to remember how mols were transformed for backward pass
  vector<mol_transform> batch_transform;

  boost::unordered_map<string, mol_info> molcache;
  mol_info mem_rec; //molecular data set programmatically with setReceptor
  mol_info mem_lig; //molecular data set programmatically with setLigand

  quaternion axial_quaternion();
  void set_mol_info(const string& file, const vector<int>& atommap, unsigned atomoffset, mol_info& minfo);
  void set_grid_ex(Dtype *grid, const example& ex, mol_transform& transform, bool gpu);
  void set_grid_minfo(Dtype *grid, const mol_info& recatoms, const mol_info& ligatoms, mol_transform& transform, bool gpu);

  void forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu);
  void backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu);
  //void Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
};


}  // namespace caffe

#endif  // CAFFE_MOLGRID_DATA_LAYER_HPP_
