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
#include <boost/unordered_set.hpp>
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
  explicit MolGridDataLayer(const LayerParameter& param) :
      BaseDataLayer<Dtype>(param), num_rotations(0), current_rotation(0),
      example_size(0), balanced(false), paired(false), inmem(false), resolution(0.5),
      dimension(23.5), radiusmultiple(1.5), fixedradius(0), randtranslate(0),
      binary(false), randrotate(false), dim(0), numgridpoints(0),
      numReceptorTypes(0), numLigandTypes(0), gpu_alloc_size(0),
      gpu_gridatoms(NULL), gpu_gridwhich(NULL) {}
  virtual ~MolGridDataLayer();
  virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "MolGridData"; }
  virtual inline int ExactNumBottomBlobs() const { return 0; }
  virtual inline int ExactNumTopBlobs() const { return 2+
      this->layer_param_.molgrid_data_param().has_affinity()+
      this->layer_param_.molgrid_data_param().has_rmsd(); }

  virtual inline void resetRotation() { current_rotation = 0; }

  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  //the following really shouldn't be recalculated each evaluation (not including gradients)
  void getReceptorAtoms(int batch_idx, vector<float4>& atoms)
  {
    atoms.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
        atoms.push_back(mol.atoms[i]);
  }

  void getLigandAtoms(int batch_idx, vector<float4>& atoms)
  {
    atoms.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
        atoms.push_back(mol.atoms[i]);
  }

  void getReceptorChannels(int batch_idx, vector<short>& whichGrid)
  {
    whichGrid.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
        whichGrid.push_back(mol.whichGrid[i]);
  }

  void getLigandChannels(int batch_idx, vector<short>& whichGrid)
  {
    whichGrid.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
        whichGrid.push_back(mol.whichGrid[i]);
  }

  void getReceptorGradient(int batch_idx, vector<float3>& gradient)
  {
    gradient.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
        gradient.push_back(-mol.gradient[i]);
  }

  void getLigandGradient(int batch_idx, vector<float3>& gradient)
  {
    gradient.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
        gradient.push_back(-mol.gradient[i]);
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
        if (fixedradius <= 0)
          ainfo.w = xs_radius(t);
        else
          ainfo.w = fixedradius;
        float3 gradient(0,0,0);

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
        if (fixedradius <= 0)
          ainfo.w = xs_radius(t);
        else
          ainfo.w = fixedradius;
        float3 gradient(0,0,0);

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

  double getDimension() const { return dimension; }
  double getResolution() const { return resolution; }

  void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top, double scale) const;

 protected:

  typedef GridMaker::quaternion quaternion;
  typedef typename boost::multi_array_ref<Dtype, 4>  Grids;

  //for memory efficiency, only store a given string once and use the const char*
  class string_cache
  {
    boost::unordered_set<string> strings;
  public:
    const char* get(const string& s)
    {
      strings.insert(s);
      //we assume even as the set is resized that strings never get allocated
      return strings.find(s)->c_str();
    }
  };

  struct example
  {
    const char* receptor;
    const char* ligand;
    Dtype label;
    Dtype affinity;
    Dtype rmsd;

    example(): receptor(NULL), ligand(NULL), label(0), affinity(0), rmsd(0) {}
    example(Dtype l, const char* r, const char* lig): receptor(r), ligand(lig), label(l), affinity(0), rmsd(0) {}
    example(Dtype l, Dtype a, Dtype rms, const char* r, const char* lig): receptor(r), ligand(lig), label(l), affinity(a), rmsd(rms) {}
    example(string_cache& cache, string line, bool hasaffinity, bool hasrmsd);
  };

  //organize examples with respect to receptor
  struct paired_examples
  {
	  vector<string> receptors;
	  vector< vector<example> > actives; //indexed by receptor index first
	  vector< pair<unsigned, vector<example> > > decoys; //indexed by receptor index fist, includes current index into vector
	  vector< pair<unsigned, unsigned> > indices; //receptor/active indices; can be shuffled
	  unsigned curr_index; //where we are indices for getting examples

	  boost::unordered_map<string, unsigned> recmap; //map to receptor indices

	  paired_examples(): curr_index(0) {}

	  void add(const example& ex);

	  void shuffle_pairs(); //randomize - only necessary at start

	  //get next pair of examples, will shuffle as necessary
	  void next(example& active, example& decoy);

  };

  struct examples
  {
    bool store_all;
    bool store_actives_decoys;
    bool store_pairs;
    int count;

    string root_folder;

    vector<example> all;
    vector<example> actives;
    vector<example> decoys;
    paired_examples pairs;

    int all_index;
    int actives_index;
    int decoys_index;
    bool shuffle_on_wrap; //TODO this doesn't apply to pairs for now

    examples():
        store_all(true), store_actives_decoys(true), store_pairs(true), count(0),
        all_index(0), actives_index(0), decoys_index(0), shuffle_on_wrap(false) {}

    examples(bool all, bool actives_decoys, bool pairs, bool shuffle, string& root):
        store_all(all), store_actives_decoys(actives_decoys), store_pairs(pairs), count(0),
        all_index(0), actives_index(0), decoys_index(0), shuffle_on_wrap(shuffle),
        root_folder(root) {}

    void add(const example& ex);
    void shuffle_();
    void next(example& ex);
    void next_active(example& ex);
    void next_decoy(example& ex);
  };

  string_cache scache;
  examples data;
  examples data2;
  float data_ratio;
  bool two_data_sources;

  unsigned num_rotations;
  unsigned current_rotation;
  unsigned example_size; //channels*numgridpoints

  vector<int> top_shape;
  bool balanced;
  bool paired;
  bool inmem;

  //batch labels
  vector<Dtype> labels;
  vector<Dtype> affinities;
  vector<Dtype> rmsds;

  //grid stuff
  GridMaker gmaker;
  double resolution;
  double dimension;
  double radiusmultiple; //extra to consider past vdw radius
  double fixedradius;
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
  void set_grid_ex(Dtype *grid, const example& ex, const string& root_folder, mol_transform& transform, bool gpu);
  void set_grid_minfo(Dtype *grid, const mol_info& recatoms, const mol_info& ligatoms, mol_transform& transform, bool gpu);

  void forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu);
  void backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu);
  void Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  //stuff for outputing dx grids
  std::string getIndexName(const vector<int>& map, unsigned index) const;
  void outputDXGrid(std::ostream& out, Grids& grids, unsigned g, double scale) const;

};


}  // namespace caffe

#endif  // CAFFE_MOLGRID_DATA_LAYER_HPP_
