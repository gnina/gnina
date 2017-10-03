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
#include "caffe/util/rng.hpp"

#include "gninasrc/lib/atom_constants.h"
#include "gninasrc/lib/gridmaker.h"

namespace caffe {

/*
 * @brief Provides data to the Net from n-dimension  files of raw floating point data.
 *
 * TODO(dox): thorough documentation for Forward and proto params.
 */
template <typename Dtype>
class MolGridDataLayer : public BaseDataLayer<Dtype> {
public:
  explicit MolGridDataLayer(const LayerParameter& param) :
      BaseDataLayer<Dtype>(param), data(NULL), data2(NULL), data_ratio(0),
      num_rotations(0), current_rotation(0),
      example_size(0), inmem(false), resolution(0.5),
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

  void getReceptorGradient(int batch_idx, vector<float3>& gradient, bool lrp = false)
  {
    gradient.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] < numReceptorTypes)
      {
        if(lrp)
        {
            gradient.push_back(mol.gradient[i]);
        }
        else
        {
            gradient.push_back(-mol.gradient[i]);
        }
      }
  }

  void getLigandGradient(int batch_idx, vector<float3>& gradient, bool lrp = false)
  {
    gradient.resize(0);
    mol_info& mol = batch_transform[batch_idx].mol;
    for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
      if (mol.whichGrid[i] >= numReceptorTypes)
      {
        if(lrp)
        {
            gradient.push_back(mol.gradient[i]);
        }
        else
        {
            gradient.push_back(-mol.gradient[i]);
        }
      }
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
      else
      {
        CHECK_LE(t, 1) << "Unsupported atom type " << smina_type_to_string(t);
      }
    }
    center /= acnt; //not ligand.size() because of hydrogens

    mem_lig.center = center;
  }

  double getDimension() const { return dimension; }
  double getResolution() const { return resolution; }

  void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top, double scale) const;

 protected:

  ///////////////////////////   PROTECTED DATA TYPES   //////////////////////////////
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

  //abstract class for storing training examples
  class example_provider
  {
  public:
    virtual void add(const example& ex) = 0;
    virtual void setup() = 0; //essentially shuffle if necessary
    virtual void next(example& ex) = 0;
    virtual unsigned size() const = 0;
    virtual ~example_provider() {}
  };

  //single array of examples, possibly shuffled
  class uniform_example_provider: public example_provider
  {
    vector<example> all;
    size_t current;
    bool randomize;

  public:
    uniform_example_provider(): current(0), randomize(false) {}
    uniform_example_provider(const MolGridDataParameter& parm): current(0)
    {
      randomize = parm.shuffle();
    }

    void add(const example& ex)
    {
      all.push_back(ex);
    }

    void setup()
    {
      current = 0;
      if(randomize) shuffle(all.begin(), all.end(), caffe::caffe_rng());
      CHECK_GT(all.size(), 0) << "Not enough examples (or at least the right kinds) in training set.";
    }

    void next(example& ex)
    {
      CHECK_LT(current, all.size()) << "Out of bounds error";
      ex = all[current];
      current++;
      if(current >= all.size())
      {
        setup(); //reset current and shuffle if necessary
      }
    }

    unsigned size() const { return all.size(); }
  };

  //sample uniformly from actives and decoys
  class balanced_example_provider: public example_provider
  {
    uniform_example_provider actives;
    uniform_example_provider decoys;
    size_t current;
    bool randomize;

  public:
    balanced_example_provider(): current(0), randomize(false) {}
    balanced_example_provider(const MolGridDataParameter& parm): actives(parm), decoys(parm), current(0)
    {
      randomize = parm.shuffle();
    }

    void add(const example& ex)
    {
      if (ex.label)
        actives.add(ex);
      else
        decoys.add(ex);
    }

    void setup()
    {
      current = 0;
      actives.setup();
      decoys.setup();
    }

    void next(example& ex)
    {
      //alternate between actives and decoys
      if(current % 2 == 0)
        actives.next(ex);
      else
        decoys.next(ex);

      current++;
    }

    unsigned size() const { return actives.size()+decoys.size(); }

    unsigned num_actives() const { return actives.size(); }
    unsigned num_decoys() const { return decoys.size(); }

    void next_active(example& ex)
    {
      actives.next(ex);
    }

    void next_decoy(example& ex)
    {
      decoys.next(ex);
    }
  };


  //partition examples by receptor and sample k times uniformly from each receptor
  //with k=2 and a balanced_provider you get paired examples from each receptor
  template<class Provider, int K=1>
  class receptor_stratified_example_provider: public example_provider
  {
    vector<Provider> examples;
    MolGridDataParameter p;
    boost::unordered_map<const char*, unsigned> recmap; //map to receptor indices

    size_t currenti, currentk; //position in array, and number of times sampling it
    bool randomize;

  public:
    receptor_stratified_example_provider(): currenti(0), currentk(0), randomize(false) {}
    receptor_stratified_example_provider(const MolGridDataParameter& parm): p(parm), currenti(0), currentk(0)
    {
      randomize = parm.shuffle();
    }

    void add(const example& ex)
    {
      if(recmap.count(ex.receptor) == 0)
      {
        //allocate
        recmap[ex.receptor] = examples.size();
        examples.push_back(Provider(p));
      }
      unsigned pos = recmap[ex.receptor];
      examples[pos].add(ex);
    }

    //NOTE: this has specializations for balanced/2 in the cpp file
    void setup()
    {
      CHECK_GT(K,0) << "Invalid sampling k for receptor_stratified_example_provider";
      currenti = 0; currentk = 0;
      for(unsigned i = 0, n = examples.size(); i < n; i++)
      {
        examples[i].setup();
      }
      //also shuffle receptors
      if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
    }

    void next(example& ex)
    {
      if(currentk >= K)
      {
        currentk = 0; //on to next receptor
        currenti++;
      }
      if(currenti >= examples.size())
      {
        currenti = 0;
        CHECK_EQ(currentk, 0) << "Invalid indices";
        if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
      }

      examples[currenti].next(ex);
      currentk++;
    }

    unsigned size() const
    {
      //no one said this had to be particularly efficient..
      unsigned ret = 0;
      for(unsigned i = 0, n = examples.size(); i < n; i++)
      {
        ret += examples[i].size();
      }
      return ret;
    }
  };

  //partition examples by affinity and sample uniformly from each affinity bin
  //affinities are binned by absolute value according to molgriddataparameters
  template<class Provider>
  class affinity_stratified_example_provider: public example_provider
  {
    vector<Provider> examples;
    size_t currenti;//position in array
    double min, max, step;

    //return bin for given affinity
    unsigned bin(double affinity) const
    {
      affinity = fabs(affinity);
      if(affinity < min) affinity = min;
      if(affinity >= max) affinity = max-FLT_EPSILON;
      affinity -= min;
      unsigned pos = affinity/step;
      return pos;
    }
  public:
    affinity_stratified_example_provider(): currenti(0), min(0), max(0), step(0) {}
    affinity_stratified_example_provider(const MolGridDataParameter& parm): currenti(0)
    {
      max = parm.stratify_affinity_max();
      min = parm.stratify_affinity_min();
      step = parm.stratify_affinity_step();
      CHECK_NE(min,max) << "Empty range for affinity stratification";
      unsigned maxbin = bin(max);
      CHECK_GT(maxbin, 0) << "Not enough bins";
      for(unsigned i = 0; i <= maxbin; i++)
      {
        examples.push_back(Provider(parm));
      }
    }

    void add(const example& ex)
    {
      unsigned i = bin(ex.affinity);
      CHECK_LT(i, examples.size()) << "Error with affinity stratification binning";
      examples[i].add(ex);
    }

    void setup()
    {
      currenti = 0;
      vector<Provider> tmp;
      for(unsigned i = 0, n = examples.size(); i < n; i++)
      {
        if(examples[i].size() > 0) {
          //eliminate empty buckets
          tmp.push_back(examples[i]);
          tmp.back().setup();
        }
	else {
	  LOG(INFO) << "Empty bucket " << i;
	}
      }
      swap(examples,tmp);
      CHECK_GT(examples.size(), 0) << "No examples in affinity stratification!";
    }

    void next(example& ex)
    {
      examples[currenti].next(ex);
      currenti = (currenti+1)%examples.size();
    }

    unsigned size() const
    {
      //no one said this had to be particularly efficient..
      unsigned ret = 0;
      for(unsigned i = 0, n = examples.size(); i < n; i++)
      {
        ret += examples[i].size();
      }
      return ret;
    }
  };

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

  ///////////////////   PROTECTED DATA   ////////////////

  string_cache scache;

  //we are manually stratifying by file, this could be made more general-purpose and flexible
  //as an example_provider subclass, but this is all we need for now
  example_provider *data;
  example_provider *data2;
  //store exampels without the root folder prefix to save memory
  //(which means they must be unique even without the prefix!)
  string root_folder;
  string root_folder2;
  float data_ratio;

  unsigned num_rotations;
  unsigned current_rotation;
  unsigned example_size; //channels*numgridpoints

  vector<int> top_shape;
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

  //need to remember how mols were transformed for backward pass
  vector<mol_transform> batch_transform;

  boost::unordered_map<string, mol_info> molcache;
  mol_info mem_rec; //molecular data set programmatically with setReceptor
  mol_info mem_lig; //molecular data set programmatically with setLigand


  ////////////////////   PROTECTED METHODS   //////////////////////
  static void remove_missing_and_setup(vector<balanced_example_provider>& examples);
  void allocateGPUMem(unsigned sz);

  example_provider* create_example_data(const MolGridDataParameter& parm);
  void populate_data(const string& root_folder, const string& source, example_provider* data, bool hasaffinity, bool hasrmsd);

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
