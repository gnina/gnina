#ifndef CAFFE_MOLGRID_DATA_LAYER_HPP_
#define CAFFE_MOLGRID_DATA_LAYER_HPP_

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <map>

#include <boost/array.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include "caffe/blob.hpp"
#include "caffe/data_transformer.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/base_data_layer.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/rng.hpp"

#include "gninasrc/lib/atom.h"
#include "gninasrc/lib/gridmaker.h"

void test_set_atom_gradients();
void test_subcube_grids();
struct atom_params;
template <typename atomT, typename MGridT, typename GridMakerT> 
  void set_cnn_grids(MGridT* mgrid, GridMakerT& gmaker, 
      std::vector<atom_params>& mol_atoms, std::vector<atomT>& mol_types);

namespace caffe {


//sample uniformly between 0 and 1
inline double unit_sample(rng_t *rng)
{
  return ((*rng)() - rng->min()) / double(rng->max() - rng->min());
}

/*
 * @brief Provides data to the Net from n-dimension  files of raw floating point data.
 * MolGridDataLayer is an abstract base class that allows us to interface with
 * other code without requiring it to know anything about the underlying
 * templating on the GridMaker type. That templating facilitates generating 
 * different kinds of grids that can be used on either the host or the device. 
 *
 * TODO(dox): thorough documentation for Forward and proto params.
 */

template<typename Dtype>
class MolGridDataLayer : public BaseDataLayer<Dtype> {
  public:
    explicit MolGridDataLayer(const LayerParameter& param) : 
      BaseDataLayer<Dtype>(param) {}
    virtual ~MolGridDataLayer() {}
    virtual vec getCenter() const = 0;
    virtual vec getCenter(unsigned mol_idx) const = 0;
    virtual double getDimension() const = 0;
    virtual double getResolution() const = 0;
    virtual unsigned getExampleSize() const = 0;
    virtual void getReceptorAtoms(int batch_idx, vector<float4>& atoms) = 0;
    virtual void getLigandAtoms(int batch_idx, vector<float4>& atoms) = 0;
    virtual void getMappedReceptorGradient(int batch_idx, unordered_map<string,
        float3>& gradient) = 0;
    virtual void getMappedLigandGradient(int batch_idx, unordered_map<string, 
        float3>& gradient) = 0;
    virtual void setReceptor(const vector<atom>& receptor, const vec& translate = 
        {}, const qt& rotate = {}) = 0;
    virtual void setLigand(const vector<atom>& ligand, const vector<vec>& coords, 
        bool calcCenter=true) = 0;
    virtual void setCenter(const vec& center) = 0;
    virtual void setLabels(Dtype pose, Dtype affinity=0, Dtype rmsd=0) = 0;
    virtual void enableAtomGradients() = 0;
    virtual void getReceptorChannels(int batch_idx, vector<short>& whichGrid) = 0;
    virtual void getLigandChannels(int batch_idx, vector<short>& whichGrid) = 0;
    virtual void getReceptorGradient(int batch_idx, vector<float3>& gradient) = 0;
    virtual void getReceptorTransformationGradient(int batch_idx, vec& force, vec& torque) = 0;
    virtual void getLigandGradient(int batch_idx, vector<float3>& gradient) = 0;
    virtual void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top, double scale) const = 0;
    virtual void dumpGridDX(const std::string& prefix, Dtype* top, double scale) const = 0;
    virtual std::vector<std::string> getRecTypes() = 0;
    virtual std::vector<std::string> getLigTypes() = 0;

    typedef qt quaternion;
    struct mol_transform;
    virtual mol_transform getMolTransform(int batch_idx) const = 0;
    struct mol_info {
      vector<float4> atoms;
      vector<short> whichGrid; //separate for better memory layout on gpu
      vector<float3> gradient;
      vec center; //precalculate centroid

      mol_info() { center[0] = center[1] = center[2] = NAN;}

      //add contents of a to this, incrementing whichGrid by offset
      void append(const mol_info& a, unsigned offset=0)
      {
        atoms.insert(atoms.end(), a.atoms.begin(), a.atoms.end());
        whichGrid.reserve(whichGrid.size()+a.whichGrid.size());
        for(auto g : a.whichGrid) {
            whichGrid.push_back(g+offset);
        }
        gradient.insert(gradient.end(), a.gradient.begin(), a.gradient.end());
      }

      //apply transformation in-place, modifying the coordinates of the mol
      //the center of the molecule is used for the rotation origin
      void apply_transform(const mol_transform& transform)
      {
        gfloat3 rcenter(center[0],center[1],center[2]);
        gfloat3 translate(-transform.center[0],-transform.center[1], -transform.center[2]);
       // LOG(INFO) << "Center: " << rcenter[0] << "," << rcenter[1] << "," << rcenter[2] << "\n";
       // LOG(INFO) << "Translate: " << translate[0] << "," << translate[1] << "," << translate[2] << "\n";
        for(unsigned i = 0, n = atoms.size(); i < n; i++) {
          float4 atom = atoms[i];
          float3 pt = transform.Q.transform(atom.x, atom.y, atom.z, rcenter, translate);
          //LOG(INFO) << "Transforming " << atom.x << "," << atom.y << "," << atom.z << " to " << pt.x << "," << pt.y << "," << pt.z << "\n";

          atom.x = pt.x;
          atom.y = pt.y;
          atom.z = pt.z;
          atoms[i] = atom;
        }

        //update center
        center += vec(translate.x,translate.y,translate.z);
      }

      //return max distance from centroid to any atom
      double radius() const
      {
        //always return relative to centroid of this molecule, not any set center
        vec c(0,0,0);
        for(unsigned i = 0, n = atoms.size(); i < n; i++) {
          float4 a = atoms[i];
          c += vec(a.x,a.y,a.z);
        }
        c /= atoms.size();

        double maxdsq = 0.0;
        for(unsigned i = 0, n = atoms.size(); i < n; i++) {
          float4 a = atoms[i];
          vec pos(a.x,a.y,a.z);
          pos -= c;
          double dsq = pos.norm_sqr();
          if(dsq > maxdsq)
            maxdsq = dsq;
        }
        return sqrt(maxdsq);
      }
    };

    struct mol_transform {
      mol_info mol;
      qt Q;  // rotation
      vec center; // translation is negative of this

      mol_transform() {
        mol = mol_info();
        Q = qt(0,0,0,0);
        center[0] = center[1] = center[2] = 0;
      }

      //zero translate, no rotate
      void reset()
      {
        Q = qt(1,0,0,0);
        center[0] = center[1] = center[2] = 0;
      }

      //add upto randtranslate in displacement (plus or minus) along each direction
      vec add_random_displacement(rng_t* rng, double randtranslate)
      {
        double offx = unit_sample(rng)*2.0-1.0;
        double offy = unit_sample(rng)*2.0-1.0;
        double offz = unit_sample(rng)*2.0-1.0;
        vec translation;
        translation[0] = offx * randtranslate;
        translation[1] = offy * randtranslate;
        translation[2] = offz * randtranslate;
        center[0] += translation[0];
        center[1] += translation[1];
        center[2] += translation[2];
        return translation;
      }

      //set random quaternion
      void set_random_quaternion(rng_t* rng)
      {
        //http://planning.cs.uiuc.edu/node198.html
        //sample 3 numbers from 0-1
        double u1 = unit_sample(rng);
        double u2 = unit_sample(rng);
        double u3 = unit_sample(rng);
        double sq1 = sqrt(1-u1);
        double sqr = sqrt(u1);
        double r1 = sq1*sin(2*M_PI*u2);
        double r2 = sq1*cos(2*M_PI*u2);
        double r3 = sqr*sin(2*M_PI*u3);
        double r4 = sqr*cos(2*M_PI*u3);

        Q = qt(r1,r2,r3,r4);
      }

    };

};

/*
 * @brief Provides data to the Net from n-dimension files of raw floating point
 * data. BaseMolGridDataLayer is derived from MolGridDataLayer but can be
 * instantiated; all other MolGridData classes derive from it and specialize
 * based on the GridMaker type. 
 *
 */
template<typename Dtype, class GridMakerT>
class BaseMolGridDataLayer : public MolGridDataLayer<Dtype> {
  public:
    explicit BaseMolGridDataLayer(const LayerParameter& param) : 
      MolGridDataLayer<Dtype>(param), data(NULL), data2(NULL), data_ratio(0),
      num_rotations(0), current_rotation(0), 
      example_size(0), current_iter(0), inmem(false), resolution(0.5),
      dimension(23.5), radiusmultiple(1.5), fixedradius(0), randtranslate(0), ligpeturb_translate(0),
      jitter(0.0), numposes(1), ligpeturb_rotate(false),
      binary(false), randrotate(false), ligpeturb(false), ignore_ligand(false),
      ignore_rec(false), use_covalent_radius(false), dim(0), numgridpoints(0),
      numchannels(0), numReceptorTypes(0), numLigandTypes(0),
      gpu_alloc_size(0), gpu_gridatoms(NULL), gpu_gridwhich(NULL),
      compute_atom_gradients(false) {}
  virtual ~BaseMolGridDataLayer();
  virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void VSLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void InitParams(const MolGridDataParameter& param);
  virtual inline const char* type() const { return "MolGridData"; }
  virtual inline int ExactNumBottomBlobs() const { return 0; }
  virtual inline int ExactNumTopBlobs() const { return 2+
      ((this->layer_param_.molgrid_data_param().subgrid_dim() != 0) || 
       (this->layer_param_.molgrid_data_param().maxgroupsize() != 1)) +
      this->layer_param_.molgrid_data_param().has_affinity()+
      this->layer_param_.molgrid_data_param().has_rmsd()+
      (this->layer_param_.molgrid_data_param().affinity_reweight_stdcut() > 0) +
      this->layer_param_.molgrid_data_param().peturb_ligand();
  }
  virtual inline void resetRotation() { current_rotation = 0; }
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void setLabels(Dtype pose, Dtype affinity=0, Dtype rmsd=0);
  virtual void setLayerSpecificDims(int number_examples, 
    vector<int>& label_shape, const vector<Blob<Dtype>*>& top);
  virtual void enableAtomGradients() { compute_atom_gradients = true; } //enable atom gradient computation

  virtual void clearLabels() {
    labels.clear();
    affinities.clear();
    rmsds.clear();
    weights.clear();
  }

  virtual void updateLabels(Dtype pose, Dtype affinity=0, Dtype rmsd=0, Dtype weight=0) {
    for(unsigned p = 0; p < numposes; p++) {
      labels.push_back(pose);
      affinities.push_back(affinity);
      rmsds.push_back(rmsd);
      weights.push_back(weight);
    }
  }

  virtual void updateTranslations(vec&& translation) {}

  virtual void copyToBlob(Dtype* src, size_t size, Blob<Dtype>* blob, bool gpu) {
    Dtype* dst = nullptr;
    if (gpu)
      dst = blob->mutable_gpu_data();
    else
      dst = blob->mutable_cpu_data();
    caffe_copy(size, src, dst);
  }

  virtual void copyToBlobs(const vector<Blob<Dtype>*>& top, bool hasaffinity, bool hasrmsd, 
      bool hasweights, bool gpu) {
    unsigned weighti = top.size()-1-ligpeturb;
    unsigned rmsdi = 2+hasaffinity;
    copyToBlob(&labels[0], labels.size(), top[1], gpu);
    if(hasaffinity) 
      copyToBlob(&affinities[0], affinities.size(), top[2], gpu);
    if(hasrmsd)
      copyToBlob(&rmsds[0], rmsds.size(), top[rmsdi], gpu);
    if(hasweights)
      copyToBlob(&weights[0], weights.size(), top[weighti], gpu);
  
    if(ligpeturb)
      copyToBlob((Dtype*)&perturbations[0], perturbations.size()*perturbations[0].size(), top.back(), gpu);
  }

  void getReceptorAtoms(int batch_idx, vector<float4>& atoms);
  void getLigandAtoms(int batch_idx, vector<float4>& atoms);

  void getReceptorChannels(int batch_idx, vector<short>& whichGrid);
  void getLigandChannels(int batch_idx, vector<short>& whichGrid);
  void getReceptorGradient(int batch_idx, vector<float3>& gradient);
  void getReceptorTransformationGradient(int batch_idx, vec& force, vec& torque);
  void getMappedReceptorGradient(int batch_idx, unordered_map<string, float3>& gradient);
  void getLigandGradient(int batch_idx, vector<float3>& gradient);
  void getMappedLigandGradient(int batch_idx, unordered_map<string, float3>& gradient);

  typename MolGridDataLayer<Dtype>::mol_transform getMolTransform(int batch_idx) const { return batch_transform[batch_idx]; }

  //set in memory buffer
  //will apply translate and rotate iff rotate is valid
  void setReceptor(const vector<atom>& receptor, const vec& translate = {}, const qt& rotate = {})
  {
    mem_rec.atoms.clear();
    mem_rec.whichGrid.clear();
    mem_rec.gradient.clear();

    float3 c = make_float3(mem_lig.center[0], mem_lig.center[1], mem_lig.center[2]);
    float3 trans = make_float3(translate[0],translate[1],translate[2]);

    //receptor atoms
    for(unsigned i = 0, n = receptor.size(); i < n; i++)
    {
      const atom& a = receptor[i];
      smt t = a.sm;
      if (rmap[t] >= 0)
      {
        float4 ainfo;
        ainfo.x = a.coords[0];
        ainfo.y = a.coords[1];
        ainfo.z = a.coords[2];

        if(rotate.real() != 0) {
          //transform receptor coordinates
          //todo: move this to the gpu
          float3 pt = rotate.transform(ainfo.x, ainfo.y, ainfo.z, c, trans);
          ainfo.x = pt.x;
          ainfo.y = pt.y;
          ainfo.z = pt.z;
        }

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

  //set center to use for memory ligand
  void setCenter(const vec& center) {
    mem_lig.center = center;
  }

  //set in memory buffer
  void setLigand(const vector<atom>& ligand, const vector<vec>& coords, bool calcCenter=true)
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
        const vec& coord = coords[i];
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
      else if(t > 1) //don't warn about hydrogens
      {
        std::cerr << "Unsupported atom type " << smina_type_to_string(t);
      }
    }
    center /= acnt; //not ligand.size() because of hydrogens

    if(calcCenter || !isfinite(mem_lig.center[0])) {
      mem_lig.center = center;
    }
  }

  vec getCenter() const {
    return mem_lig.center;
  }

  vec getCenter(unsigned mol_idx) const {
    assert(mol_idx < batch_transform.size());
    return batch_transform[mol_idx].center;
  }

  double getDimension() const { return dimension; }
  double getResolution() const { return resolution; }
  unsigned getNumGridPoints() const { return numgridpoints; }
  unsigned getExampleSize() const { return example_size; }

  virtual void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top, double scale) const;
  virtual void dumpGridDX(const std::string& prefix, Dtype* top, double scale=1.0) const;
  friend void ::test_set_atom_gradients();
  friend void ::test_subcube_grids();
  template <typename atomT, typename MGridT, typename GridMakerU> 
    friend void ::set_cnn_grids(MGridT* mgrid, GridMakerU& gmaker, 
        std::vector<atom_params>& mol_atoms, std::vector<atomT>& mol_types);

 public:

  ///////////////////////////   PUBLIC DATA TYPES   //////////////////////////////
  typedef typename boost::multi_array_ref<Dtype, 4>  Grids;
  typedef typename boost::multi_array_ref<Dtype, 6>  RNNGrids;

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
    vector<const char*> ligands;
    Dtype label;
    Dtype affinity;
    Dtype rmsd;
    Dtype affinity_weight;
    int group;

    example(): receptor(NULL), label(0), affinity(0), rmsd(0), affinity_weight(1.0), group(-1) {}
    example(Dtype l, const char* r, const vector<const char*>& ligs): receptor(r), ligands(ligs), label(l), affinity(0), rmsd(0), group(-1)  {}
    example(Dtype l, Dtype a, Dtype rms, int gr, const char* r, const vector<const char*>& ligs, Dtype weight=1.0): receptor(r), ligands(ligs), label(l), affinity(a), rmsd(rms), affinity_weight(weight), group(gr)  {}
    example(string_cache& cache, string line, const MolGridDataParameter& param);
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
      CHECK_GT(examples.size(), 0) << "No valid stratified examples.";
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

      CHECK_GT(examples[currenti].size(), 0) << "No valid sub-stratified examples.";
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

  //group multiple grids into a single example, e.g. multiple frames of an MD
  //simulation that will be processed by an RNN. next() returns the next frame
  //for the next example in the batch; traversal is row-major with layout TxN
  template<class Provider> 
  class grouped_example_provider : public example_provider {

    typedef std::pair<typename std::vector<example>::iterator, typename std::vector<example>::iterator> location;

    Provider examples;
    unsigned batch_size;
    unsigned maxgroupsize;
    bool continuing;
    boost::unordered_map<int, vector<example>> frame_groups;
    std::deque<location> current_locations;

  public:
    grouped_example_provider(): examples(), batch_size(1), maxgroupsize(1), continuing(false) {}
    grouped_example_provider(const MolGridDataParameter& parm): examples(parm), 
                                                batch_size(parm.batch_size()), 
                                                maxgroupsize(parm.maxgroupsize()),
                                                continuing(false) {}
    //only add the first example for each group to examples; after that just
    //the filenames to the frame_groups map
    void add(const example& ex) {
      auto& group = ex.group;
      //TODO: when using c++17 switch to try_emplace
      auto ret = frame_groups.emplace(group, vector<example>());
      if (ret.second) 
        examples.add(ex);
      else {
        CHECK(frame_groups[group].size() <= (maxgroupsize - 1)) << "Frame group " << group << " size " << frame_groups[group].size()+1 << " exceeds max group size " << maxgroupsize << "."; //this could be handled, but it'd be messy the way things are now and it's proven to be a useful sanity check
        frame_groups[group].push_back(ex);
      }
    }

    void setup() {
      examples.setup();
      for (auto& group : frame_groups) {
        //if we have fewer than maxgroupsize examples for this group, pad with
        //ignore_labels
        for (unsigned idx=group.second.size(); idx<(maxgroupsize-1); ++idx) {
          group.second.push_back(example(-1, -1, -1, -1, NULL, {NULL}));
        }
      }
      CHECK_EQ(current_locations.size(), 0) << "All timesteps were not traversed prior to new provider setup";
      continuing = false;
    }

    void next(example& ex) {
      if (current_locations.size() < batch_size && !continuing) {
        examples.next(ex);
        current_locations.push_back(location(frame_groups[ex.group].begin(), 
              frame_groups[ex.group].end()));
      }
      else {
        auto progress = current_locations[0];
        current_locations.pop_front();
        ex = *progress.first;
        progress.first++;
        if (progress.first != progress.second) {
          current_locations.push_back(progress);
        }

        if (current_locations.size())
          continuing = true;
        else
          continuing = false;
      }
    }

    unsigned size() const
    {
      return examples.size();
    }
  };

  //N numbers representing a transformation
  struct output_transform {
    //contains only the N numbers representing the transformation, N is returned by size
    Dtype x;
    Dtype y;
    Dtype z;

    //three different representations of rotation, for experiments sake
    Dtype a;
    Dtype b;
    Dtype c;
    Dtype d;

    Dtype roll;
    Dtype pitch;
    Dtype yaw;

    output_transform(): x(0), y(0), z(0), a(0), b(0), c(0), d(0), roll(0), pitch(0), yaw(0) {}

    output_transform(Dtype X, Dtype Y, Dtype, Dtype Z, const qt& Q): x(X), y(Y), z(Z) {
      set_from_quaternion(Q);
    }


    //modify in-place to values are class labels instead of actual values
    void discretize(double maxtrans, int bins) {
      x = convert_to_label(x,-maxtrans,maxtrans,bins);
      y = convert_to_label(y,-maxtrans,maxtrans,bins);
      z = convert_to_label(z,-maxtrans,maxtrans,bins);

      a = convert_to_label(a,-1.0,1.0,bins);
      b = convert_to_label(b,-1.0,1.0,bins);
      c = convert_to_label(c,-1.0,1.0,bins);
      d = convert_to_label(d,-1.0,1.0,bins);

      roll = convert_to_label(roll,-M_PI,M_PI,bins);
      pitch = convert_to_label(pitch,-M_PI_2,M_PI_2,bins);
      yaw = convert_to_label(yaw,-M_PI,M_PI,bins);
    }
    static unsigned size() { return sizeof(output_transform)/sizeof(Dtype); }
    void set_from_quaternion(const qt& Q) {
      //convert to euler angles
      //https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Quaternion_to_Euler_Angles_Conversion
      // roll (x-axis rotation)
      a = Q.R_component_1();
      b = Q.R_component_2();
      c = Q.R_component_3();
      d = Q.R_component_4();

      double sinr = 2.0 * (a*b + c*d);
      double cosr = 1.0 - 2.0 * (b*b + c*c);
      roll = atan2(sinr, cosr);

      // pitch (y-axis rotation)
      double sinp = 2.0 * (a*c - d*b);

      pitch = 0.0;
      if (fabs(sinp) >= 1)
        pitch = copysign(M_PI / 2, sinp);// use 90 degrees if out of range
      else
        pitch = asin(sinp);

      // yaw (z-axis rotation)
      double siny = 2.0 * (a*d + b*c);
      double cosy = 1.0 - 2.0 * (c*c + d*d);
      yaw = atan2(siny, cosy);
    }

    private:
      //discretize single value
      Dtype convert_to_label(Dtype value, double min, double max, int bins)
      {
        int bin = bins*(value-min)/(max-min);
        if(bin < 0) bin = 0;
        if(bin >= bins) bin = bins-1;
        return bin;
      }
  };

  ///////////////////   PROTECTED DATA   ////////////////
 protected:
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
  unsigned current_iter;

  vector<int> top_shape;
  bool inmem;

  //batch labels
  vector<Dtype> labels;
  vector<Dtype> affinities;
  vector<Dtype> rmsds;
  vector<Dtype> weights;
  vector<output_transform> perturbations;

  //grid stuff
  GridMakerT gmaker;
  double resolution;
  double dimension;
  double radiusmultiple; //extra to consider past vdw radius
  double fixedradius;
  double randtranslate;
  double ligpeturb_translate;
  double jitter;
  unsigned numposes;
  bool ligpeturb_rotate;
  bool binary; //produce binary occupancies
  bool randrotate;
  bool ligpeturb; //for spatial transformer
  bool ignore_ligand; //for debugging
  bool ignore_rec;
  bool use_covalent_radius;

  unsigned dim; //grid points on one side
  unsigned numgridpoints; //dim*dim*dim
  unsigned numchannels;

  vector<int> rmap; //map atom types to position in grid vectors
  vector<int> lmap;
  
  unsigned numReceptorTypes;
  unsigned numLigandTypes;


  unsigned gpu_alloc_size;
  float4 *gpu_gridatoms;
  short *gpu_gridwhich;
  bool compute_atom_gradients;

  //need to remember how mols were transformed for backward pass
  vector<typename MolGridDataLayer<Dtype>::mol_transform> batch_transform;

  typedef boost::unordered_map<string, typename MolGridDataLayer<Dtype>::mol_info> MolCache;
  static MolCache recmolcache; //the cache is shared GLOBALLY
  static MolCache ligmolcache; //the cache is shared GLOBALLY

  typename MolGridDataLayer<Dtype>::mol_info mem_rec; //molecular data set programmatically with setReceptor
  typename MolGridDataLayer<Dtype>::mol_info mem_lig; //molecular data set programmatically with setLigand

  ////////////////////   PROTECTED METHODS   //////////////////////
  static void remove_missing_and_setup(vector<balanced_example_provider>& examples);
  void allocateGPUMem(unsigned sz);

  example_provider* create_example_data(const MolGridDataParameter& parm);
  void populate_data(const string& root_folder, const string& source, example_provider* data, bool hasaffinity, bool hasrmsd, bool hasgroup);

  typename MolGridDataLayer<Dtype>::quaternion axial_quaternion();

  bool add_to_minfo(const string& file, const vector<int>& atommap, unsigned mapoffset, smt t, float x, float y, float z,  typename MolGridDataLayer<Dtype>::mol_info& minfo);
  void load_cache(const string& file, const vector<int>& atommap, unsigned atomoffset, MolCache& molcache);
  void set_mol_info(const string& file, const vector<int>& atommap, unsigned atomoffset, typename MolGridDataLayer<Dtype>::mol_info& minfo);
  void set_grid_ex(Dtype *grid, const example& ex, const string& root_folder,
                    typename MolGridDataLayer<Dtype>::mol_transform& transform, 
                    int pose, output_transform& pertub, bool gpu);
  virtual void set_grid_minfo(Dtype *grid, 
      const typename MolGridDataLayer<Dtype>::mol_info& recatoms, 
      const typename MolGridDataLayer<Dtype>::mol_info& ligatoms,
                    typename MolGridDataLayer<Dtype>::mol_transform& transform, 
                    output_transform& peturb, bool gpu);
  void setAtomGradientsGPU(GridMakerT& gmaker, Dtype *diff, unsigned batch_size);

  virtual void forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu);
  virtual void backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu);
  virtual void Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  //stuff for outputing dx grids
  std::string getIndexName(const vector<int>& map, unsigned index) const;
  void outputDXGrid(std::ostream& out, Grids& grids, unsigned g, double scale, unsigned n) const;

  virtual std::vector<std::string> getRecTypes() {
    //map recmap idx to all the types that it corresponds to
    std::map<int, std::string> typemap;
    //final list of names, in recmap order with 1:1 correspondence to recmap
    std::vector<std::string> rectypes;
    //i is smina atom type
    for (int i=0; i<rmap.size(); ++i) {
      //at is recmap idx
      int at = rmap[i];
      if (at >= 0) {
        //it's in the map
        if (typemap.find(at) == typemap.end())
          typemap[at] = smina_type_to_string(static_cast<smt>(i));
        else {
          std::string other_types = typemap[at];
          other_types.push_back('_');
          typemap[at] = other_types + smina_type_to_string(static_cast<smt>(i));
        }
      }
    }
    for (auto it=typemap.begin(); it!=typemap.end(); ++it)
      rectypes.push_back("Rec_" + it->second);
    return rectypes;
  }

  virtual std::vector<std::string>getLigTypes() {
    //map ligmap idx to all the types that it corresponds to
    std::map<int, std::string> typemap;
    //final list of names, in ligmap order with 1:1 correspondence to ligmap
    std::vector<std::string> ligtypes;
    //i is smina atom type
    for (int i=0; i<lmap.size(); ++i) {
      //at is ligmap idx
      int at = lmap[i];
      if (at >= 0) {
        //it's in the map
        if (typemap.find(at) == typemap.end())
          typemap[at] = smina_type_to_string(static_cast<smt>(i));
        else {
          std::string other_types = typemap[at];
          other_types.push_back('_');
          typemap[at] = other_types + smina_type_to_string(static_cast<smt>(i));
        }
      }
    }
    for (auto it=typemap.begin(); it!=typemap.end(); ++it)
      ligtypes.push_back("Lig_" + it->second);
    return ligtypes;
  }
};

/*
 * @brief Provides data to the Net from n-dimension files of raw floating point
 * data. GenericMolGridData is used for the vanilla case of single inputs
 * processed as a single cubic grid and producing a single output (though there
 * is support for multiple poses, but not for processing them with an RNN -
 * they are treated as static alternate poses rather than dynamic snapshots)
 *
 */
template <typename Dtype>
class GenericMolGridDataLayer : public BaseMolGridDataLayer<Dtype, GridMaker> {
  public:
    explicit GenericMolGridDataLayer(const LayerParameter& param) : 
      BaseMolGridDataLayer<Dtype, GridMaker>(param) {}
    virtual ~GenericMolGridDataLayer() {};
    friend void ::test_set_atom_gradients();
    friend void ::test_subcube_grids();
    template <typename atomT, typename MGridT, typename GridMakerU> 
      friend void ::set_cnn_grids(MGridT* mgrid, GridMakerU& gmaker, 
          std::vector<atom_params>& mol_atoms, std::vector<atomT>& mol_types);

};

/*
 * @brief Provides data to the Net from n-dimension files of raw floating point
 * data. GroupedMolGridData is used to process multiple frames, each as a
 * separate grid (to be processed e.g. by an LSTM). 
 *
 */
template <typename Dtype>
class GroupedMolGridDataLayer : public BaseMolGridDataLayer<Dtype, GridMaker> {
  public:
    explicit GroupedMolGridDataLayer(const LayerParameter& param) : 
      BaseMolGridDataLayer<Dtype, GridMaker>(param), 
      maxgroupsize(param.molgrid_data_param().maxgroupsize()), 
      batch_size(param.molgrid_data_param().batch_size()), 
      example_idx(0), translations(batch_size) {
        CHECK_EQ(param.molgrid_data_param().subgrid_dim(), 0) << 
          "Subgrids and groups are mutually exclusive";
        unsigned input_chunksize = param.molgrid_data_param().maxchunksize();
        if (input_chunksize > 0) {
          CHECK_EQ(maxgroupsize % input_chunksize, 0) << 
            "maxchunksize must evenly divide maxgroupsize; pad if necessary";
          maxchunksize = input_chunksize;
        }
        else
          maxchunksize = maxgroupsize;
      }
    virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) { 
      NOT_IMPLEMENTED;
    }
    virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
      NOT_IMPLEMENTED;
    }

    virtual void Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
      NOT_IMPLEMENTED;
    }
    virtual void clearLabels() {
      seqcont.clear();
      BaseMolGridDataLayer<Dtype, GridMaker>::clearLabels();
    }

    virtual void updateLabels(Dtype pose, Dtype affinity=0, Dtype rmsd=0, Dtype weight=0) {
      for(unsigned p = 0; p < this->numposes; p++) {
        this->seqcont.push_back(example_idx < batch_size ? 0 : 1);
      }
      BaseMolGridDataLayer<Dtype, GridMaker>::updateLabels(pose, affinity, rmsd, weight);
    }

    virtual void updateTranslations(vec&& translation) {
      CHECK_LT(example_idx, batch_size) << "Only first frame should update translations";
      translations[example_idx] = translation;
    }

    protected:
    unsigned maxgroupsize;
    unsigned maxchunksize;
    unsigned batch_size;
    unsigned example_idx;
    vector<vec> translations; 
    vector<int> seqcont_shape;
    vector<Dtype> seqcont; //necessary for LSTM layer; indicates if a batch instance 
                           //is a continuation of a previous example sequence or 
                           //the beginning of a new one
    virtual void setLayerSpecificDims(int number_examples, 
      vector<int>& label_shape, const vector<Blob<Dtype>*>& top);

    virtual void copyToBlobs(const vector<Blob<Dtype>*>& top, bool hasaffinity, bool hasrmsd, 
        bool hasweights, bool gpu) {
      int idx = this->ExactNumTopBlobs() - this->ligpeturb - hasweights - 1;
      this->copyToBlob(&seqcont[0], seqcont.size(), top[idx], gpu);
      BaseMolGridDataLayer<Dtype, GridMaker>::copyToBlobs(top, hasaffinity, hasrmsd, hasweights, gpu);
    }

    virtual void set_grid_minfo(Dtype *data, 
        const typename MolGridDataLayer<Dtype>::mol_info& recatoms, 
        const typename MolGridDataLayer<Dtype>::mol_info& ligatoms,
        typename MolGridDataLayer<Dtype>::mol_transform& transform, 
        typename BaseMolGridDataLayer<Dtype, GridMaker>::output_transform& peturb, bool gpu);
};

/*
 * @brief Provides data to the Net from n-dimension files of raw floating point
 * data. SubcubeMolGridData is intended for use with the main branch LSTM layer
 * implementation. It decomposes the static inputs (such as those you would
 * process with GenericMolGridData) into disjoint subcubes and traverses them
 * as separate timesteps. 
 *
 */
template <typename Dtype>
class SubcubeMolGridDataLayer : public BaseMolGridDataLayer<Dtype, SubcubeGridMaker> {
  public:
    explicit SubcubeMolGridDataLayer(const LayerParameter& param) : 
      BaseMolGridDataLayer<Dtype, SubcubeGridMaker>(param) {
        CHECK_EQ(param.molgrid_data_param().maxgroupsize(), 1) << 
          "Subgrids and groups are mutually exclusive";
      }

    virtual ~SubcubeMolGridDataLayer() {};
      
    virtual void clearLabels() {
      seqcont.clear();
      BaseMolGridDataLayer<Dtype, SubcubeGridMaker>::clearLabels();
    }

    virtual void updateLabels(Dtype pose, Dtype affinity=0, Dtype rmsd=0, Dtype weight=0) {
      unsigned grids_per_dim = this->gmaker.grids_per_dim;
      unsigned ncubes = grids_per_dim * grids_per_dim * grids_per_dim;
      unsigned batch_size = this->gmaker.batch_size;
      unsigned batch_idx = this->gmaker.batch_idx;
      int nexamples = ncubes * batch_size;
      bool duplicate = this->layer_param_.molgrid_data_param().duplicate_poses();
      if(duplicate) nexamples = batch_size*this->numposes;
      if (seqcont.size() < nexamples)
        seqcont.resize(nexamples);
      if (this->labels.size() < nexamples) 
        this->labels.resize(nexamples);
      if (this->affinities.size() < nexamples)
        this->affinities.resize(nexamples);
      if (this->rmsds.size() < nexamples)
        this->rmsds.resize(nexamples);
      if (this->weights.size() < nexamples)
        this->weights.resize(nexamples);

      //this need to be TxN, so we end up having to write a column at a time, 
      //sadly
      for (size_t cube_id = 0; cube_id < ncubes; ++cube_id) {
        unsigned idx = cube_id * batch_size + batch_idx;
        for(unsigned p = 0; p < this->numposes; p++) {
          unsigned pose_idx = idx + p;
          if (cube_id == 0)
            seqcont[pose_idx] = 0;
          else
            seqcont[pose_idx] = 1;
          this->labels[pose_idx] = pose;
          this->affinities[pose_idx] = affinity;
          this->rmsds[pose_idx] = rmsd;
        }
      }
    }
    
    virtual void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top, double scale) const;

    friend void ::test_set_atom_gradients();
    friend void ::test_subcube_grids();
    template <typename atomT, typename MGridT, typename GridMakerU> 
      friend void ::set_cnn_grids(MGridT* mgrid, GridMakerU& gmaker, 
          std::vector<atom_params>& mol_atoms, std::vector<atomT>& mol_types);
  protected:
  vector<int> seqcont_shape;
  vector<Dtype> seqcont; //necessary for LSTM layer; indicates if a batch instance 
                         //is a continuation of a previous example sequence or 
                         //the beginning of a new one
  virtual void setLayerSpecificDims(int number_examples, 
    vector<int>& label_shape, const vector<Blob<Dtype>*>& top);

  virtual void copyToBlobs(const vector<Blob<Dtype>*>& top, bool hasaffinity, bool hasrmsd, 
      bool hasweights, bool gpu) {
    int idx = this->ExactNumTopBlobs() - this->ligpeturb - hasweights - 1;
    this->copyToBlob(&seqcont[0], seqcont.size(), top[idx], gpu);
    BaseMolGridDataLayer<Dtype, SubcubeGridMaker>::copyToBlobs(top, hasaffinity, hasrmsd, hasweights, gpu);
  }

  virtual void backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, 
      bool gpu) {
    NOT_IMPLEMENTED;
  }
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) { 
    // backward(top, bottom, false);
    NOT_IMPLEMENTED;
  }
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    // backward(top, bottom, true);
    NOT_IMPLEMENTED;
  }

  virtual void Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    NOT_IMPLEMENTED;
  }
};

}  // namespace caffe

//round coordinates to same precision as pdb
//for identifying atoms
template<typename T>
static string xyz_to_string(T x, T y, T z)
{
    //avoid negative zeros in string representation
    if(x == 0) x = 0;
    if(y == 0) y = 0;
    if(z == 0) z = 0;

    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << x;
    std::string rounded_x = ss.str();
    ss.str("");
    ss << std::fixed << std::setprecision(3) << y;
    std::string rounded_y = ss.str();
    ss.str("");
    ss << std::fixed << std::setprecision(3) << z;
    std::string rounded_z = ss.str();

    string xyz = rounded_x + rounded_y + rounded_z;
    return xyz;
}


#endif  // CAFFE_MOLGRID_DATA_LAYER_HPP_
