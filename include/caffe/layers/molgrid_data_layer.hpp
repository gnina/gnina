#ifndef CAFFE_MOLGRID_DATA_LAYER_HPP_
#define CAFFE_MOLGRID_DATA_LAYER_HPP_

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>

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

void test_set_atom_gradients();

namespace caffe {


//sample uniformly between 0 and 1
inline double unit_sample(rng_t *rng)
{
  return ((*rng)() - rng->min()) / double(rng->max() - rng->min());
}

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
      dimension(23.5), radiusmultiple(1.5), fixedradius(0), randtranslate(0), ligpeturb_translate(0),
      binary(false), randrotate(false), ligpeturb(false), dim(0), numgridpoints(0),
      numReceptorTypes(0), numLigandTypes(0), gpu_alloc_size(0),
      gpu_gridatoms(NULL), gpu_gridwhich(NULL), compute_atom_gradients(false) {}
  virtual ~MolGridDataLayer();
  virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "MolGridData"; }
  virtual inline int ExactNumBottomBlobs() const { return 0; }
  virtual inline int ExactNumTopBlobs() const { return 2+
      this->layer_param_.molgrid_data_param().has_affinity()+
      this->layer_param_.molgrid_data_param().has_rmsd()+
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

  void setLabels(Dtype pose, Dtype affinity=0, Dtype rmsd=0);
  void enableAtomGradients() { compute_atom_gradients = true; } //enable atom gradient computation

  void getReceptorAtoms(int batch_idx, vector<float4>& atoms);
  void getLigandAtoms(int batch_idx, vector<float4>& atoms);

  void getReceptorChannels(int batch_idx, vector<short>& whichGrid);
  void getLigandChannels(int batch_idx, vector<short>& whichGrid);
  void getReceptorGradient(int batch_idx, vector<float3>& gradient);
  void getReceptorTransformationGradient(int batch_idx, vec& force, vec& torque);
  void getMappedReceptorGradient(int batch_idx, unordered_map<string ,float3>& gradient);
  void getLigandGradient(int batch_idx, vector<float3>& gradient);
  void getMappedLigandGradient(int batch_idx, unordered_map<string, float3>& gradient);

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

  //set center to use for memory ligand
  template<typename Vec3>
  void setCenter(const Vec3& center) {
    mem_lig.center = center;
  }

  vec getCenter() const {
    return mem_lig.center;
  }

  //set in memory buffer
  template<typename Atom, typename Vec3>
  void setLigand(const vector<Atom>& ligand, const vector<Vec3>& coords, bool calcCenter=true)
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

    if(calcCenter || isnan(mem_lig.center[0])) {
      mem_lig.center = center;
    }
  }

  double getDimension() const { return dimension; }
  double getResolution() const { return resolution; }

  void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top, double scale) const;
  friend void ::test_set_atom_gradients();

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

  struct mol_transform;
  struct mol_info {
    vector<float4> atoms;
    vector<short> whichGrid; //separate for better memory layout on gpu
    vector<float3> gradient;
    vec center; //precalculate centroid, includes any random translation

    mol_info() { center[0] = center[1] = center[2] = NAN;}

    void append(const mol_info& a)
    {
      atoms.insert(atoms.end(), a.atoms.begin(), a.atoms.end());
      whichGrid.insert(whichGrid.end(), a.whichGrid.begin(), a.whichGrid.end());
      gradient.insert(gradient.end(), a.gradient.begin(), a.gradient.end());
    }

    void transform_and_append(const mol_info& a, const mol_transform& transform)
    {
      //copy atoms from a into this, transforming the coordinates according to transform
     // LOG(INFO) << "About to transform " << a.atoms.size() << " atoms";
      for(unsigned i = 0, n = a.atoms.size(); i < n; i++) {
        //non-coordinate stuff
        whichGrid.push_back(a.whichGrid[i]);
        gradient.push_back(a.gradient[i]); //NOT rotating, but that shouldn't matter, right?

        float4 atom = a.atoms[i];
        quaternion p(0, atom.x-a.center[0], atom.y-a.center[1], atom.z-a.center[2]);
        p = transform.Q * p * (conj(transform.Q) / norm(transform.Q));
        atom.x = p.R_component_2() + a.center[0] + transform.center[0];
        atom.y = p.R_component_3() + a.center[1] + transform.center[1];
        atom.z = p.R_component_4() + a.center[2] + transform.center[2];
        atoms.push_back(atom);

        //LOG(INFO) << "Transforming " << a.atoms[i].x<<","<<a.atoms[i].y<<","<<a.atoms[i].z<<" to "<<atom.x<<","<<atom.y<<","<<atom.z;
      }
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

  //6 numbers representing a transformation
  struct output_transform {
    Dtype x;
    Dtype y;
    Dtype z;
    Dtype pitch;
    Dtype yaw;
    Dtype roll;

    output_transform(): x(0), y(0), z(0), pitch(0), yaw(0), roll(0) {}

    output_transform(Dtype X, Dtype Y, Dtype, Dtype Z, const quaternion& Q): x(X), y(Y), z(Z) {
      set_from_quaternion(Q);
    }

    void set_from_quaternion(const quaternion& Q) {
      //convert to euler angles
      //https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Quaternion_to_Euler_Angles_Conversion
      // roll (x-axis rotation)
      double w = Q.R_component_1();
      double x = Q.R_component_2();
      double y = Q.R_component_3();
      double z = Q.R_component_4();

      double sinr = 2.0 * (w*x + y*z);
      double cosr = 1.0 - 2.0 * (x*x + y*y);
      roll = atan2(sinr, cosr);

      // pitch (y-axis rotation)
      double sinp = 2.0 * (w*y - z*x);
      if (fabs(sinp) >= 1)
        pitch = copysign(M_PI / 2, sinp);// use 90 degrees if out of range
      else
        pitch = asin(sinp);

      // yaw (z-axis rotation)
      double siny = 2.0 * (w*z + x*y);
      double cosy = 1.0 - 2.0 * (y*y + z*z);
      yaw = atan2(siny, cosy);
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

    //add upto randtranslate in displacement (plus or minus) along each direction
    void add_random_displacement(rng_t* rng, double randtranslate)
    {
      double offx = unit_sample(rng)*2.0-1.0;
      double offy = unit_sample(rng)*2.0-1.0;
      double offz = unit_sample(rng)*2.0-1.0;
      center[0] += offx * randtranslate;
      center[1] += offy * randtranslate;
      center[2] += offz * randtranslate;
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

      Q = quaternion(r1,r2,r3,r4);
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
  vector<output_transform> perturbations;

  //grid stuff
  GridMaker gmaker;
  double resolution;
  double dimension;
  double radiusmultiple; //extra to consider past vdw radius
  double fixedradius;
  double randtranslate;
  double ligpeturb_translate;
  bool ligpeturb_rotate;
  bool binary; //produce binary occupancies
  bool randrotate;
  bool ligpeturb; //for spatial transformer

  unsigned dim; //grid points on one side
  unsigned numgridpoints; //dim*dim*dim

  vector<int> rmap; //map atom types to position in grid vectors
  vector<int> lmap;
  unsigned numReceptorTypes;
  unsigned numLigandTypes;


  unsigned gpu_alloc_size;
  float4 *gpu_gridatoms;
  short *gpu_gridwhich;
  bool compute_atom_gradients;

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
  void set_grid_ex(Dtype *grid, const example& ex, const string& root_folder,
                    mol_transform& transform, output_transform& pertub, bool gpu);
  void set_grid_minfo(Dtype *grid, const mol_info& recatoms, const mol_info& ligatoms,
                    mol_transform& transform, output_transform& peturb, bool gpu);

  void setAtomGradientsGPU(GridMaker& gmaker, Dtype *diff, unsigned batch_size);
  void forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu);
  void backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu);
  void Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  //stuff for outputing dx grids
  std::string getIndexName(const vector<int>& map, unsigned index) const;
  void outputDXGrid(std::ostream& out, Grids& grids, unsigned g, double scale) const;

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
