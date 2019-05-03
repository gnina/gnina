#ifndef CAFFE_MOLGRID_DATA_LAYER_HPP_
#define CAFFE_MOLGRID_DATA_LAYER_HPP_

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <map>
#include <cmath>

#include <boost/array.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include "caffe/blob.hpp"
#include "caffe/data_transformer.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/base_data_layer.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/rng.hpp"

#include "gninasrc/lib/quaternion.h"
#include "gninasrc/lib/atom.h"
#include <libmolgrid/example_provider.h>
#include <libmolgrid/atom_typer.h>
#include <libmolgrid/grid_maker.h>

namespace caffe {


/*
 * @brief Provides data to the Net from files of examples.
 * libmolgrid is used to sample and select from provided example files
 * and to do the gridding.
 */

template<typename Dtype>
class MolGridDataLayer : public BaseDataLayer<Dtype> {
  public:
    typedef qt quaternion;

    explicit MolGridDataLayer(const LayerParameter& param) :
        BaseDataLayer<Dtype>(param) {
    }
    virtual ~MolGridDataLayer();
    virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual inline const char* type() const {
      return "MolGridData";
    }
    virtual inline int ExactNumBottomBlobs() const {
      return 0;
    }
    virtual inline int ExactNumTopBlobs() const {
      return 2 +
          ((this->layer_param_.molgrid_data_param().subgrid_dim() != 0) ||
              (this->layer_param_.molgrid_data_param().maxgroupsize() != 1)) +
          this->layer_param_.molgrid_data_param().has_affinity() +
          this->layer_param_.molgrid_data_param().has_rmsd() +
          this->layer_param_.molgrid_data_param().peturb_ligand();
    }

    virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
    virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
    virtual void setLabels(Dtype pose, Dtype affinity = 0, Dtype rmsd = 0);
    virtual void setLayerSpecificDims(int number_examples,
        vector<int>& label_shape, const vector<Blob<Dtype>*>& top);
    virtual void enableAtomGradients() {
      compute_atom_gradients = true;
    } //enable atom gradient computation

    virtual void clearLabels();
    void updateLabels(const std::vector<float>& l, bool hasaffinity, bool hasrmsd);
    virtual void updateLabels(Dtype pose, Dtype affinity = 0, Dtype rmsd = 0);

    virtual void copyToBlob(Dtype* src, size_t size, Blob<Dtype>* blob, bool gpu);
    virtual void copyToBlobs(const vector<Blob<Dtype>*>& top, bool hasaffinity, bool hasrmsd, bool gpu);

    void getReceptorAtoms(int batch_idx, vector<float4>& atoms);
    void getLigandAtoms(int batch_idx, vector<float4>& atoms);

    void getReceptorChannels(int batch_idx, vector<short>& whichGrid);
    void getLigandChannels(int batch_idx, vector<short>& whichGrid);
    void getReceptorGradient(int batch_idx, vector<gfloat3>& gradient);
    void getReceptorTransformationGradient(int batch_idx, vec& force,
        vec& torque);
    void getMappedReceptorGradient(int batch_idx,
        std::unordered_map<string, gfloat3>& gradient);
    void getLigandGradient(int batch_idx, vector<gfloat3>& gradient);
    void getMappedLigandGradient(int batch_idx,
        std::unordered_map<string, gfloat3>& gradient);

    virtual void setReceptor(const vector<atom>& receptor, const vec& translate =
        {}, const qt& rotate = {});
    virtual void setLigand(const vector<atom>& ligand, const vector<vec>& coords,
        bool calcCenter=true);

    //set center to use for memory ligand
    void setCenter(const vec& center) {
      CHECK_GT(batch_info.size(), 0) << "Empty batch info";
      batch_info[0].transform.set_rotation_center(gfloat3{center[0],center[1],center[2]});
      center_set = true;
    }

    vec getCenter() const {
      gfloat3 c= batch_info[0].transform.get_rotation_center();
      return vec(c.x,c.y,c.z);
    }

    vec getCenter(unsigned mol_idx) const {
      CHECK_GT(batch_info.size(), 0) << "Empty batch info";
      assert(mol_idx < batch_info.size());
      gfloat3 c = batch_info[mol_idx].transform.get_rotation_center();
      return vec(c.x,c.y,c.z);
    }

    double getDimension() const {
      return gmaker.get_dimension();
    }
    double getResolution() const {
      return gmaker.get_resolution();
    }

    gfloat3 getGridDims() const {
      return gmaker.get_grid_dims();
    }

    unsigned getNumChannels() const {
      return numchannels;
    }

    virtual void dumpDiffDX(const std::string& prefix, Blob<Dtype>* top,
        double scale) const;
    virtual void dumpGridDX(const std::string& prefix, Dtype* top,
        double scale = 1.0) const;

    virtual std::vector<std::string> getRecTypes() {
      return recTypes->get_type_names();
    }

    virtual std::vector<std::string> getLigTypes() {
      return ligTypes->get_type_names();
    }

    ///////////////////////////   PUBLIC DATA TYPES   //////////////////////////////

    struct mol_info {
      libmolgrid::CoordinateSet rec_atoms;
      libmolgrid::CoordinateSet lig_atoms;
      libmolgrid::Transform transform;
      libmolgrid::Transform ligand_perturbation;
      vector<gfloat3> rec_gradient; //todo: change to mgrid
      vector<gfloat3> lig_gradient;

      void setReceptor(const libmolgrid::CoordinateSet& c) {
        rec_atoms = c;
        rec_gradient.assign(c.size(), gfloat3(0,0,0));
      }

      void setLigand(const libmolgrid::CoordinateSet& c) {
        lig_atoms = c;
        lig_gradient.assign(c.size(), gfloat3(0,0,0));
      }

      //return max distance from centroid to any ligand atom
      double ligandRadius() const
      {
        //always return relative to centroid of this molecule, not any set center
        gfloat3 c3 = transform.get_rotation_center();
        vec c(c3.x,c3.y,c3.z);

        double maxdsq = 0.0;
        for (unsigned i = 0, n = lig_atoms.size(); i < n; i++) {
          vec pos(lig_atoms.coord[i][0],lig_atoms.coord[i][1],lig_atoms.coord[i][2]);
          pos -= c;
          double dsq = pos.norm_sqr();
          if (dsq > maxdsq)
            maxdsq = dsq;
        }
        return sqrt(maxdsq);
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

        output_transform() :
            x(0), y(0), z(0), a(0), b(0), c(0), d(0), roll(0), pitch(0), yaw(0) {
        }

        output_transform(Dtype X, Dtype Y, Dtype, Dtype Z, const qt& Q) :
            x(X), y(Y), z(Z) {
          set_from_quaternion(Q);
        }

        //modify in-place to values are class labels instead of actual values
        void discretize(double maxtrans, int bins);

        static unsigned size() {
          return sizeof(output_transform) / sizeof(Dtype);
        }

        void set_from_quaternion(const qt& Q);

      private:
        //discretize single value
        Dtype convert_to_label(Dtype value, double min, double max, int bins);
    };


    mol_info& getMolInfo(int batch_idx) {
      return batch_info[batch_idx];
    }

    const mol_info& getMolInfo(int batch_idx) const {
      return batch_info[batch_idx];
    }

    virtual void forward(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top, bool gpu);
    virtual void backward(const vector<Blob<Dtype>*>& top,
        const vector<Blob<Dtype>*>& bottom, bool gpu);
    virtual void Backward_relevance(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

    ///////////////////   PROTECTED DATA   ////////////////
  protected:

    //we are manually stratifying by file, this could be made more general-purpose and flexible
    //as an example_provider subclass, but this is all we need for now
    libmolgrid::ExampleProvider data;
    libmolgrid::ExampleProvider data2;

    float data_ratio = 0.0;
    unsigned example_size = 0; //channels*numgridpoints

    vector<int> top_shape;
    bool inmem = false;
    bool center_set = false;

    //batch labels split into individual vectors
    vector<Dtype> labels;
    vector<Dtype> affinities;
    vector<Dtype> rmsds;
    vector<output_transform> perturbations;

    //grid stuff
    libmolgrid::GridMaker gmaker;
    double randtranslate = 0;
    double ligpeturb_translate = 0;
    double jitter = 0;
    unsigned numposes = 1;
    bool ligpeturb_rotate = false;
    bool randrotate = false;
    bool ligpeturb = false; //for spatial transformer
    bool ignore_ligand = false; //for debugging

    unsigned dim = 0; //grid points on one side
    unsigned numgridpoints = 0; //dim*dim*dim
    unsigned numchannels = 0;

    unsigned numReceptorTypes = 0; //these are cached from rec/lig types
    unsigned numLigandTypes = 0;
    std::shared_ptr<libmolgrid::AtomTyper> recTypes;
    std::shared_ptr<libmolgrid::AtomTyper> ligTypes;

    bool compute_atom_gradients = false;

    //need to remember how mols were transformed for backward pass; store gradient as well
    vector<typename MolGridDataLayer<Dtype>::mol_info> batch_info;

    ////////////////////   PROTECTED METHODS   //////////////////////
    void set_grid_ex(Dtype *grid, const libmolgrid::Example& ex,
        typename MolGridDataLayer<Dtype>::mol_info& minfo,
        int pose, output_transform& pertub, bool gpu);
    virtual void set_grid_minfo(Dtype *grid,
        typename MolGridDataLayer<Dtype>::mol_info& minfo,
        output_transform& peturb, bool gpu);
    void setAtomGradientsGPU(libmolgrid::GridMaker& gmaker, Dtype *diff,
        unsigned batch_size);

    //stuff for outputing dx grids
    std::string getIndexName(const vector<int>& map, unsigned index) const;

};

}  // namespace caffe

//round coordinates to same precision as pdb
//for identifying atoms
template<typename T>
inline std::string xyz_to_string(T x, T y, T z)
    {
  //avoid negative zeros in string representation
  if (x == 0) x = 0;
  if (y == 0) y = 0;
  if (z == 0) z = 0;

  std::stringstream ss;
  ss << std::fixed << std::setprecision(3) << x;
  std::string rounded_x = ss.str();
  ss.str("");
  ss << std::fixed << std::setprecision(3) << y;
  std::string rounded_y = ss.str();
  ss.str("");
  ss << std::fixed << std::setprecision(3) << z;
  std::string rounded_z = ss.str();

  std::string xyz = rounded_x + rounded_y + rounded_z;
  return xyz;
}

#endif  // CAFFE_MOLGRID_DATA_LAYER_HPP_
