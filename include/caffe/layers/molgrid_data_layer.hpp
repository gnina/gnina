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

#include "gnina/src/lib/atom_constants.h"

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
        example_size(0),balanced(false),inmem(false),data_avail(0),
				resolution(0.5), dimension(23.5), radiusmultiple(1.5), randtranslate(0),
				binary(false), randrotate(false), dim(0), numgridpoints(0),
				numReceptorTypes(0),numLigandTypes(0) {}
  virtual ~MolGridDataLayer();
  virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "NDimData"; }
  virtual inline int ExactNumBottomBlobs() const { return 0; }
  virtual inline int ExactNumTopBlobs() const { return 2; }

  virtual inline vector<Dtype>& getMemoryData() { return memdata; }
  virtual void memoryIsSet(); //safe to read memory data now

  virtual inline void resetRotation() { current_rotation = 0; }

  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top);
//  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
//      const vector<Blob<Dtype>*>& top);

 protected:

  typedef typename boost::math::quaternion<Dtype> quaternion;
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

  vector<Dtype> memdata; //if inmemory is set, get data from here
  int data_avail; //can be more than 1 if rotating
  boost::mutex mem_mutex; //for guarding condition variable
  boost::condition_variable mem_cond;

  //grid stuff
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


	vector<Dtype> labels;

	struct atom_info {
	  //the thought is this will map to a float4 on the gpu
	  vec coord;
	  float radius;
	};

	struct mol_info {
	  vector<atom_info> atoms;
	  vector<short> whichGrid; //separate for better memory layout on gpu
	  vec center; //precalculate centroid, includes any random translation
	  boost::array< pair<float, float>, 3> dims;

	  mol_info() { center[0] = center[1] = center[2] = 0;}

	  void append(const mol_info& a)
	  {
	    atoms.insert(atoms.end(), a.atoms.begin(), a.atoms.end());
	    whichGrid.insert(whichGrid.end(), a.whichGrid.begin(), a.whichGrid.end());
	  }

	  void setCenter(Dtype dimension, const vec& c)
	  {
	    center = c;
	    Dtype half = dimension/2.0;
	    for(unsigned i = 0; i < 3; i++)
	    {
	      dims[i].first = c[i] - half;
	      dims[i].second = c[i] + half;
	    }
	  }
	};

	boost::unordered_map<string, mol_info> molcache;

	Dtype calcPoint(const vec& coords, double ar, const vec& pt);

	pair<unsigned, unsigned> getrange(const pair<float, float>& dim, double c, double r);
  static void zeroGrids(Grids& grids);
  //set the relevant grid points for atom
  void set_atom(const mol_info& mol, const atom_info& atom, int whichgrid, const quaternion& Q, Grids& grids);

  //set the relevant grid points for passed info
  void set_atoms(const mol_info& mol, const quaternion& Q, Grids& grids);


	quaternion axial_quaternion();
	void set_mol_info(const string& file, const vector<int>& atommap, unsigned atomoffset, mol_info& minfo);
  void set_grid(Dtype *grid, example ex, bool gpu);

  void forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu);
};


}  // namespace caffe

#endif  // CAFFE_MOLGRID_DATA_LAYER_HPP_
