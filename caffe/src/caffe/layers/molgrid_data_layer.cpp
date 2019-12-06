#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>

#include "caffe/data_transformer.hpp"
#include "caffe/layers/base_data_layer.hpp"
#include "caffe/layers/molgrid_data_layer.hpp"
#include "caffe/util/benchmark.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/math_functions.hpp"

#include <boost/timer/timer.hpp>

#include <libmolgrid/grid_io.h>

using namespace libmolgrid;
using namespace std;

namespace caffe {


template<typename Dtype>
double MolGridDataLayer<Dtype>::mol_info::ligandRadius() const {
  //always return relative to centroid of this molecule, not any set center
  gfloat3 c3 = transform.get_rotation_center();
  vec c(c3.x,c3.y,c3.z);

  double maxdsq = 0.0;
  for (unsigned i = 0, n = orig_lig_atoms.size(); i < n; i++) {
    vec pos(orig_lig_atoms.coords(i,0),orig_lig_atoms.coords(i,1),orig_lig_atoms.coords(i,2));
    pos -= c;
    double dsq = pos.norm_sqr();
    if (dsq > maxdsq)
      maxdsq = dsq;
  }
  return sqrt(maxdsq);
}

//modify in-place to values are class labels instead of actual values
template<typename Dtype>
void MolGridDataLayer<Dtype>::output_transform::discretize(double maxtrans, int bins) {
  x = convert_to_label(x, -maxtrans, maxtrans, bins);
  y = convert_to_label(y, -maxtrans, maxtrans, bins);
  z = convert_to_label(z, -maxtrans, maxtrans, bins);

  a = convert_to_label(a, -1.0, 1.0, bins);
  b = convert_to_label(b, -1.0, 1.0, bins);
  c = convert_to_label(c, -1.0, 1.0, bins);
  d = convert_to_label(d, -1.0, 1.0, bins);

  roll = convert_to_label(roll, -M_PI, M_PI, bins);
  pitch = convert_to_label(pitch, -M_PI_2, M_PI_2, bins);
  yaw = convert_to_label(yaw, -M_PI, M_PI, bins);
}


template<typename Dtype>
void MolGridDataLayer<Dtype>::output_transform::set_from_quaternion(const qt& Q) {
  //convert to euler angles
  //https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Quaternion_to_Euler_Angles_Conversion
  // roll (x-axis rotation)
  a = Q.R_component_1();
  b = Q.R_component_2();
  c = Q.R_component_3();
  d = Q.R_component_4();

  double sinr = 2.0 * (a * b + c * d);
  double cosr = 1.0 - 2.0 * (b * b + c * c);
  roll = atan2(sinr, cosr);

  // pitch (y-axis rotation)
  double sinp = 2.0 * (a * c - d * b);

  pitch = 0.0;
  if (fabs(sinp) >= 1)
    pitch = copysign(M_PI / 2, sinp);  // use 90 degrees if out of range
  else
    pitch = asin(sinp);

  // yaw (z-axis rotation)
  double siny = 2.0 * (a * d + b * c);
  double cosy = 1.0 - 2.0 * (c * c + d * d);
  yaw = atan2(siny, cosy);
}

//discretize single value
template<typename Dtype>
Dtype MolGridDataLayer<Dtype>::output_transform::convert_to_label(Dtype value, double min, double max, int bins)
    {
  int bin = bins * (value - min) / (max - min);
  if (bin < 0) bin = 0;
  if (bin >= bins) bin = bins - 1;
  return bin;
}


template <typename Dtype>
MolGridDataLayer<Dtype>::~MolGridDataLayer<Dtype>() {
  //this->StopInternalThread();
}


//for in-memory inputs, set the desired label (for gradient computation)
//note that zero affinity/rmsd means to ignore these
template<typename Dtype>
void MolGridDataLayer<Dtype>::setLabels(Dtype pose, Dtype affinity, Dtype rmsd)
{
  clearLabels();
  labels.push_back(pose);
  affinities.push_back(affinity);
  rmsds.push_back(rmsd);
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getReceptorAtoms(int batch_idx, vector<gfloat3>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorAtoms";
  CoordinateSet& ratoms = batch_info[batch_idx].transformed_rec_atoms;
  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    gfloat3 a; //todo: change to float3
    a.x = ratoms.coords(i,0);
    a.y = ratoms.coords(i,1);
    a.z = ratoms.coords(i,2);
    atoms.push_back(a);
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getLigandAtoms(int batch_idx, vector<gfloat3>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getLigandAtoms";
  CoordinateSet& latoms = batch_info[batch_idx].transformed_lig_atoms;
  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    gfloat3 a; //todo: change to float3
    a.x = latoms.coords(i,0);
    a.y = latoms.coords(i,1);
    a.z = latoms.coords(i,2);
    atoms.push_back(a);
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getReceptorChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorChannels";
  CoordinateSet& ratoms = batch_info[batch_idx].transformed_rec_atoms;
  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    whichGrid.push_back(ratoms.type_index[i]);
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getLigandChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getLigandChannels";
  CoordinateSet& latoms = batch_info[batch_idx].transformed_lig_atoms;
  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    whichGrid.push_back(latoms.type_index[i]);
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getReceptorGradient(int batch_idx, vector<gfloat3>& gradient)
{
  gradient.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorGradient";
  CHECK(compute_receptor_gradients) << "Receptor gradients requested but not computed";
  CoordinateSet& ratoms = batch_info[batch_idx].transformed_rec_atoms;
  ManagedGrid<Dtype, 2>& grads = batch_info[batch_idx].rec_gradient;
  CHECK_EQ(ratoms.size(), grads.dimension(0));

  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    gradient.push_back(gfloat3(grads(i,0),grads(i,1),grads(i,2)));
  }
}

/*
 * Compute the transformation gradient of a rigid receptor around the center.
 * The first three numbers are the translation.  The next are the torque.
 */
template<typename Dtype>
void MolGridDataLayer<Dtype>::getReceptorTransformationGradient(int batch_idx, vec& force, vec& torque)
{
  force = vec(0,0,0);
  torque = vec(0,0,0);

  CHECK(compute_ligand_gradients) << "Ligand gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorTransformationGradient";
  CoordinateSet& ratoms = batch_info[batch_idx].transformed_rec_atoms;
  ManagedGrid<Dtype, 2>& grads = batch_info[batch_idx].rec_gradient;
  CHECK_EQ(ratoms.size(), grads.dimension(0));

  gfloat3 tc = batch_info[batch_idx].transform.get_rotation_center();
  vec c(tc.x,tc.y,tc.z);

  for (unsigned i = 0, n = ratoms.size(); i < n; ++i)
  {
    gfloat3 g = gfloat3(grads(i,0), grads(i,1), grads(i,2));
    gfloat3 a {ratoms.coords(i,0), ratoms.coords(i,1), ratoms.coords(i,2)};
    vec v(g.x,g.y,g.z);
    vec pos(a.x,a.y,a.z);

    force += v;
    torque += cross_product(pos - c, v);
  }
}


template<typename Dtype>
void MolGridDataLayer<Dtype>::getMappedReceptorGradient(int batch_idx, unordered_map<string, gfloat3>& gradient)
{
  CHECK(compute_receptor_gradients) << "Receptor gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getMappedReceptorGradient";

  CoordinateSet& ratoms = batch_info[batch_idx].transformed_rec_atoms;
  ManagedGrid<Dtype, 2>& grads = batch_info[batch_idx].rec_gradient;
  CHECK_EQ(ratoms.size(), grads.dimension(0));

  gradient.clear();
  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    string xyz = xyz_to_string(ratoms.coords(i,0), ratoms.coords(i,1), ratoms.coords(i,2));
    gradient[xyz] = gfloat3(grads(i,0), grads(i,1), grads(i,2));
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getMappedReceptorRelevance(int batch_idx, unordered_map<string, float>& relmap)
{
  CHECK(compute_receptor_gradients) << "Receptor gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getMappedReceptorRelevance";

  CoordinateSet& atoms = batch_info[batch_idx].transformed_rec_atoms;
  ManagedGrid<Dtype, 1>& relevance = batch_info[batch_idx].rec_relevance;
  CHECK_EQ(atoms.size(), relevance.size()) << "Receptor relevances incorrect size (not computed?)";
  relmap.clear();

  for (unsigned i = 0, n = atoms.size(); i < n; ++i) {
    string xyz = xyz_to_string(atoms.coords(i,0), atoms.coords(i,1), atoms.coords(i,2));
    relmap[xyz] = relevance[i];
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getLigandGradient(int batch_idx, vector<gfloat3>& gradient)
{
  CHECK(compute_ligand_gradients) << "Ligand gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getLigandGradient";

  gradient.resize(0);
  CoordinateSet& latoms = batch_info[batch_idx].transformed_lig_atoms;
  ManagedGrid<Dtype, 2>& grads = batch_info[batch_idx].lig_gradient;

  CHECK_EQ(latoms.size(), grads.dimension(0));

  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    gfloat3 g = gfloat3(grads(i, 0), grads(i, 1), grads(i, 2));
    gradient.push_back(g);
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getMappedLigandGradient(int batch_idx, unordered_map<string, gfloat3>& gradient)
{
  CHECK(compute_ligand_gradients) << "Ligand gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getMappedLigandGradient";

  CoordinateSet& latoms = batch_info[batch_idx].transformed_lig_atoms;
  ManagedGrid<Dtype, 2>& grads = batch_info[batch_idx].lig_gradient;
  CHECK_EQ(latoms.size(), grads.dimension(0));
  gradient.clear();

  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    gfloat3 g = gfloat3(grads(i,0), grads(i,1), grads(i,2));
    string xyz = xyz_to_string(latoms.coords(i,0), latoms.coords(i,1), latoms.coords(i,2));
    gradient[xyz] = g;
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getMappedLigandRelevance(int batch_idx, unordered_map<string, float>& relmap)
{
  CHECK(compute_ligand_gradients) << "Ligand gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getMappedLigandRelevance";

  CoordinateSet& latoms = batch_info[batch_idx].transformed_lig_atoms;
  ManagedGrid<Dtype, 1>& relevance = batch_info[batch_idx].lig_relevance;
  CHECK_EQ(latoms.size(), relevance.size()) << "Relevances incorrect size (not computed?)";
  relmap.clear();

  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    string xyz = xyz_to_string(latoms.coords(i,0), latoms.coords(i,1), latoms.coords(i,2));
    relmap[xyz] = relevance[i];
  }
}


//make sure can append to root_folder path
static string sanitize_path(const string& p)
{
  //make sure root folder(s) have trailing slash
  if (p.length() > 0 && p[p.length()-1] != '/')
    return p + "/";
  else
    return p;
}


//convert proto parameters to ExampleProviderSettings struct
//these are very similar because the same person wrote both
static ExampleProviderSettings settings_from_param(const MolGridDataParameter& param) {
  ExampleProviderSettings ret;
  ret.shuffle = param.shuffle();
  ret.balanced = param.balanced();
  ret.stratify_receptor = param.stratify_receptor();
  ret.labelpos = 0;
  ret.stratify_pos = 1;
  ret.stratify_abs = true;
  ret.stratify_min = param.stratify_affinity_min();
  ret.stratify_max = param.stratify_affinity_max();
  ret.stratify_step = param.stratify_affinity_step();
  ret.max_group_size = param.max_group_size();
  ret.group_batch_size = param.batch_size();
  ret.cache_structs = param.cache_structs();
  ret.add_hydrogens = param.addh();
  ret.duplicate_first = false;  // at least for now we do any duplication ourselves
  ret.data_root = param.root_folder(); //will have to overwrite if using root_folder2
  ret.recmolcache = param.recmolcache();
  ret.ligmolcache = param.ligmolcache();
  ret.num_copies = param.num_copies();
  return ret;
}

//wrap gninatyper to return fixed radius
class FixedRadiusTyper : public FileMappedGninaTyper {
    float radius = 0;
  public:
    FixedRadiusTyper(const FileMappedGninaTyper& b, double f_radius): FileMappedGninaTyper(b), radius(f_radius) {}

    virtual std::pair<int,float> get_atom_type_index(OpenBabel::OBAtom* a) const {
      auto res_rad = FileMappedGninaTyper::get_atom_type_index(a);
      res_rad.second = radius;
      return res_rad;
    }

    //map the type
    virtual std::pair<int,float> get_int_type(int t) const {
      auto res_rad = FileMappedGninaTyper::get_int_type(t);
      res_rad.second = radius;
      return res_rad;
    }
};

//create the appropriate atom typer given the specified strings
static std::shared_ptr<AtomTyper> create_atom_typer(const string& mapstr, const string& mapfile, const string& defaultmap, bool use_covalent_radius, double fixed_radius) {

  std::shared_ptr<FileMappedGninaTyper> ret;
  if(mapstr.size()) {
    //specify map programatically
    stringstream map(mapstr);
    ret = make_shared<FileMappedGninaTyper>(map, use_covalent_radius);
  }
  else if(mapfile.size()) {
    ret = make_shared<FileMappedGninaTyper>(mapfile, use_covalent_radius);
  }
  else {
    stringstream map(defaultmap.c_str());
    ret = make_shared<FileMappedGninaTyper>(map, use_covalent_radius);
  }

  if(fixed_radius > 0) {
    ret = make_shared<FixedRadiusTyper>(*ret, fixed_radius);
  }
  return ret;
}

const string default_recmap(R"(AliphaticCarbonXSHydrophobe 
AliphaticCarbonXSNonHydrophobe 
AromaticCarbonXSHydrophobe 
AromaticCarbonXSNonHydrophobe
Bromine Iodine Chlorine Fluorine
Nitrogen NitrogenXSAcceptor 
NitrogenXSDonor NitrogenXSDonorAcceptor
Oxygen OxygenXSAcceptor 
OxygenXSDonorAcceptor OxygenXSDonor
Sulfur SulfurAcceptor
Phosphorus 
Calcium
Zinc
GenericMetal Boron Manganese Magnesium Iron
)");

const string default_ligmap(R"(AliphaticCarbonXSHydrophobe 
AliphaticCarbonXSNonHydrophobe 
AromaticCarbonXSHydrophobe 
AromaticCarbonXSNonHydrophobe
Bromine Iodine
Chlorine
Fluorine
Nitrogen NitrogenXSAcceptor 
NitrogenXSDonor NitrogenXSDonorAcceptor
Oxygen OxygenXSAcceptor 
OxygenXSDonorAcceptor OxygenXSDonor
Sulfur SulfurAcceptor
Phosphorus
GenericMetal Boron Manganese Magnesium Zinc Calcium Iron
)");

//read in structure input and atom type maps
template <typename Dtype>
void MolGridDataLayer<Dtype>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {

  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();

  inmem = param.inmemory();
  double dimension = param.dimension();
  double resolution = param.resolution();
  bool binary = param.binary_occupancy();
  randtranslate = param.random_translate();
  randrotate = param.random_rotation();
  ligpeturb = param.peturb_ligand();
  ligpeturb_translate = param.peturb_ligand_translate();
  ligpeturb_rotate = param.peturb_ligand_rotate();
  jitter = param.jitter();
  ignore_ligand = param.ignore_ligand();
  double fixedradius = param.fixed_radius();
  bool use_covalent_radius = param.use_covalent_radius();
  bool hasaffinity = param.has_affinity();
  bool hasrmsd = param.has_rmsd();
  data_ratio = param.source_ratio();
  numposes = param.num_poses();
  group_size = param.max_group_size();
  chunk_size = param.max_group_chunk_size();
  if(chunk_size == 0) chunk_size = group_size;

  CHECK_EQ(group_size % chunk_size, 0) << "Chunk size must evenly divide group size.";
  CHECK_LE(fabs(remainder(dimension,resolution)), 0.001) << "Resolution does not evenly divide dimension.";

  gmaker.initialize(resolution, dimension, binary, param.radius_scaling(), param.gaussian_radius_multiple());

  if(param.radius_multiple() != 0) {
    LOG(WARNING) << "radius_multiple option is deprecated and ignored.  Use radius_scaling and gaussian_radius_multiple instead.";
  }

  dim = gmaker.get_grid_dims().x; //number of grid points on a side
  numgridpoints = dim*dim*dim;
  if(numgridpoints % 512 != 0)
    LOG(INFO) << "Total number of grid points (" << numgridpoints << ") is not evenly divisible by 512.";

  //shape must come from parameters
  int batch_size = param.batch_size();

  //initialize atom type maps
  string recmapfile = param.recmap();  //these are file names
  string ligmapfile = param.ligmap();

  //these are the actual contents
  string recmapstr = param.mem_recmap();
  string ligmapstr = param.mem_ligmap();

  recTypes = create_atom_typer(recmapstr, recmapfile, default_recmap, use_covalent_radius, fixedradius);
  ligTypes = create_atom_typer(ligmapstr, ligmapfile, default_ligmap, use_covalent_radius, fixedradius);

  numReceptorTypes = recTypes->num_types();
  numLigandTypes = ligTypes->num_types();
  ExampleProviderSettings settings = settings_from_param(param);
  data = ExampleProvider(settings, recTypes, ligTypes); //initialize types in data

  if(!inmem)
  {
    const string& source = param.source();
    const string& source2 = param.source2();
    CHECK_GT(source.length(), 0) << "No data source file provided";

    // Read source file(s) with labels and structures,
    // each line is label [group] [affinity] [rmsd] receptor_file ligand_file  (labels not optional)
    data.populate(source, 1+hasaffinity+hasrmsd);

    if(source2.length() > 0)
    {
      string root_folder = sanitize_path(param.root_folder());
      string root_folder2 = param.root_folder2();
      if(root_folder2.length() > 0) root_folder2 = sanitize_path(root_folder2);
      else root_folder2 = root_folder; //fall back on first

      CHECK_GE(data_ratio, 0) << "Must provide non-negative ratio for two data sources";
      settings.data_root = root_folder2;
      data2 = ExampleProvider(settings, recTypes, ligTypes);
      data2.populate(source2, 1+hasaffinity+hasrmsd);
    }

    LOG(INFO) << "Total examples: " << data.size() + data2.size();

    // Check if we would need to randomly skip a few data points
    if (param.rand_skip())
    {
      unsigned int skip = caffe_rng_rand() %  param.rand_skip();

      LOG(INFO) << "Skipping first " << skip << " data points from each source.";

      data.skip(skip);
      if(data2.size()) data2.skip(skip);
    }
  }
  else //in memory always batch size of 1
  {
    batch_size = 1;
  }


  CHECK_GT(batch_size, 0) << "Positive batch size required";
  //keep track of atoms and transformations for each example in batch
  //note that for grouped inputs, all the frames of the group share an info,
  //the molecular data of which gets overwritten
  batch_info.resize(batch_size);

  int number_examples = batch_size;
  bool duplicate = this->layer_param_.molgrid_data_param().duplicate_poses();

  if(duplicate) number_examples = batch_size*numposes;
  numchannels = numReceptorTypes+numLigandTypes;
  if(!duplicate && numposes > 1) numchannels = numReceptorTypes+numposes*numLigandTypes;

  //setup shape
  top_shape.clear();
  if(group_size > 1) top_shape.push_back(chunk_size); //output is T x B x C x W x H x D
  top_shape.push_back(number_examples);
  top_shape.push_back(numchannels);
  top_shape.push_back(dim);
  top_shape.push_back(dim);
  top_shape.push_back(dim);

  example_size = (numchannels)*numgridpoints;

  vector<int> label_shape;
  if(group_size > 1) label_shape.push_back(chunk_size);
  label_shape.push_back(number_examples); // [batch_size]

  top[0]->Reshape(top_shape);
  // Reshape label, affinity, rmsds, seqcont, libpeturb
  top[1]->Reshape(label_shape);

  unsigned idx = 2;
  if (hasaffinity) {
    top[idx]->Reshape(label_shape);
    idx++;
  }
  if (hasrmsd) {
      top[idx]->Reshape(label_shape);
      idx++;
  }

  if(group_size > 1) {
    vector<int> seqcont_shape{chunk_size, batch_size};
    top[idx]->Reshape(seqcont_shape);
    idx++;
  }

  if(ligpeturb) {
    vector<int> peturbshape(2);
    peturbshape[0] = batch_size;
    peturbshape[1] = output_transform::size(); //trans+orient
    top[idx]->Reshape(peturbshape);
    idx++;
  }
  CHECK_EQ(idx,top.size()) << "Inconsistent top size!";
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_ex(Dtype *data, const Example& ex,
    typename MolGridDataLayer<Dtype>::mol_info& minfo,
    int pose, output_transform& peturb, bool gpu, bool keeptransform)
{
  //set grid values for example
  //cache atom info
  //pose specifies which ligand pose to use (relevant if num_poses > 1)
  //if it is negative, use them all (but with distinct channels
  //data should be positioned at the start of the example
  bool doall = false;
  if(pose < 0) {
      doall = true;
      pose = 0;
  }

  CHECK_LT(pose, ex.sets.size()-1) << "Incorrect pose index";  //pose index doesn't include rec
  minfo.setReceptor(ex.sets[0], gpu);

  if (doall) { //first ligand pose already set, add others as discrete channels
    CHECK_EQ(pose,0)<< "Invalid pose specifier (internal error)";
    //merge non-receptor poses
    minfo.setLigand(ex.merge_coordinates(1), gpu);
  } else {
    //normal single ligand case
    minfo.setLigand(ex.sets[pose+1], gpu);
  }

  try {
    set_grid_minfo(data, minfo, peturb, gpu, keeptransform);
  } catch(...) {
    LOG(WARNING) << "Error processing";
    for(unsigned i = 0, n = ex.sets.size(); i < n; i++) {
      if(ex.sets[i].src) LOG(WARNING) << " " << ex.sets[i].src;
    }
    LOG(WARNING) << "\n";
    throw;
  }

}

//sample uniformly between 0 and 1
inline double unit_sample(rng_t *rng)
{
  return ((*rng)() - rng->min()) / double(rng->max() - rng->min());
}

//apply random noise (in place) to coordinates
static void apply_jitter(CoordinateSet& c, float jitter) {
  if(jitter <= 0) return;
  rng_t* rng = caffe_rng();
  bool ongpu = c.coords.ongpu();
  c.coords.tocpu();
  for (unsigned i = 0, n = c.size(); i < n; i++) {
    for(unsigned j = 0; j < 3; j++) {
      float diff = jitter * (unit_sample(rng)*2.0-1.0);
      c.coords(i,j) += diff;
    }
  }
  if(ongpu) c.coords.togpu();
}

//take a mol info, which includes receptor and ligand atoms
//and generate the appropriate grids into data
//applies jitter, peturbation, transformations as needed
//transform in minfo will be updated as appropriate
template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data,
    typename MolGridDataLayer<Dtype>::mol_info& minfo,
    output_transform& peturb, bool gpu, bool keeptransform)
{
  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
  bool fixcenter = param.fix_center_to_origin();

  //figure out transformation

  //what is the rotational center?
  gfloat3 rot_center{0,0,0};
  if (param.use_rec_center())
    rot_center = minfo.orig_rec_atoms.center();
  else
    rot_center = minfo.orig_lig_atoms.center();

  //limit translation by size of ligand/box
  float rtranslate = randtranslate;
  if (rtranslate)
  {
    float radius = minfo.ligandRadius();
    //don't let ligand atoms translate out of sphere inscribed in box
    if (ignore_ligand) radius = 0;
    double maxtrans = max(gmaker.get_dimension() / 2.0 - radius, 0.0);
    rtranslate = min(randtranslate, maxtrans);
  }

  if(!keeptransform) //for groups, only the first frame should set the transform
    minfo.transform = Transform(rot_center, rtranslate, randrotate);

  CoordinateSet& rec_atoms = minfo.transformed_rec_atoms;
  CoordinateSet& lig_atoms = minfo.transformed_lig_atoms;

  if(!minfo.transform.is_identity()) { //transformed are initialized to orig
    minfo.transform.forward(rec_atoms, rec_atoms);
    minfo.transform.forward(lig_atoms, lig_atoms);
  } 

  if (ligpeturb) { //apply ligand specific peturbation
    Transform P(lig_atoms.center(), ligpeturb_translate, ligpeturb_rotate);
    minfo.ligand_perturbation = P;
    P.forward(lig_atoms, lig_atoms); //peturb

    //store the inverse transformation
    gfloat3 trans = P.get_translation();
    peturb.x = -trans.x;
    peturb.y = -trans.y;
    peturb.z = -trans.z;
    Quaternion mgQ = P.get_quaternion();
    qt Q(mgQ.R_component_1(),mgQ.R_component_2(),mgQ.R_component_3(),mgQ.R_component_4());

    qt qinv = conj(Q) / norm(Q); //not Cayley, not euclidean norm - already squared
    peturb.set_from_quaternion(qinv);
  }

  //CoordinateSet atoms = ex.merge_coordinates(); //properly offsets types
  //apply jitter - this is on cpu and so inefficient, if it every turns out to be useful,
  //incorporate it into transform
  apply_jitter(rec_atoms, jitter);
  apply_jitter(lig_atoms, jitter);

  //set the grid center
  if(fixcenter) {
    minfo.grid_center = gfloat3(0,0,0);
  } else if(std::isfinite(grid_center.x)) {
    minfo.grid_center = grid_center;
  } else {
    minfo.grid_center = rot_center;
  }

  //compute grid from atoms
  //this would be slightly faster (10%) if we did rec and lig at the same time,
  //BUT allocating and copying into a combined buffer is significantly 
  //slower (50%) and for flexibility I want to keep them separate
  //if the buffer is preallocated and we mergeInto, it's only 10% slower, but still slower
  unsigned dim = gmaker.get_grid_dims().x;
  if (gpu)
  {
    Grid<Dtype, 4, true> recgrid(data, numReceptorTypes, dim, dim, dim);
    gmaker.forward(minfo.grid_center, rec_atoms, recgrid);
    if(!ignore_ligand) {
      Grid<Dtype, 4, true> liggrid(data+numgridpoints*numReceptorTypes, numchannels-numReceptorTypes, dim, dim, dim);
      gmaker.forward(minfo.grid_center, lig_atoms, liggrid);
    }
  }
  else
  {
    Grid<Dtype, 4, false> recgrid(data, numReceptorTypes, dim, dim, dim);
    gmaker.forward(minfo.grid_center, rec_atoms, recgrid);
    if(!ignore_ligand) {
      Grid<Dtype, 4, false> liggrid(data+numgridpoints*numReceptorTypes, numchannels-numReceptorTypes, dim, dim, dim);
      gmaker.forward(minfo.grid_center, lig_atoms, liggrid);
    }
  }

}


//dump dx files for every atom type, with files names starting with prefix
//only does first example
template<typename Dtype>
void MolGridDataLayer<Dtype>::dumpDiffDX(const std::string& prefix,  Blob<Dtype>* top, double scale) const {
  Dtype* diff = top->mutable_cpu_diff();
  Grid<Dtype, 4, false> grid(diff, numchannels, dim, dim, dim); //ignore rest of batch
  write_dx_grids(prefix, data.get_type_names(), grid, batch_info[0].grid_center, gmaker.get_resolution(), scale);
}


template <typename Dtype>
inline bool grid_empty(Dtype* data, size_t n_elements) {
  for (unsigned i=0; i<n_elements; ++i) {
    if (data[i] != 0.0) return false;
  }
  return true;
}

template bool grid_empty(float* data, size_t n_elements);
template bool grid_empty(double* data, size_t n_elements);

template <typename Dtype>
void MolGridDataLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top)
{
  forward(bottom, top, false);
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::dumpGridDX(const std::string& prefix, Dtype* d, double scale) const
{
  Grid<Dtype, 4, false> grid(d, numchannels, dim, dim, dim);
  write_dx_grids(prefix, data.get_type_names(), grid, batch_info[0].grid_center, gmaker.get_resolution(), scale);
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu)
{
  bool hasaffinity = this->layer_param_.molgrid_data_param().has_affinity();
  bool hasrmsd = this->layer_param_.molgrid_data_param().has_rmsd();
  bool duplicate = this->layer_param_.molgrid_data_param().duplicate_poses();
  int peturb_bins = this->layer_param_.molgrid_data_param().peturb_bins();
  double peturb_translate = this->layer_param_.molgrid_data_param().peturb_ligand_translate();

  Dtype *top_data = NULL;
  if(gpu)
    top_data = top[0]->mutable_gpu_data();
  else
    top_data = top[0]->mutable_cpu_data();

  perturbations.clear();

  unsigned batch_size;
  if (group_size>1) {
    batch_size = top_shape[1];
  }
  else
    batch_size = top_shape[0];

  unsigned div = 1;
  if(numposes > 1 && duplicate) div = numposes;
  batch_size /= div;
  if(duplicate) CHECK_EQ(top_shape[0] % numposes,0) << "Batch size not multiple of numposes??";
  output_transform peturb;

  //if in memory must be set programmatically
  //TODO: support for groups in cnn_scorer and elsewhere as needed
  if(inmem)
  {
    CHECK_GT(batch_info.size(), 0) << "Empty batch info";
    CHECK_EQ(group_size, 1) << "Groups not currently supported with structure in memory";
    if(batch_info[0].orig_rec_atoms.size() == 0) LOG(WARNING) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(batch_info[0].orig_lig_atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(top_data, batch_info[0], peturb, gpu, false);
    perturbations.push_back(peturb);

    CHECK_GT(labels.size(),0) << "Did not set labels in memory based molgrid";
  }
  else
  {
    clearLabels();

    //percent of batch from first data source
    unsigned dataswitch = batch_size;
    if (data2.size())
      dataswitch = batch_size*data_ratio/(data_ratio+1);

    for (int idx = 0, n = chunk_size*batch_size; idx < n; ++idx)
    {
      int batch_idx = idx % batch_size;
      Example ex;
      if (batch_idx < dataswitch) {
        data.next(ex);
      } else {
        data2.next(ex);
      }

      int step = idx / batch_size;
      int offset = ((batch_size * step) + batch_idx) * example_size;

      if(!duplicate) {
        updateLabels(ex.labels, hasaffinity, hasrmsd, ex.seqcont);
        set_grid_ex(top_data+offset, ex, batch_info[batch_idx], numposes > 1 ? -1 : 0, peturb, gpu, ex.seqcont);
        perturbations.push_back(peturb);
      }
      else {
        for(unsigned p = 0; p < numposes; p++) {
          updateLabels(ex.labels, hasaffinity, hasrmsd, ex.seqcont);
          int p_offset = batch_idx*(example_size*numposes)+example_size*p;
          set_grid_ex(top_data+p_offset, ex, batch_info[batch_idx], p, peturb, gpu, ex.seqcont);
          perturbations.push_back(peturb);
        }
      }
    }

  }

  copyToBlobs(top, hasaffinity, hasrmsd, gpu);

  if(peturb_bins > 0) {
    //discretize
    for(unsigned i = 0, n = perturbations.size(); i < n; i++) {
      perturbations[i].discretize(peturb_translate, peturb_bins);
    }

  }
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  backward(top, bottom, false);
}


/* backpropagates gradients onto atoms (note there is not actual bottom)
 * Only performed when compute_atom_gradients is true.
 */
template<typename Dtype>
void MolGridDataLayer<Dtype>::backward(const vector<Blob<Dtype>*>& top,
    const vector<Blob<Dtype>*>& bottom, bool gpu)
    {

  //propagate gradient grid onto atom positions
  if (compute_ligand_gradients || compute_receptor_gradients) {
    CHECK(numposes == 1) << "Atomic gradient calculation not supported with numposes != 1";
    unsigned batch_size = top_shape[0];
    CHECK(batch_info.size() == batch_size) << "Inconsistent batch sizes in backward";
    for (unsigned i = 0; i < batch_size; i++) {
      if (gpu) {
        Dtype *ptr = top[0]->mutable_gpu_diff()+example_size*i;
        mol_info& minfo = batch_info[i];
        //make sure grid is properly offset into types - do ligand and receptor separately
        //for more efficient computation if receptor isn't needed
        if(compute_ligand_gradients) {
          Grid<Dtype, 4, true> diff(ptr+numgridpoints*numReceptorTypes, numLigandTypes, dim, dim, dim); //offset by receptor types
          gmaker.backward(minfo.grid_center, minfo.transformed_lig_atoms, diff, minfo.lig_gradient.gpu());
          minfo.transform.backward(minfo.lig_gradient.gpu(), minfo.lig_gradient.gpu(), false);
        }
        if(compute_receptor_gradients) {
          Grid<Dtype, 4, true> diff(ptr, numReceptorTypes, dim, dim, dim);
          gmaker.backward(minfo.grid_center, minfo.transformed_rec_atoms, diff, minfo.rec_gradient.gpu());
          minfo.transform.backward(minfo.rec_gradient.gpu(), minfo.rec_gradient.gpu(), false);
        }
      } else {
        Dtype *ptr = top[0]->mutable_cpu_diff()+example_size*i;
        mol_info& minfo = batch_info[i];
        if(compute_ligand_gradients) {
          Grid<Dtype, 4, false> diff(ptr+numgridpoints*numReceptorTypes, numLigandTypes, dim, dim, dim);
          gmaker.backward(minfo.grid_center, minfo.transformed_lig_atoms, diff, minfo.lig_gradient.cpu());
          minfo.transform.backward(minfo.lig_gradient.cpu(), minfo.lig_gradient.cpu(), false);
        }
        if(compute_receptor_gradients) {
          Grid<Dtype, 4, false> diff(ptr, numReceptorTypes, dim, dim, dim);
          gmaker.backward(minfo.grid_center, minfo.transformed_rec_atoms, diff, minfo.rec_gradient.cpu());
          minfo.transform.backward(minfo.rec_gradient.cpu(), minfo.rec_gradient.cpu(), false);
        }
      }
    }
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  forward(bottom, top, true);
}

//since these are included when class is instantiated, manual instantiate
template void MolGridDataLayer<float>::Forward_gpu(const vector<Blob<float>*>& bottom, const vector<Blob<float>*>& top);
template void MolGridDataLayer<double>::Forward_gpu(const vector<Blob<double>*>& bottom, const vector<Blob<double>*>& top);

template<typename Dtype>
void MolGridDataLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
  backward(top, bottom, true);
}

template void MolGridDataLayer<float>::Backward_gpu(const vector<Blob<float>*>& top, const vector<bool>& propagate_down, const vector<Blob<float>*>& bottom);
template void MolGridDataLayer<double>::Backward_gpu(const vector<Blob<double>*>& top, const vector<bool>& propagate_down, const vector<Blob<double>*>& bottom);


template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  CHECK(numposes == 1) << "Relevance calculations not supported with numposes != 1";

  //propagate gradient grid onto atom positions
  unsigned batch_size = top_shape[0];
  CHECK_EQ(batch_size, batch_info.size()) << "Inconsistent batch sizes";

  Grid<Dtype, 5, false> diff(top[0]->mutable_cpu_diff(), batch_size, numchannels, dim, dim, dim);
  Grid<Dtype, 5, false> density(top[0]->mutable_cpu_data(), batch_size, numchannels, dim, dim, dim);

  for (int item_id = 0; item_id < batch_size; ++item_id) {
    mol_info& minfo = batch_info[item_id];
    minfo.rec_relevance = minfo.rec_relevance.resized(minfo.orig_rec_atoms.size());
    minfo.rec_relevance.fill_zero();
    minfo.lig_relevance = minfo.lig_relevance.resized(minfo.orig_lig_atoms.size());
    minfo.lig_relevance.fill_zero();

    gmaker.backward_relevance(minfo.grid_center, minfo.transformed_lig_atoms, density[item_id], diff[item_id], minfo.lig_relevance.cpu());
    gmaker.backward_relevance(minfo.grid_center, minfo.transformed_rec_atoms, density[item_id], diff[item_id], minfo.rec_relevance.cpu());
  }

}

template <typename Dtype>
void MolGridDataLayer<Dtype>::clearLabels() {
  labels.clear();
  affinities.clear();
  rmsds.clear();
  seqcont.clear();
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::updateLabels(const std::vector<float>& l, bool hasaffinity,
    bool hasrmsd, bool seq_continued) {
  float pose = 0, affinity = 0, rmsd = 0;
  unsigned n = l.size();

  //caffe has a cannonical ordering of labels
  if (n > 0) {
    pose = l[0];
    if (n > 1) {
      if (hasaffinity) {
        affinity = l[1];
        if (hasrmsd && n > 2)
          rmsd = l[2];
      } else if (hasrmsd) {
        rmsd = l[1];
      }
    }
  }
  labels.push_back(pose);
  affinities.push_back(affinity);
  rmsds.push_back(rmsd);
  seqcont.push_back(seq_continued); //zero for first member of group
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::copyToBlob(Dtype* src, size_t size, Blob<Dtype>* blob,
    bool gpu) {
  Dtype* dst = nullptr;
  if (gpu)
    dst = blob->mutable_gpu_data();
  else
    dst = blob->mutable_cpu_data();
  caffe_copy(size, src, dst);
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::copyToBlobs(const vector<Blob<Dtype>*>& top, bool hasaffinity,
    bool hasrmsd, bool gpu) {
  unsigned idx = 1;
  copyToBlob(&labels[0], labels.size(), top[idx], gpu);
  idx++;
  if (hasaffinity) {
    copyToBlob(&affinities[0], affinities.size(), top[idx], gpu);
    idx++;
  }
  if (hasrmsd) {
    copyToBlob(&rmsds[0], rmsds.size(), top[idx], gpu);
    idx++;
  }

  if(group_size > 1) {
    copyToBlob(&seqcont[0], seqcont.size(), top[idx], gpu);
    idx++;
  }
  if (ligpeturb)
    copyToBlob((Dtype*) &perturbations[0], perturbations.size() * perturbations[0].size(), top.back(), gpu);
}

//set in memory buffer
//will apply translate and rotate iff rotate is valid
template <typename Dtype>
void MolGridDataLayer<Dtype>::setReceptor(const vector<float3>& coords, const vector<smt>& smtypes, const vec& translate, const qt& rotate ) {
  CHECK_GT(batch_info.size(), 0) << "Empty batch info in setReceptor";

  vector<float> types; types.reserve(smtypes.size());
  vector<float> radii; radii.reserve(smtypes.size());

  CHECK_EQ(coords.size(), smtypes.size()) << "Size mismatch between receptor coords and smtypes";


  //receptor atoms
  for (unsigned i = 0, n = smtypes.size(); i < n; i++) {
    smt origt = smtypes[i];
    auto t_r = recTypes->get_int_type(origt);
    int t = t_r.first;
    types.push_back(t);
    radii.push_back(t_r.second);

    if(t < 0 && origt > 1) { //don't warn about hydrogens
      std::cerr << "Unsupported receptor atom type " << GninaIndexTyper::gnina_type_name(origt) << "\n";
    }
  }

  CoordinateSet rec(coords, types, radii, recTypes->num_types());
  if(rotate.real() != 0) {
    //apply transformation
    vec c  = getGridCenter();
    float3 center{c[0],c[1],c[2]};
    float3 trans = {translate[0], translate[1], translate[2]};
    Quaternion Q(rotate.a, rotate.b, rotate.c, rotate.d);
    Transform rectrans(Q, center, trans);

    rectrans.forward(rec, rec);
  }

  batch_info[0].setReceptor(rec);
}

//set in memory buffer, will set grid_Center if it isn't set, but will only overwrite set grid_center if calcCenter
template <typename Dtype>
void MolGridDataLayer<Dtype>::setLigand(const vector<float3>& coords, const vector<smt>& smtypes, bool calcCenter)  {

  CHECK_GT(batch_info.size(), 0) << "Empty batch info in setLigand";

  vector<float> types; types.reserve(coords.size());
  vector<float> radii; radii.reserve(coords.size());

  CHECK_EQ(coords.size(), smtypes.size()) << "Size mismatch between ligand coords and smtypes";

  vec center(0, 0, 0);
  for (unsigned i = 0, n = smtypes.size(); i < n; i++) {
    smt origt = smtypes[i];
    auto t_r = ligTypes->get_int_type(origt);
    int t = t_r.first;
    types.push_back(t);
    radii.push_back(t_r.second);

    if(t < 0 && origt > 1) { //don't warn about hydrogens
      std::cerr << "Unsupported ligand atom type " << GninaIndexTyper::gnina_type_name(origt) << "\n";
    }
  }

  CoordinateSet ligatoms(coords, types, radii, ligTypes->num_types());
  batch_info[0].setLigand(ligatoms);

  if (calcCenter || !isfinite(grid_center[0])) {
    gfloat3 c = batch_info[0].orig_lig_atoms.center();
    setGridCenter(vec(c.x,c.y,c.z));
  }
}

INSTANTIATE_CLASS(MolGridDataLayer);
REGISTER_LAYER_CLASS(MolGridData);

}  // namespace caffe
