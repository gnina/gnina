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

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <boost/timer/timer.hpp>

using namespace libmolgrid;
using namespace std;

namespace caffe {



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
void MolGridDataLayer<Dtype>::getReceptorAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorAtoms";
  CoordinateSet& ratoms = batch_info[batch_idx].rec_atoms;
  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    float4 a; //todo: change to float3
    if(ratoms.type_index[i] >= 0) {
      a.x = ratoms.coord[i][0];
      a.y = ratoms.coord[i][1];
      a.z = ratoms.coord[i][2];
      atoms.push_back(a);
    }
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getLigandAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getLigandAtoms";
  CoordinateSet& latoms = batch_info[batch_idx].lig_atoms;
  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    float4 a; //todo: change to float3
    if(latoms.type_index[i] >= 0) {
      a.x = latoms.coord[i][0];
      a.y = latoms.coord[i][1];
      a.z = latoms.coord[i][2];
      atoms.push_back(a);
    }
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getReceptorChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorChannels";
  CoordinateSet& ratoms = batch_info[batch_idx].rec_atoms;
  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    if(ratoms.type_index[i] >= 0) {
      whichGrid.push_back(ratoms.type_index[i]);
    }
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getLigandChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getLigandChannels";
  CoordinateSet& latoms = batch_info[batch_idx].lig_atoms;
  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    if(latoms.type_index[i] >= 0) {
      whichGrid.push_back(latoms.type_index[i]);
    }
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getReceptorGradient(int batch_idx, vector<gfloat3>& gradient)
{
  gradient.resize(0);
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorGradient";
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CoordinateSet& ratoms = batch_info[batch_idx].rec_atoms;
  vector<gfloat3>& grads = batch_info[batch_idx].rec_gradient;
  CHECK_EQ(ratoms.size(), grads.size());

  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    if(ratoms.type_index[i] >= 0) {
      gradient.push_back(grads[i]);
    }
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

  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getReceptorTransformationGradient";
  CoordinateSet& ratoms = batch_info[batch_idx].rec_atoms;
  vector<gfloat3>& grads = batch_info[batch_idx].rec_gradient;
  CHECK_EQ(ratoms.size(), grads.size());

  gfloat3 tc = batch_info[batch_idx].transform.get_rotation_center();
  vec c(tc.x,tc.y,tc.z);

  for (unsigned i = 0, n = ratoms.size(); i < n; ++i)
  {
    if(ratoms.type_index[i] >= 0) {
      gfloat3 g = grads[i];
      gfloat3 a {ratoms.coord[i][0], ratoms.coord[i][1], ratoms.coord[i][2]};
      vec v(g.x,g.y,g.z);
      vec pos(a.x,a.y,a.z);

      force += v;
      torque += cross_product(pos - c, v);
    }
  }
}


template<typename Dtype>
void MolGridDataLayer<Dtype>::getMappedReceptorGradient(int batch_idx, unordered_map<string, gfloat3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getMappedReceptorGradient";

  CoordinateSet& ratoms = batch_info[batch_idx].rec_atoms;
  vector<gfloat3>& grads = batch_info[batch_idx].rec_gradient;
  CHECK_EQ(ratoms.size(), grads.size());

  for (unsigned i = 0, n = ratoms.size(); i < n; ++i) {
    if(ratoms.type_index[i] >= 0) {
      string xyz = xyz_to_string(ratoms.coord[i][0], ratoms.coord[i][1], ratoms.coord[i][2]);
      gradient[xyz] = grads[i];
    }
  }
}


template<typename Dtype>
void MolGridDataLayer<Dtype>::getLigandGradient(int batch_idx, vector<gfloat3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getLigandGradient";

  gradient.resize(0);
  CoordinateSet& latoms = batch_info[batch_idx].lig_atoms;
  vector<gfloat3>& grads = batch_info[batch_idx].lig_gradient;
  CHECK_EQ(latoms.size(), grads.size());

  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    if(latoms.type_index[i] >= 0) {
      gradient.push_back(grads[i]);
    }
  }
}

template<typename Dtype>
void MolGridDataLayer<Dtype>::getMappedLigandGradient(int batch_idx, unordered_map<string, gfloat3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_info.size()) << "Incorrect batch size in getMappedLigandGradient";

  CoordinateSet& latoms = batch_info[batch_idx].lig_atoms;
  vector<gfloat3>& grads = batch_info[batch_idx].lig_gradient;
  CHECK_EQ(latoms.size(), grads.size());

  for (unsigned i = 0, n = latoms.size(); i < n; ++i) {
    if(latoms.type_index[i] >= 0) {
      string xyz = xyz_to_string(latoms.coord[i][0], latoms.coord[i][1], latoms.coord[i][2]);
      gradient[xyz] = grads[i];
    }
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


//reshape any layer-specific blobs, and set layer-specific dims for any
//universal blobs
template <typename Dtype>
void MolGridDataLayer<Dtype>::setLayerSpecificDims(int number_examples,
    vector<int>& label_shape, const vector<Blob<Dtype>*>& top) {
  top_shape.clear();
  top_shape.push_back(number_examples);
  top_shape.push_back(numchannels);
  top_shape.push_back(dim);
  top_shape.push_back(dim);
  top_shape.push_back(dim);

  example_size = (numchannels)*numgridpoints;
  label_shape.push_back(number_examples); // [batch_size]
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
  ret.max_group_size = param.maxgroupsize();
  ret.group_batch_size = param.maxchunksize();
  ret.cache_structs = param.cache_structs();
  ret.add_hydrogens = param.addh();
  ret.duplicate_first = false;  // at least for now we do any duplication ourselves
  ret.data_root = param.root_folder(); //will have to overwrite if using root_folder2
  ret.recmolcache = param.recmolcache();
  ret.ligmolcache = param.ligmolcache();

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
  double radiusmultiple = param.radius_multiple();
  double fixedradius = param.fixed_radius();
  bool use_covalent_radius = param.use_covalent_radius();
  bool hasaffinity = param.has_affinity();
  bool hasrmsd = param.has_rmsd();
  bool hasgroup = param.maxgroupsize()-1;
  data_ratio = param.source_ratio();
  numposes = param.num_poses();

  if(binary) radiusmultiple = 1.0;
  CHECK_LE(fabs(remainder(dimension,resolution)), 0.001) << "Resolution does not evenly divide dimension.";

  gmaker.initialize(resolution, dimension, radiusmultiple, binary);

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


  if(!inmem)
  {
    const string& source = param.source();
    const string& source2 = param.source2();

    CHECK_GT(source.length(), 0) << "No data source file provided";

    ExampleProviderSettings settings = settings_from_param(param);

    // Read source file(s) with labels and structures,
    // each line is label [group] [affinity] [rmsd] receptor_file ligand_file  (labels not optional)
    data = ExampleProvider(settings, recTypes, ligTypes);
    data.populate(source, 1+hasaffinity+hasrmsd, hasgroup);

    if(source2.length() > 0)
    {
      string root_folder = sanitize_path(param.root_folder());
      string root_folder2 = param.root_folder2();
      if(root_folder2.length() > 0) root_folder2 = sanitize_path(root_folder2);
      else root_folder2 = root_folder; //fall back on first

      CHECK_GE(data_ratio, 0) << "Must provide non-negative ratio for two data sources";
      settings.data_root = root_folder2;
      data2 = ExampleProvider(settings, recTypes, ligTypes);
      data2.populate(source2, 1+hasaffinity+hasrmsd, hasgroup);
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
  batch_info.resize(batch_size * param.maxgroupsize());

  int number_examples = batch_size;
  bool duplicate = this->layer_param_.molgrid_data_param().duplicate_poses();
  if(duplicate) number_examples = batch_size*numposes;
  numchannels = numReceptorTypes+numLigandTypes;
  if(!duplicate && numposes > 1) numchannels = numReceptorTypes+numposes*numLigandTypes;
  vector<int> label_shape;
  setLayerSpecificDims(number_examples, label_shape, top);
  top[0]->Reshape(top_shape);

  // Reshape label, affinity, rmsds
  top[1]->Reshape(label_shape);

  if (hasaffinity)
  {
    top[2]->Reshape(label_shape);
    if (hasrmsd)
    {
      top[3]->Reshape(label_shape);
    }
  }
  else if(hasrmsd)
  {
    top[2]->Reshape(label_shape);
  }

  if(ligpeturb) {
    vector<int> peturbshape(2);
    peturbshape[0] = batch_size;
    peturbshape[1] = output_transform::size(); //trans+orient
    top.back()->Reshape(peturbshape);
  }
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_ex(Dtype *data, const Example& ex,
    typename MolGridDataLayer<Dtype>::mol_info& minfo,
    int pose, output_transform& peturb, bool gpu)
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
  minfo.setReceptor(ex.sets[0]);

  if (doall) { //first ligand pose already set, add others as discrete channels
    CHECK_EQ(pose,0)<< "Invalid pose specifier (internal error)";
    //merge non-receptor poses
    minfo.setLigand(ex.merge_coordinates(1));
  } else {
    //normal single ligand case
    minfo.setLigand(ex.sets[pose+1]);
  }

  set_grid_minfo(data, minfo, peturb, gpu);

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
  c.coord.tocpu();
  for (unsigned i = 0, n = c.size(); i < n; i++) {
    for(unsigned j = 0; j < 3; j++) {
      float diff = jitter * (unit_sample(rng)*2.0-1.0);
      c.coord[i][j] += diff;
    }
  }
}

//take a mol info, which includes receptor and ligand atoms
//and generate the appropriate grids into data
//applies jitter, peturbation, transformations as needed
//transform in minfo will be updated as appropriate
template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data,
    typename MolGridDataLayer<Dtype>::mol_info& minfo,
    output_transform& peturb, bool gpu)
{
  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
  bool fixcenter = param.fix_center_to_origin();

  //figure out transformation

  //what is the rotational center?
  gfloat3 rot_center{0,0,0};
  if (param.use_rec_center())
    rot_center = minfo.rec_atoms.center();
  else
    rot_center = minfo.lig_atoms.center();

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

  minfo.transform = Transform(rot_center, rtranslate, randrotate);

  //copy atoms to void modifying in place
  CoordinateSet& rec_atoms = minfo.rec_atoms;
  CoordinateSet& lig_atoms = minfo.lig_atoms;

  if(!minfo.transform.is_identity()) {
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

  Example ex;
  if(ignore_ligand) {
    ex.sets = {rec_atoms};
  } else {
    ex.sets = {rec_atoms, lig_atoms};
  }
  CoordinateSet atoms = ex.merge_coordinates(); //properly offsets types
  //apply jitter
  apply_jitter(atoms, jitter);

  //set the grid center
  gfloat3 grid_center = rot_center;
  if(fixcenter) {
    grid_center = gfloat3(0,0,0);
  }

  if (atoms.size() == 0) {
    std::cerr << "ERROR: No atoms in molecule.  I can't deal with this.\n";
    exit(-1); //presumably you never actually want this and it results in a cuda error
  }

  //compute grid from atoms
  unsigned dim = gmaker.get_grid_dims().x;
  if (gpu)
  {
    Grid<Dtype, 4, true> outgrid(data, numchannels, dim, dim, dim);
    gmaker.forward(grid_center, atoms, outgrid);
  }
  else
  {
    Grid<Dtype, 4, false> outgrid(data, numchannels, dim,dim,dim);
    gmaker.forward(grid_center, atoms, outgrid);
  }
}


//dump dx files for every atom type, with files names starting with prefix
//only does the very first grid for now
template<typename Dtype>
void MolGridDataLayer<Dtype>::dumpDiffDX(
    const std::string& prefix,
    Blob<Dtype>* top, double scale) const
    {
 abort(); //todo: implement
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
void MolGridDataLayer<Dtype>::dumpGridDX(const std::string& prefix, Dtype* data, double scale) const
{
  abort(); //TODO: implement with grid_io -> molgridder?
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu)
{
  bool hasaffinity = this->layer_param_.molgrid_data_param().has_affinity();
  bool hasrmsd = this->layer_param_.molgrid_data_param().has_rmsd();
  float subgrid_dim = this->layer_param_.molgrid_data_param().subgrid_dim();
  unsigned maxgroupsize = this->layer_param_.molgrid_data_param().maxgroupsize();
  unsigned maxchunksize = this->layer_param_.molgrid_data_param().maxchunksize();

  if ((maxgroupsize >1) && !maxchunksize)
    maxchunksize = maxgroupsize;

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
  if (subgrid_dim || (maxgroupsize-1)) {
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
    CHECK_EQ(maxgroupsize, 1) << "Groups not currently supported with structure in memory";
    if(batch_info[0].rec_atoms.size() == 0) LOG(WARNING) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(batch_info[0].lig_atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(top_data, batch_info[0], peturb, gpu); //TODO how do we know what batch position?
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

    for (int idx = 0; idx < batch_size + (maxchunksize-1)*batch_size; ++idx)
    {
      int batch_idx = idx % batch_size;
      Example ex;
      if (batch_idx < dataswitch)
      {
        data.next(ex);
      }
      else
      {
        data2.next(ex);
      }

      int step = idx / batch_size;
      int offset = ((batch_size * step) + batch_idx) * example_size;

      //if label == -1 then this is a padding example for grouped data; just
      //memset data to 0
      if (ex.labels[0] == -1) {
          unsigned gsize = (numReceptorTypes + numLigandTypes) * dim *
            dim * dim;
          if (gpu) {
            CUDA_CHECK(cudaMemset(top_data+offset, 0, gsize * sizeof(Dtype)));
          }
          else
            memset(top_data+offset, 0, gsize*sizeof(Dtype));
      }

      else {
        if(!duplicate) {
          updateLabels(ex.labels, hasaffinity, hasrmsd);
          set_grid_ex(top_data+offset, ex, batch_info[batch_idx], numposes > 1 ? -1 : 0, peturb, gpu);
          perturbations.push_back(peturb);
        }
        else {
      	  for(unsigned p = 0; p < numposes; p++) {
      	    updateLabels(ex.labels, hasaffinity, hasrmsd);
            int p_offset = batch_idx*(example_size*numposes)+example_size*p;
            set_grid_ex(top_data+p_offset, ex, batch_info[batch_idx], p, peturb, gpu);
            perturbations.push_back(peturb);
          }
        }
        //NOTE: batch_transform contains transformation of last pose only - don't use unless numposes == 1
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
template <typename Dtype>
void MolGridDataLayer<Dtype>::backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu)
{
  abort();
  /*
  //propagate gradient grid onto atom positions
  if(compute_atom_gradients) {
    CHECK(numposes == 1) << "Atomic gradient calculation not supported with numposes != 1";
    unsigned batch_size = top_shape[0];
    Dtype *diff = NULL;
    if(gpu) {
      diff = top[0]->mutable_gpu_diff();
      setAtomGradientsGPU(gmaker, diff, batch_size);
    }
    else {
      diff = top[0]->mutable_cpu_diff();
      for (int item_id = 0; item_id < batch_size; ++item_id) {

        int offset = item_id*example_size;
        typename MolGridDataLayer<Dtype>::mol_transform& transform = batch_transform[item_id];
        gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);
        gmaker.setAtomGradientsCPU(transform.mol.atoms, transform.mol.whichGrid, 
                transform.Q.boost(), diff, transform.mol.gradient, offset, numchannels);
      }
    }
  }
  */
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  abort();
  /*
  typedef typename boost::multi_array_ref<Dtype, 4>  Grids;
  CHECK(numposes == 1) << "Relevance calculations not supported with numposes != 1";

  Dtype *diff = top[0]->mutable_cpu_diff(); //TODO: implement gpu
  Dtype *data = top[0]->mutable_cpu_data();

  //propagate gradient grid onto atom positions
  unsigned batch_size = top_shape[0];
  for (int item_id = 0; item_id < batch_size; ++item_id) {

    int offset = item_id*example_size;
    Grids diffgrids(diff+offset, boost::extents[numchannels][dim][dim][dim]);
    Grids densegrids(data+offset, boost::extents[numchannels][dim][dim][dim]);
    typename MolGridDataLayer<Dtype>::mol_transform& transform = batch_transform[item_id];
    gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);
    gmaker.setAtomRelevanceCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q.boost(),
        densegrids, diffgrids, transform.mol.gradient);
  }
*/
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::clearLabels() {
  labels.clear();
  affinities.clear();
  rmsds.clear();
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::updateLabels(const std::vector<float>& l, bool hasaffinity,
    bool hasrmsd) {
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
  updateLabels(pose, affinity, rmsd);
}

// add labels for a _single_ example (call multiple times with duplicated multiple poses)
template <typename Dtype>
void MolGridDataLayer<Dtype>::updateLabels(Dtype pose, Dtype affinity, Dtype rmsd) {
  labels.push_back(pose);
  affinities.push_back(affinity);
  rmsds.push_back(rmsd);
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
  unsigned rmsdi = 2 + hasaffinity;
  copyToBlob(&labels[0], labels.size(), top[1], gpu);
  if (hasaffinity)
    copyToBlob(&affinities[0], affinities.size(), top[2], gpu);
  if (hasrmsd)
    copyToBlob(&rmsds[0], rmsds.size(), top[rmsdi], gpu);

  if (ligpeturb)
    copyToBlob((Dtype*) &perturbations[0],
        perturbations.size() * perturbations[0].size(), top.back(), gpu);
}

//set in memory buffer
//will apply translate and rotate iff rotate is valid
template <typename Dtype>
void MolGridDataLayer<Dtype>::setReceptor(const vector<atom>& receptor, const vec& translate, const qt& rotate ) {
  CHECK_GT(batch_info.size(), 0) << "Empty batch info in setReceptor";

  vector<float3> cs; cs.reserve(receptor.size());
  vector<float> types; types.reserve(receptor.size());
  vector<float> radii; radii.reserve(receptor.size());

  //receptor atoms
  for (unsigned i = 0, n = receptor.size(); i < n; i++) {
    const atom& a = receptor[i];
    smt t = a.sm;
    auto t_r = recTypes->get_int_type(t);
    if (t_r.first >= 0) {
      const vec& coord = a.coords;
      cs.push_back(gfloat3(coord[0],coord[1],coord[2]));
      types.push_back(t_r.first);
      radii.push_back(t_r.second);
    }
  }

  CoordinateSet rec(cs, types, radii, recTypes->num_types());
  if(rotate.real() != 0) {
    //apply transformation
    vec c  = getCenter();
    float3 center{c[0],c[1],c[2]};
    float3 trans = {translate[0], translate[1], translate[2]};
    Quaternion Q(rotate.a, rotate.b, rotate.c, rotate.d);
    Transform rectrans(Q, center, trans);

    rectrans.forward(rec, rec);
  }

  batch_info[0].setReceptor(rec);
}

//set in memory buffer
template <typename Dtype>
void MolGridDataLayer<Dtype>::setLigand(const vector<atom>& ligand, const vector<vec>& coords, bool calcCenter)  {

  CHECK_GT(batch_info.size(), 0) << "Empty batch info in setLigand";

  vector<float3> cs; cs.reserve(ligand.size());
  vector<float> types; types.reserve(ligand.size());
  vector<float> radii; radii.reserve(ligand.size());

  CHECK_EQ(ligand.size(), coords.size()) << "Size mismatch between ligand atoms and coords";
  //ligand atoms, grid positions offset and coordinates are specified separately
  vec center(0, 0, 0);
  for (unsigned i = 0, n = ligand.size(); i < n; i++) {
    int t = ligand[i].sm;
    auto t_r = ligTypes->get_int_type(t);
    t = t_r.first;
    float r = t_r.second;
    if (t >= 0) {
      const vec& coord = coords[i];
      cs.push_back(gfloat3(coord[0],coord[1],coord[2]));
      types.push_back(t);
      radii.push_back(r);
    }
    else if (t > 1) { //don't warn about hydrogens
      if(t >= ligTypes->num_types()) {
        std::cerr << "Unsupported atom type " << t << " >= " << ligTypes->num_types() << "\n";
      } else {
        std::cerr << "Unsupported atom type " << ligTypes->get_type_names()[t] << "\n";
      }
    }
  }

  CoordinateSet ligatoms(cs, types, radii, ligTypes->num_types());
  batch_info[0].setLigand(ligatoms);

  if (calcCenter || !center_set) {
    gfloat3 c = batch_info[0].lig_atoms.center();
    setCenter(vec(c.x,c.y,c.z));
  }
}

INSTANTIATE_CLASS(MolGridDataLayer);
REGISTER_LAYER_CLASS(MolGridData);

}  // namespace caffe
