#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>
#include <vector>

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

//allocate and initialize atom type data
namespace smina_atom_type
{
  info data[NumTypes] = {{},};

  struct atom_data_initializer {
    atom_data_initializer() {
      for(size_t i = 0u; i < smina_atom_type::NumTypes; ++i)
        smina_atom_type::data[i] = smina_atom_type::default_data[i];
    }
  };

  atom_data_initializer initialize_defaults;
}

namespace caffe {

template <typename Dtype, class GridMakerT>
BaseMolGridDataLayer<Dtype, GridMakerT>::~BaseMolGridDataLayer<Dtype, GridMakerT>() {
  //this->StopInternalThread();

  if(gpu_gridatoms) {
    cudaFree(gpu_gridatoms);
    gpu_gridatoms = NULL;
  }
  if(gpu_gridwhich) {
    cudaFree(gpu_gridwhich);
    gpu_gridwhich = NULL;
  }
}


//for in-memory inputs, set the desired label (for gradient computation)
//note that zero affinity/rmsd means to ignore these
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::setLabels(Dtype pose, Dtype affinity, Dtype rmsd)
{
  clearLabels();
  labels.push_back(pose);
  affinities.push_back(affinity);
  rmsds.push_back(rmsd);
}


//the following really shouldn't be recalculated each evaluation (not including gradients)
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getReceptorAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getReceptorAtoms";
  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
      atoms.push_back(mol.atoms[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getLigandAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getLigandAtoms";
  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
      atoms.push_back(mol.atoms[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getReceptorChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getReceptorChannels";
  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
      whichGrid.push_back(mol.whichGrid[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getLigandChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getLigandChannels";

  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
      whichGrid.push_back(mol.whichGrid[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getReceptorGradient(int batch_idx, vector<float3>& gradient)
{
  gradient.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getReceptorGradient";
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
    {
      gradient.push_back(mol.gradient[i]);
    }
}

/*
 * Compute the transformation gradient of a rigid receptor around the center.
 * The first three numbers are the translation.  The next are the torque.
 */
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getReceptorTransformationGradient(int batch_idx, vec& force, vec& torque)
{
  force = vec(0,0,0);
  torque = vec(0,0,0);

  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getReceptorTransformationGradient";

  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;

  CHECK(mol.center == mem_lig.center) << "Centers not equal; receptor transformation gradient only supported in-mem";


  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
  {
    if (mol.whichGrid[i] < numReceptorTypes)
    {
      float3 g = mol.gradient[i];
      float4 a = mol.atoms[i];
      vec v(g.x,g.y,g.z);
      vec pos(a.x,a.y,a.z);

      force += v;
      torque += cross_product(pos - mol.center, v);
    }
  }
}


template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getMappedReceptorGradient(int batch_idx, unordered_map<string, float3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getMappedReceptorGradient";

  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
  {
    if (mol.whichGrid[i] < numReceptorTypes)
    {
        string xyz = xyz_to_string(mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z);
        gradient[xyz] = mol.gradient[i];
    }
  }
}


template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getLigandGradient(int batch_idx, vector<float3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getLigandGradient";

  gradient.resize(0);
  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
    {
      gradient.push_back(mol.gradient[i]);
    }
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getMappedLigandGradient(int batch_idx, unordered_map<string, float3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getMappedLigandGradient";

  typename MolGridDataLayer<Dtype>::mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
  {
    if (mol.whichGrid[i] >= numReceptorTypes)
    {
        string xyz = xyz_to_string(mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z);
        gradient[xyz] = mol.gradient[i];
    }
  }
}



//ensure gpu memory is of sufficient size
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::allocateGPUMem(unsigned sz)
{
  if(sz > gpu_alloc_size) {
    //deallocate
    if(gpu_gridatoms) {
      cudaFree(gpu_gridatoms);
    }
    if(gpu_gridwhich) {
      cudaFree(gpu_gridwhich);
    }
    //allocate larger
    CUDA_CHECK(cudaMalloc(&gpu_gridatoms, sz*sizeof(float4)));
    CUDA_CHECK(cudaMalloc(&gpu_gridwhich, sz*sizeof(short)));
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
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::setLayerSpecificDims(int number_examples, 
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

template <typename Dtype>
void SubcubeMolGridDataLayer<Dtype>::setLayerSpecificDims(int number_examples, 
    vector<int>& label_shape, const vector<Blob<Dtype>*>& top) {
  //layer shape is TxNx... for RNN
  unsigned grids_per_dim = this->gmaker.grids_per_dim;
  unsigned n_timesteps = grids_per_dim * grids_per_dim * grids_per_dim;

  //setup shape of layer
  this->top_shape.clear();
  if (this->gmaker.stride) {
    this->top_shape.push_back(number_examples);
    this->top_shape.push_back(this->numchannels);
    this->top_shape.push_back(this->dim);
    this->top_shape.push_back(this->dim);
    this->top_shape.push_back(this->dim);

    this->example_size = this->numchannels * this->numgridpoints;
  }
  else {
    this->top_shape.push_back(n_timesteps);
    this->top_shape.push_back(number_examples);
    this->top_shape.push_back(this->numchannels);
    unsigned subgrid_dim_pts = this->gmaker.subgrid_dim_in_points;
    this->top_shape.push_back(subgrid_dim_pts);
    this->top_shape.push_back(subgrid_dim_pts);
    this->top_shape.push_back(subgrid_dim_pts);

    this->example_size = 0;
  }

  label_shape.push_back(n_timesteps);
  label_shape.push_back(number_examples);

  //RNN layer requires a TxN "sequence continuation" blob
  seqcont_shape.clear();
  seqcont_shape.push_back(n_timesteps);
  seqcont_shape.push_back(number_examples);
  int idx = this->ExactNumTopBlobs() - this->ligpeturb - 1;
  top[idx]->Reshape(seqcont_shape);
}

template <typename Dtype>
void GroupedMolGridDataLayer<Dtype>::setLayerSpecificDims(int number_examples, 
    vector<int>& label_shape, const vector<Blob<Dtype>*>& top) {
  this->top_shape.clear();
  this->top_shape.push_back(maxchunksize);
  this->top_shape.push_back(number_examples);
  this->top_shape.push_back(this->numchannels);
  this->top_shape.push_back(this->dim);
  this->top_shape.push_back(this->dim);
  this->top_shape.push_back(this->dim);

  this->example_size = this->numchannels * this->numgridpoints;

  // Reshape label, affinity, rmsds
  label_shape.push_back(maxchunksize);
  label_shape.push_back(batch_size);

  //RNN layer requires a TxN "sequence continuation" blob
  seqcont_shape.clear();
  seqcont_shape.push_back(maxchunksize);
  seqcont_shape.push_back(batch_size);
  int idx = this->ExactNumTopBlobs() - this->ligpeturb - 1;
  top[idx]->Reshape(seqcont_shape);
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

//create the appropriate atom typer given the specified strings
static std::shared_ptr<AtomTyper> create_atom_typer(const string& mapstr, const string& mapfile, const string& defaultmap, bool use_covalent_radius) {

  if(mapstr.size()) {
    //specify map programatically
    stringstream map(mapstr);
    return make_shared<FileMappedGninaTyper>(map, use_covalent_radius);
  }
  else if(mapfile.size()) {
    return make_shared<FileMappedGninaTyper>(mapfile, use_covalent_radius);
  }
  else {
    stringstream map(defaultmap.c_str());
    return make_shared<FileMappedGninaTyper>(map, use_covalent_radius);
  }
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
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {

  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();

  root_folder = param.root_folder();
  num_rotations = param.rotate();
  inmem = param.inmemory();
  dimension = param.dimension();
  resolution = param.resolution();
  binary = param.binary_occupancy();
  randtranslate = param.random_translate();
  randrotate = param.random_rotation();
  ligpeturb = param.peturb_ligand();
  ligpeturb_translate = param.peturb_ligand_translate();
  ligpeturb_rotate = param.peturb_ligand_rotate();
  jitter = param.jitter();
  ignore_ligand = param.ignore_ligand();
  radiusmultiple = param.radius_multiple();
  fixedradius = param.fixed_radius();
  use_covalent_radius = param.use_covalent_radius();
  bool hasaffinity = param.has_affinity();
  bool hasrmsd = param.has_rmsd();
  bool hasgroup = param.maxgroupsize()-1;
  data_ratio = param.source_ratio();
  root_folder2 = param.root_folder2();
  numposes = param.num_poses();

  if(binary) radiusmultiple = 1.0;
  CHECK_LE(fabs(remainder(dimension,resolution)), 0.001) << "Resolution does not evenly divide dimension.";

  gmaker.initialize(param);

  dim = round(dimension/resolution)+1; //number of grid points on a side
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

  recTypes = create_atom_typer(recmapstr, recmapfile, default_recmap, use_covalent_radius);
  ligTypes = create_atom_typer(ligmapstr, ligmapfile, default_ligmap, use_covalent_radius);

  numReceptorTypes = recTypes->num_types();
  numLigandTypes = ligTypes->num_types();


  if(!inmem)
  {
    const string& source = param.source();
    const string& source2 = param.source2();
    root_folder = sanitize_path(param.root_folder());
    root_folder2 = param.root_folder2();
    if(root_folder2.length() > 0) root_folder2 = sanitize_path(root_folder2);
    else root_folder2 = root_folder; //fall back on first

    CHECK_GT(source.length(), 0) << "No data source file provided";

    ExampleProviderSettings settings = settings_from_param(param);

    // Read source file(s) with labels and structures,
    // each line is label [group] [affinity] [rmsd] receptor_file ligand_file  (labels not optional)
    data = ExampleProvider(settings, recTypes, ligTypes);
    data.populate(source, 1+hasaffinity+hasrmsd, hasgroup);

    if(source2.length() > 0)
    {
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
  batch_transform.resize(batch_size * param.maxgroupsize());

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

//return quaternion representing one of 24 distinct axial rotations
template <typename Dtype, class GridMakerT>
typename MolGridDataLayer<Dtype>::quaternion BaseMolGridDataLayer<Dtype, GridMakerT>::axial_quaternion()
{
  unsigned rot = current_rotation;
  qt ret;
  //first rotate to a face
  switch(rot%6) {
    case 0:
      ret = typename MolGridDataLayer<Dtype>::quaternion(1,0,0,0); //identity
      break;
    case 1: //rotate z 90
      ret = typename MolGridDataLayer<Dtype>::quaternion(sqrt(0.5),0,0,sqrt(0.5));
      break;
    case 2: //rotate z 180
      ret = typename MolGridDataLayer<Dtype>::quaternion(0,0,0,1);
      break;
    case 3: //rotate z 270 (-90)
      ret = typename MolGridDataLayer<Dtype>::quaternion(sqrt(0.5),0,0,-sqrt(0.5));
      break;
    case 4: //rotate y 90
      ret = typename MolGridDataLayer<Dtype>::quaternion(sqrt(0.5),0,sqrt(0.5),0);
      break;
    case 5: //rotate y -90
      ret = typename MolGridDataLayer<Dtype>::quaternion(sqrt(0.5),0,-sqrt(0.5),0);
      break;
  }

  //now four rotations around x axis
  rot /= 6;
  switch(rot%4) {
    case 0:
      break; //identity
    case 1: //90 degrees
      ret *= typename MolGridDataLayer<Dtype>::quaternion(sqrt(0.5),sqrt(0.5),0,0);
      break;
    case 2: //180
      ret *= typename MolGridDataLayer<Dtype>::quaternion(0,1,0,0);
      break;
    case 3: //270
      ret *= typename MolGridDataLayer<Dtype>::quaternion(sqrt(0.5),-sqrt(0.5),0,0);
      break;
  }
  return ret;
}


template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::set_mol_info(const CoordinateSet& c,
    unsigned mapoffset, typename MolGridDataLayer<Dtype>::mol_info& minfo)
{
  minfo.atoms.clear();
  minfo.whichGrid.clear();
  minfo.gradient.clear();

  unsigned N = c.size();
  minfo.atoms.reserve(N);
  minfo.whichGrid.reserve(N);
  minfo.gradient.reserve(N);
  unsigned cnt = 0;
  vec center(0,0,0);

  for(unsigned i = 0; i < N; i++) {
    float4 atom;
    atom.x = c.coord[i][0];
    atom.y = c.coord[i][1];
    atom.z = c.coord[i][2];
    atom.w = c.radius[i];

    int t = c.type_index[i];
    if(t >= 0) {
      cnt++;
      minfo.atoms.push_back(atom);
      minfo.whichGrid.push_back(t);
      minfo.gradient.push_back(make_float3(0,0,0));
      center += vec(atom.x,atom.y,atom.z);
    }
  }

  if(cnt == 0) {
    if(c.src && c.src != string("none")) std::cerr << "WARNING: No atoms in " << c.src <<"\n";
  }
  else {
    center /= cnt;
  }
  minfo.center = center;

}

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::set_grid_ex(Dtype *data, const Example& ex,
    typename MolGridDataLayer<Dtype>::mol_transform& transform,
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

  typename MolGridDataLayer<Dtype>::mol_info rec;
  typename MolGridDataLayer<Dtype>::mol_info lig;

  CHECK_LT(pose, ex.sets.size()-1) << "Incorrect pose index";  //pose index doesn't include rec

  set_mol_info(ex.sets[0], 0, rec);  //could be "none"
  set_mol_info(ex.sets[pose+1], numReceptorTypes, lig);

  if(doall) { //first ligand pose already set, add others as discrete channels
      for(unsigned p = 1, np = ex.sets.size()-1; p < np; p++) {
        typename MolGridDataLayer<Dtype>::mol_info tmplig;
        set_mol_info(ex.sets[p+1], numReceptorTypes+numLigandTypes*p, tmplig);
        lig.append(tmplig);
      }
  }

  set_grid_minfo(data, rec, lig, transform, peturb, gpu);

}


template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::set_grid_minfo(Dtype *data, 
    const typename MolGridDataLayer<Dtype>::mol_info& recatoms,
    const typename MolGridDataLayer<Dtype>::mol_info& ligatoms, 
    typename MolGridDataLayer<Dtype>::mol_transform& transform,
    output_transform& peturb, bool gpu)
{
  bool fixcenter = this->layer_param_.molgrid_data_param().fix_center_to_origin();
  //set grid values from mol info
  //first clear transform from the previous batch
  rng_t* rng = caffe_rng();
  transform = typename MolGridDataLayer<Dtype>::mol_transform();
  typename MolGridDataLayer<Dtype>::mol_transform ligtrans;

  //figure out transformation
  //note there are currently two code paths - one where setAtoms performs the transformation
  //and one where the transformation is applied here; the first is faster, as it can be done
  //on the GPU, but the second let's us support "weird" transformation like peturbations
  //in the first case, the coordinates of the mol don't change; in the other they are mogrified
  //TODO: unify this - I think transformation should be treated separately from gridding
  //I don't think random rotation and backwards gradients are currently working
  transform.Q = typename MolGridDataLayer<Dtype>::quaternion(1, 0, 0, 0);

  if (current_rotation == 0 && !randrotate)
    transform.Q = typename MolGridDataLayer<Dtype>::quaternion(1, 0, 0, 0); //check real part to avoid mult

  if (randrotate)
  {
    transform.set_random_quaternion(rng);
  }

  if (randtranslate)
  {
    double radius = ligatoms.radius();
    //don't let ligand atoms translate out of sphere inscribed in box
    if (ignore_ligand) radius = 0;
    double maxtrans = max(dimension / 2.0 - radius, 0.0);
    transform.add_random_displacement(rng, min(randtranslate, maxtrans));
  }

  if (current_rotation > 0) {
    transform.Q *= axial_quaternion();
  }

  //include receptor and ligand atoms
  transform.mol.append(recatoms);
  //unless otherwise specified, set rotation center and grid center to ligand center
  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
  if (param.use_rec_center())
    transform.mol.center = recatoms.center;
  else
    transform.mol.center = ligatoms.center;

  //GPU transformation will use Q and grid_center to apply transformation
  typename MolGridDataLayer<Dtype>::quaternion Q = transform.Q;
  vec grid_center(0, 0, 0);

  if(!fixcenter) {
  //center on ligand, offset by random translate
    if (param.use_rec_center())
      transform.center += recatoms.center;
    else
      transform.center += ligatoms.center;
    grid_center = transform.center;
  }

  //both ligand_peturbation and zero center need to apply mol transform here,
  //since they don't fit the current framework for GPU transformation
  //TODO move this into gridmaker.setAtoms, have it take separate receptor and ligand transformations (or have addAtoms)
  typename MolGridDataLayer<Dtype>::mol_info ligmol = ligatoms;
  if (fixcenter || ligpeturb) {
    transform.mol.apply_transform(transform); //mogrify the coordinates -- these are recatoms only, rotate around ligand center
    ligmol.apply_transform(transform); //modify ligand coordinates
    //Q is already applied
    Q = qt(1, 0, 0, 0);
    grid_center = vec(0,0,0);
  }

  //add ligatoms to transform.mol
  if (ligpeturb) {
    if (ligpeturb_rotate)
    {
      ligtrans.set_random_quaternion(rng);
    }
    else
    {
      ligtrans.Q = typename MolGridDataLayer<Dtype>::quaternion(1, 0, 0, 0); //identity
    }
    ligtrans.add_random_displacement(rng, ligpeturb_translate);
    ligmol.apply_transform(ligtrans); //peturb
    transform.mol.append(ligmol); //append

    //store the inverse transformation
    peturb.x = ligtrans.center[0];
    peturb.y = ligtrans.center[1];
    peturb.z = ligtrans.center[2];

    qt qinv = conj(ligtrans.Q) / norm(ligtrans.Q); //not Cayley, not euclidean norm - already squared
    peturb.set_from_quaternion(qinv);

    //set the center to the translated value
    transform.mol.center = ligmol.center + ligtrans.center;
  } else if (ignore_ligand) {
    //do nothing - ligand is only used to set center
  } else {
    transform.mol.append(ligmol);
  }


  //with fix_center, this should be zero
  gmaker.setCenter(grid_center[0], grid_center[1], grid_center[2]);

  if (transform.mol.atoms.size() == 0) {
    std::cerr << "ERROR: No atoms in molecule.  I can't deal with this.\n";
    exit(-1); //presumably you never actually want this and it results in a cuda error
  }
  if (jitter > 0) {
    //add small random displacement (in-place) to atoms
    for (unsigned i = 0, n = transform.mol.atoms.size(); i < n; i++) {
      float4& atom = transform.mol.atoms[i];
      float xdiff = jitter * (unit_sample(rng) * 2.0 - 1.0);
      atom.x += xdiff;
      float ydiff = jitter * (unit_sample(rng) * 2.0 - 1.0);
      atom.y += ydiff;
      float zdiff = jitter * (unit_sample(rng) * 2.0 - 1.0);
      atom.z += zdiff;
    }
  }

  //compute grid from atom info arrays
  if (gpu)
  {
    unsigned natoms = transform.mol.atoms.size();
    allocateGPUMem(natoms);
    CUDA_CHECK(cudaMemcpy(gpu_gridatoms, &transform.mol.atoms[0], natoms*sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gpu_gridwhich, &transform.mol.whichGrid[0], natoms*sizeof(short), cudaMemcpyHostToDevice));

    gmaker.template setAtomsGPU<Dtype>(natoms, gpu_gridatoms, gpu_gridwhich, Q, numchannels, data);
  }
  else
  {
    gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, Q.boost(), data, numchannels);
  }
}

template <typename Dtype>
void GroupedMolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data, 
    const typename MolGridDataLayer<Dtype>::mol_info& recatoms,
    const typename MolGridDataLayer<Dtype>::mol_info& ligatoms, 
    typename MolGridDataLayer<Dtype>::mol_transform& transform,
    typename BaseMolGridDataLayer<Dtype, GridMaker>::output_transform& peturb, bool gpu)
{
  //if it's the first frame, call the base class method
  if (example_idx < batch_size) {
    BaseMolGridDataLayer<Dtype, GridMaker>::set_grid_minfo(data, recatoms, ligatoms, 
        transform, peturb, gpu);
  }
  //otherwise use the same transformation from frame 0
  else {
    rng_t* rng = caffe_rng();
    bool fixcenter = this->layer_param_.molgrid_data_param().fix_center_to_origin();
    unsigned transform_idx = example_idx % batch_size;
    transform.Q = this->batch_transform[transform_idx].Q;
    transform.center = translations[transform_idx];
    typename MolGridDataLayer<Dtype>::mol_transform ligtrans;

    //include receptor and ligand atoms
    transform.mol.append(recatoms);
    //unless otherwise specified, set rotation center and grid center to ligand center
    const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
    if (param.use_rec_center())
      transform.mol.center = recatoms.center;
    else
      transform.mol.center = ligatoms.center;

    //GPU transformation will use Q and grid_center to apply transformation
    typename MolGridDataLayer<Dtype>::quaternion Q = transform.Q;
    vec grid_center(0, 0, 0);

    if(!fixcenter) {
    //center on ligand, offset by random translate
      if (param.use_rec_center())
        transform.center += recatoms.center;
      else
        transform.center += ligatoms.center;
      grid_center = transform.center;
    }

    typename MolGridDataLayer<Dtype>::mol_info ligmol = ligatoms;
    if(this->ligpeturb || fixcenter) {
      transform.mol.apply_transform(transform);
      ligmol.apply_transform(transform);

      Q = qt(1, 0, 0, 0);
      grid_center = vec(0,0,0);
    }
    if (this->ligpeturb) {
      ligtrans.center[0] = this->perturbations[example_idx % batch_size].x;
      ligtrans.center[1] = this->perturbations[example_idx % batch_size].y;
      ligtrans.center[2] = this->perturbations[example_idx % batch_size].z;
      typename BaseMolGridDataLayer<Dtype, GridMaker>::output_transform& group_transform = this->perturbations[example_idx % batch_size];
      ligtrans.Q = typename MolGridDataLayer<Dtype>::quaternion(group_transform.a, group_transform.b, 
                              group_transform.c, group_transform.d);
      ligmol.apply_transform(ligtrans);
      transform.mol.append(ligmol);

      //store the inverse transformation
      peturb = this->perturbations[example_idx % batch_size];

      //set the center to the translated value
      transform.mol.center = ligmol.center + ligtrans.center;
    } else if(this->ignore_ligand) {
      //do nothing - ligand is only used to set center
    } else {
      transform.mol.append(ligmol);
    }

    transform.center += transform.mol.center;

    //TODO move this into gridmaker.setAtoms, have it just take the mol_transform as input - separate receptor transform as well
    this->gmaker.setCenter(grid_center[0], grid_center[1], grid_center[2]);

    if(transform.mol.atoms.size() == 0) {
       std::cerr << "ERROR: No atoms in molecule.  I can't deal with this.\n";
       exit(-1); //presumably you never actually want this and it results in a cuda error
    } 
    if(this->jitter > 0) {
      //add small random displacement (in-place) to atoms
      for(unsigned i = 0, n = transform.mol.atoms.size(); i < n; i++) {
        float4& atom = transform.mol.atoms[i];
        float xdiff = this->jitter*(unit_sample(rng)*2.0-1.0);
        atom.x += xdiff;
        float ydiff = this->jitter*(unit_sample(rng)*2.0-1.0);
        atom.y += ydiff;
        float zdiff = this->jitter*(unit_sample(rng)*2.0-1.0);
        atom.z += zdiff;
      }
    }

    //compute grid from atom info arrays
    if(gpu)
    {
      unsigned natoms = transform.mol.atoms.size();
      BaseMolGridDataLayer<Dtype, GridMaker>::allocateGPUMem(natoms);
      CUDA_CHECK(cudaMemcpy(this->gpu_gridatoms, &transform.mol.atoms[0], natoms*sizeof(float4), cudaMemcpyHostToDevice));
      CUDA_CHECK(cudaMemcpy(this->gpu_gridwhich, &transform.mol.whichGrid[0], natoms*sizeof(short), cudaMemcpyHostToDevice));

      this->gmaker.template setAtomsGPU<Dtype>(natoms, this->gpu_gridatoms, this->gpu_gridwhich, transform.Q, this->numReceptorTypes+this->numLigandTypes, data);
    }
    else
    {
      this->gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q.boost(), data, this->numReceptorTypes + this->numLigandTypes);
    }
  }
  //either way, done with this example
  ++example_idx;
  example_idx = example_idx % (maxgroupsize * batch_size);
}

//return a string representation of the atom type(s) represented by index
//in map - this isn't particularly efficient, but is only for debug purposes
template <typename Dtype, class GridMakerT>
string BaseMolGridDataLayer<Dtype, GridMakerT>::getIndexName(const vector<int>& map, unsigned index) const
		{
	stringstream ret;
	stringstream altret;
	for (unsigned at = 0; at < smina_atom_type::NumTypes; at++)
	{
		if (map[at] == index)
		{
			ret << smina_type_to_string((smt) at);
			altret << "_" << at;
		}
	}

	if (ret.str().length() > 32) //there are limits on file name lengths
		return altret.str();
	else
		return ret.str();
}


//output a grid the file in dx format (for debug)
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::outputDXGrid(std::ostream& out, Grids& grid, unsigned g, double scale, unsigned n) const
{
  out.precision(5);
  setprecision(5);
  out << fixed;
  out << "object 1 class gridpositions counts " << n << " " << n << " " << " " << n << "\n";
  out << "origin";
  for (unsigned i = 0; i < 3; i++)
  {
    out << " " << mem_lig.center[i]-dimension/2.0;
  }
  out << "\n";
  out << "delta " << resolution << " 0 0\ndelta 0 " << resolution << " 0\ndelta 0 0 " << resolution << "\n";
  out << "object 2 class gridconnections counts " << n << " " << n << " " << " " << n << "\n";
  out << "object 3 class array type double rank 0 items [ " << n*n*n << "] data follows\n";
  //now coordinates - x,y,z
  out << scientific;
  out.precision(6);
  unsigned total = 0;
  for (unsigned i = 0; i < n; i++)
  {
    for (unsigned j = 0; j < n; j++)
    {
      for (unsigned k = 0; k < n; k++)
      {
        out << grid[g][i][j][k]*scale;
        total++;
        if(total % 3 == 0) out << "\n";
        else out << " ";
      }
    }
  }
}

//dump dx files for every atom type, with files names starting with prefix
//only does the very first grid for now
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::dumpDiffDX(const std::string& prefix,
		Blob<Dtype>* top, double scale) const
{
  Dtype* diff = top->mutable_cpu_diff();
	Grids grids(diff,
			boost::extents[numReceptorTypes + numLigandTypes][dim][dim][dim]);
    CHECK_GT(mem_lig.atoms.size(),0) << "DX dump only works with in-memory ligand";
    CHECK_EQ(randrotate, false) << "DX dump requires no rotation";
    CHECK_LE(numposes, 1) << "DX dump requires numposes == 1";
	for (unsigned a = 0, na = numReceptorTypes; a < na; a++) {
		string name = getIndexName(rmap, a);
		string fname = prefix + "_rec_" + name + ".dx";
		ofstream out(fname.c_str());
		outputDXGrid(out, grids, a, scale, dim);
	}
	for (unsigned a = 0, na = numLigandTypes; a < na; a++) {
			string name = getIndexName(lmap, a);
			string fname = prefix + "_lig_" + name + ".dx";
			ofstream out(fname.c_str());
			outputDXGrid(out, grids, numReceptorTypes+a, scale, dim);
	}

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

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top)
{
  forward(bottom, top, false);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype,GridMakerT>::dumpGridDX(const std::string& prefix, Dtype* data, double scale) const
{
  typedef boost::multi_array_types::index_range range_t;
  typename Grids::index_gen indices;
  Grids grid(data,
			boost::extents[numReceptorTypes + numLigandTypes][dim][dim][dim]);
	for (unsigned a = 0, na = numReceptorTypes; a < na; a++) {
    if (!(grid_empty(grid[indices[a]
        [range_t(0, dim)] [range_t(0, dim)][range_t(0,dim)]].origin(), dim*dim*dim))) {
		  string name = getIndexName(rmap, a);
		  string fname = prefix + "_rec_" + name + to_string(current_iter) + ".dx";
		  ofstream out(fname.c_str());
		  outputDXGrid(out, grid, a, scale, dim);
    }
	}
	for (unsigned a = 0, na = numLigandTypes; a < na; a++) {
    if (!(grid_empty(grid[indices[numReceptorTypes+a]
        [range_t(0, dim)] [range_t(0, dim)][range_t(0,dim)]].origin(), dim*dim*dim))) {
			string name = getIndexName(lmap, a);
			string fname = prefix + "_lig_" + name + to_string(current_iter) + ".dx";
			ofstream out(fname.c_str());
			outputDXGrid(out, grid, numReceptorTypes+a, scale, dim);
    }
	}
}

//if doing subcubes, output a separate file for each subcube (for now) to
//confirm that they look reasonable
template<typename Dtype>
void SubcubeMolGridDataLayer<Dtype>::dumpDiffDX(const std::string& prefix,
		Blob<Dtype>* top, double scale) const
{
  Dtype* diff = top->mutable_cpu_diff();
  CHECK_GT(this->mem_lig.atoms.size(),0) << "DX dump only works with in-memory ligand";
  CHECK_EQ(this->randrotate, false) << "DX dump requires no rotation";
  unsigned ncubes = this->gmaker.grids_per_dim * this->gmaker.grids_per_dim * 
    this->gmaker.grids_per_dim;
  const unsigned& subgrid_dim_in_points = this->gmaker.subgrid_dim_in_points;
  unsigned instance_size = (this->numReceptorTypes + this->numLigandTypes) * subgrid_dim_in_points * subgrid_dim_in_points * subgrid_dim_in_points;
	for (unsigned a = 0, na = this->numReceptorTypes; a < na; a++) {
		string name = this->getIndexName(this->rmap, a);
    for (size_t i=0; i<ncubes; ++i) {
		  string fname = prefix + "_cube" + std::to_string(i) + "_rec_" + name + ".dx";
		  ofstream out(fname.c_str());
      unsigned offset = i * instance_size; 
      typename BaseMolGridDataLayer<Dtype, SubcubeGridMaker>::Grids subgrids(diff+offset, boost::extents[this->numReceptorTypes+this->numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);
		  this->outputDXGrid(out, subgrids, a, scale, subgrid_dim_in_points);
      }
	}
	for (unsigned a = 0, na = this->numLigandTypes; a < na; a++) {
			string name = this->getIndexName(this->lmap, a);
        for (size_t i=0; i<ncubes; ++i) {
		      string fname = prefix + "_cube" + std::to_string(i) + "_lig_" + name + ".dx";
		      ofstream out(fname.c_str());
          unsigned offset = i * instance_size; 
          typename BaseMolGridDataLayer<Dtype, SubcubeGridMaker>::Grids subgrids(diff+offset, boost::extents[this->numReceptorTypes+this->numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);
		      this->outputDXGrid(out, subgrids, this->numReceptorTypes+a, scale, 
              subgrid_dim_in_points);
        }
	}

}

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu)
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
    CHECK_EQ(maxgroupsize, 1) << "Groups not currently supported with structure in memory";
    if(mem_rec.atoms.size() == 0) LOG(WARNING) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(mem_lig.atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(top_data, mem_rec, mem_lig, batch_transform[0], peturb, gpu); //TODO how do we know what batch position?
    perturbations.push_back(peturb);

    if (num_rotations > 0) {
      current_rotation = (current_rotation+1)%num_rotations;
    }

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

      updateLabels(ex.labels, hasaffinity, hasrmsd);
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
          set_grid_ex(top_data+offset, ex, batch_transform[batch_idx], numposes > 1 ? -1 : 0, peturb, gpu);
          perturbations.push_back(peturb);
        }
        else {
      	  for(unsigned p = 0; p < numposes; p++) {
            int p_offset = batch_idx*(example_size*numposes)+example_size*p;
            set_grid_ex(top_data+p_offset, ex, batch_transform[batch_idx], p, peturb, gpu);
            perturbations.push_back(peturb);
            //NOTE: num_rotations not actually implemented!
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

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  backward(top, bottom, false);
}


/* backpropagates gradients onto atoms (note there is not actual bottom)
 * Only performed when compute_atom_gradients is true.
 */
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu)
{
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
  ++current_iter;
}

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{

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

  //float bottom_sum = 0.0;
  //for(int i = 0; i < bottom[0]->count(); i++)
  //{
  //        bottom_sum += bottom[0]->cpu_diff()[i];
  //}
  //std::cout << "MOLGRID BOTTOM: " << bottom_sum << '\n';


}

template <typename Dtype>
shared_ptr<Layer<Dtype> > GetMolGridDataLayer(const LayerParameter& param) {
  const MolGridDataParameter& mgrid_param = param.molgrid_data_param();
  if (mgrid_param.subgrid_dim()) {
    return shared_ptr<Layer<Dtype> >(new SubcubeMolGridDataLayer<Dtype>(param));
  }
  else if(mgrid_param.maxgroupsize() > 1) {
    return shared_ptr<Layer<Dtype> >(new GroupedMolGridDataLayer<Dtype>(param));
  }
  else
    return shared_ptr<Layer<Dtype> >(new GenericMolGridDataLayer<Dtype>(param));
}

INSTANTIATE_CLASS(GroupedMolGridDataLayer);
INSTANTIATE_CLASS(SubcubeMolGridDataLayer);
INSTANTIATE_CLASS(GenericMolGridDataLayer);
REGISTER_LAYER_CREATOR(MolGridData, GetMolGridDataLayer);

}  // namespace caffe
