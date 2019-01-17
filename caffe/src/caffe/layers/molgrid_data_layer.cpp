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

  if(data) delete data;
  if(data2) delete data2;
}

template <typename Dtype, class GridMakerT>
BaseMolGridDataLayer<Dtype, GridMakerT>::example::example(BaseMolGridDataLayer<Dtype, GridMakerT>::string_cache& cache, string line, const MolGridDataParameter& param)
  : label(0), affinity(0.0), rmsd(0.0), affinity_weight(1.0), group(-1)
{
  stringstream stream(line);
  string tmp;

  bool hasaffinity = param.has_affinity();
  bool hasrmsd = param.has_rmsd();
  bool hasgroup = param.maxgroupsize() - 1;
  unsigned numposes = param.num_poses();
  double affinity_reweight_mean = param.affinity_reweight_mean();
  double affinity_reweight_std = param.affinity_reweight_std();
  double affinity_reweight_stdcut = param.affinity_reweight_stdcut();

  //first the label
  stream >> label;
  if(hasaffinity)
   stream >> affinity;
  if(hasrmsd)
   stream >> rmsd;
  if(hasgroup) {
    stream >> group;
  }
  //receptor
  stream >> tmp;
  CHECK(tmp.length() > 0) << "Empty receptor, missing affinity/rmsd? Line:\n" << line;
  receptor = cache.get(tmp);
  //ligand(s)

  for(unsigned i = 0; i < numposes; i++) {
    tmp.clear();
    stream >> tmp;
    CHECK(tmp.length() > 0) << "Empty ligand, missing affinity/rmsd? Line:\n" << line;
    ligands.push_back(cache.get(tmp));
  }

  if(affinity_reweight_stdcut > 0 && affinity != 0) {
    //weight the affinities inversely to a normal distribution
    double x = fabs(fabs(affinity)-affinity_reweight_mean);
    x = min(x,affinity_reweight_stdcut*affinity_reweight_std);
    x = x*x; //square, but no negative since we want the inverse
    affinity_weight = exp(x/(2.0*affinity_reweight_std*affinity_reweight_std));
  }

}


//for in-memory inputs, set the desired label (for gradient computation)
//note that zero affinity/rmsd means to ignore these
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::setLabels(Dtype pose, Dtype affinity, Dtype rmsd)
{
  clearLabels();
  updateLabels(pose, affinity, rmsd);
}


//the following really shouldn't be recalculated each evaluation (not including gradients)
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getReceptorAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getReceptorAtoms";
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
      atoms.push_back(mol.atoms[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getLigandAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getLigandAtoms";
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
      atoms.push_back(mol.atoms[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getReceptorChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getReceptorChannels";
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
      whichGrid.push_back(mol.whichGrid[i]);
}

template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::getLigandChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  CHECK_LT(batch_idx,batch_transform.size()) << "Incorrect batch size in getLigandChannels";

  mol_info& mol = batch_transform[batch_idx].mol;
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
  mol_info& mol = batch_transform[batch_idx].mol;
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

  mol_info& mol = batch_transform[batch_idx].mol;

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

  mol_info& mol = batch_transform[batch_idx].mol;
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
  mol_info& mol = batch_transform[batch_idx].mol;
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

  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
  {
    if (mol.whichGrid[i] >= numReceptorTypes)
    {
        string xyz = xyz_to_string(mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z);
        gradient[xyz] = mol.gradient[i];
    }
  }
}


//modify examples to remove any without both actives an inactives
//factored this into its own function due to the need to fully specialize setup below
template<typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::remove_missing_and_setup(vector<typename BaseMolGridDataLayer<Dtype, GridMakerT>::balanced_example_provider>& examples)
{
  vector<balanced_example_provider> tmp;
  for(unsigned i = 0, n = examples.size(); i < n; i++)
  {
    if(examples[i].num_actives() > 0 && examples[i].num_decoys() > 0) {
      //eliminate empty buckets
      tmp.push_back(examples[i]);
      tmp.back().setup();
    }
    else if(examples[i].num_actives() > 0)
    {
      example tmp;
      examples[i].next_active(tmp);
      LOG(INFO) << "Dropping receptor " << tmp.receptor << " with no decoys.";
    }
    else if(examples[i].num_decoys() > 0)
    {
      example tmp;
      examples[i].next_decoy(tmp);
      LOG(INFO) << "Dropping receptor " << tmp.receptor << " with no actives.";
    }
  }

  swap(examples,tmp);
}

//specialized version for balanced data that remove receptors without any actives or decoys
//annoyingly, have to specialize Dtype and GridMakerT - TODO: replace with a macro
template <>
template <>
void BaseMolGridDataLayer<float, GridMaker>::receptor_stratified_example_provider<typename BaseMolGridDataLayer<float, GridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template<>
template<>
void BaseMolGridDataLayer<double, GridMaker>::receptor_stratified_example_provider<typename BaseMolGridDataLayer<double, GridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template <>
template <>
void BaseMolGridDataLayer<float, SubcubeGridMaker>::receptor_stratified_example_provider<typename BaseMolGridDataLayer<float, SubcubeGridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template<>
template<>
void BaseMolGridDataLayer<double, SubcubeGridMaker>::receptor_stratified_example_provider<typename BaseMolGridDataLayer<double, SubcubeGridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
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

//allocate and return an example provider to the specifications of the parm object
template <typename Dtype, class GridMakerT>
typename BaseMolGridDataLayer<Dtype, GridMakerT>::example_provider* BaseMolGridDataLayer<Dtype, GridMakerT>::create_example_data(const MolGridDataParameter& parm)
{
  bool balanced  = parm.balanced();
  bool strat_receptor  = parm.stratify_receptor();
  bool strat_aff = parm.stratify_affinity_max() != parm.stratify_affinity_min();
  bool grouped = parm.maxgroupsize()-1;

  //strat_aff > strat_receptor > balanced
  if(strat_aff)
  {
    if(strat_receptor)
    {
      if(balanced) // sample 2 from each receptor
      {
        if(grouped)
          return new grouped_example_provider<affinity_stratified_example_provider<receptor_stratified_example_provider<balanced_example_provider, 2> > >(parm);
        else
          return new affinity_stratified_example_provider<receptor_stratified_example_provider<balanced_example_provider, 2> >(parm);
      }
      else //sample 1 from each receptor
      {
        if(grouped)
          return new grouped_example_provider<affinity_stratified_example_provider<receptor_stratified_example_provider<uniform_example_provider, 1> > >(parm);
        else
          return new affinity_stratified_example_provider<receptor_stratified_example_provider<uniform_example_provider, 1> >(parm);
      }
    }
    else
    {
      if(balanced)
      {
        if(grouped)
          return new grouped_example_provider<affinity_stratified_example_provider<balanced_example_provider> >(parm);
        else
          return new affinity_stratified_example_provider<balanced_example_provider>(parm);
      }
      else //sample 1 from each receptor
      {
        if(grouped)
          return new grouped_example_provider<affinity_stratified_example_provider<uniform_example_provider> >(parm);
        else
          return new affinity_stratified_example_provider<uniform_example_provider>(parm);
      }
    }
  }
  else if(strat_receptor)
  {
    if(balanced) // sample 2 from each receptor
    {
      if(grouped)
        return new grouped_example_provider<receptor_stratified_example_provider<balanced_example_provider, 2> >(parm);
      else
        return new receptor_stratified_example_provider<balanced_example_provider, 2>(parm);
    }
    else //sample 1 from each receptor
    {
      if(balanced)
        return new grouped_example_provider<receptor_stratified_example_provider<uniform_example_provider, 1> >(parm);
      else
        return new receptor_stratified_example_provider<uniform_example_provider, 1>(parm);
    }
  }
  else if(balanced)
  {
    if(grouped)
      return new grouped_example_provider<balanced_example_provider>(parm);
    else
      return new balanced_example_provider(parm);
  }
  else
  {
    if(grouped)
      return new grouped_example_provider<uniform_example_provider>(parm);
    else
      return new uniform_example_provider(parm);
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

//fill in training examples
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::populate_data(const string& root_folder, const string& source,
    BaseMolGridDataLayer<Dtype, GridMakerT>::example_provider* data, bool hasaffinity, bool hasrmsd, bool hasgroup)
{
  LOG(INFO) << "Opening file " << source;
  std::ifstream infile(source.c_str());
  CHECK((bool)infile) << "Could not open " << source;
  string line;
  while (getline(infile, line))
  {
    example ex(scache, line, this->layer_param_.molgrid_data_param());
    data->add(ex);
  }
  CHECK_GT(data->size(),0) << "No examples provided in source: " << source;

  data->setup(); //done adding

}

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::setBlobShape(const vector<Blob<Dtype>*>& top, 
    bool hasrmsd, bool hasaffinity) {
  int batch_size = batch_transform.size();
  //setup shape of layer
  top_shape.clear();
  top_shape.push_back(batch_size);
  top_shape.push_back(numReceptorTypes+numLigandTypes);
  top_shape.push_back(dim);
  top_shape.push_back(dim);
  top_shape.push_back(dim);

  example_size = (numReceptorTypes+numLigandTypes)*numgridpoints;

  // Reshape prefetch_data and top[0] according to the batch_size.
  top[0]->Reshape(top_shape);

  // Reshape label, affinity, rmsds
  vector<int> label_shape(1, batch_size); // [batch_size]

  top[1]->Reshape(label_shape);

  if (hasaffinity)
  {
    top[2]->Reshape(label_shape);
    if (hasrmsd)
      top[3]->Reshape(label_shape);
  else if(hasrmsd)
    top[2]->Reshape(label_shape);
  }

  if(ligpeturb) {
    vector<int> peturbshape(2);
    peturbshape[0] = batch_size;
    peturbshape[1] = 6; //trans+orient
    top.back()->Reshape(peturbshape);
  }

}

template <typename Dtype>
void SubcubeMolGridDataLayer<Dtype>::setBlobShape(const vector<Blob<Dtype>*>& top, 
    bool hasrmsd, bool hasaffinity) {
  //layer shape is TxNx... for RNN
  int batch_size = this->batch_transform.size();
  unsigned grids_per_dim = this->gmaker.grids_per_dim;
  unsigned n_timesteps = grids_per_dim * grids_per_dim * grids_per_dim;

  //setup shape of layer
  this->top_shape.clear();
  this->top_shape.push_back(n_timesteps);
  this->top_shape.push_back(batch_size);
  this->top_shape.push_back(this->numReceptorTypes+this->numLigandTypes);
  this->top_shape.push_back(this->gmaker.subgrid_dim_in_points);
  this->top_shape.push_back(this->gmaker.subgrid_dim_in_points);
  this->top_shape.push_back(this->gmaker.subgrid_dim_in_points);

  //the subcube recurrence manually indexes into the correct location; this
  //value is effectively a dummy value set for compatibility with parent class
  //functions
  this->example_size = 0;

  // Reshape prefetch_data and top[0] according to the batch_size.
  top[0]->Reshape(this->top_shape);

  // Reshape label, affinity, rmsds
  vector<int> label_shape;
  label_shape.push_back(n_timesteps);
  label_shape.push_back(batch_size);

  top[1]->Reshape(label_shape);
  if (hasaffinity)
  {
    top[2]->Reshape(label_shape);
    if (hasrmsd)
      top[3]->Reshape(label_shape);
  else if(hasrmsd)
    top[2]->Reshape(label_shape);
  }

  //RNN layer requires a TxN "sequence continuation" blob
  seqcont_shape.clear();
  seqcont_shape.push_back(n_timesteps);
  seqcont_shape.push_back(batch_size);
  int idx = this->ExactNumTopBlobs() - this->ligpeturb - 1;
  top[idx]->Reshape(seqcont_shape);

  if(this->ligpeturb) {
    vector<int> peturbshape(2);
    peturbshape[0] = batch_size;
    peturbshape[1] = 6; //trans+orient
    top.back()->Reshape(peturbshape);
  }
}

//read in structure input and atom type maps
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {

  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
  bool duplicate = param.duplicate_poses();

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

  if(!inmem)
  {

    const string& source = param.source();
    const string& source2 = param.source2();
    root_folder = sanitize_path(param.root_folder());
    root_folder2 = param.root_folder2();
    if(root_folder2.length() > 0) root_folder2 = sanitize_path(root_folder2);
    else root_folder2 = root_folder; //fall back on first

    CHECK_GT(source.length(), 0) << "No data source file provided";

    // Read source file(s) with labels and structures,
    // each line is label [group] [affinity] [rmsd] receptor_file ligand_file
    data = create_example_data(param);
    populate_data(root_folder, source, data, hasaffinity, hasrmsd, hasgroup);

    if(source2.length() > 0)
    {
      CHECK_GE(data_ratio, 0) << "Must provide non-negative ratio for two data sources";
      data2 = create_example_data(param);
      populate_data(root_folder2, source2, data2, hasaffinity, hasrmsd, hasgroup);
    }

    LOG(INFO) << "Total examples: " << data->size() + (data2 ? data2->size() : 0);

    // Check if we would need to randomly skip a few data points
    if (param.rand_skip())
    {
      unsigned int skip = caffe_rng_rand() %  param.rand_skip();

      LOG(INFO) << "Skipping first " << skip << " data points from each source.";

      example dummy;
      for(unsigned i = 0; i < skip; i++) {
        data->next(dummy);
      }
      if(data2)
      {
        for(unsigned i = 0; i < skip; i++) {
          data2->next(dummy);
        }
      }
    }
  }
  else //in memory always batch size of 1
  {
    batch_size = 1;
  }

  //initialize atom type maps
  string recmapfile = param.recmap();  //these are file names
  string ligmapfile = param.ligmap();

  //these are the actual contents
  string recmapstr = param.mem_recmap();
  string ligmapstr = param.mem_ligmap();

  //can specify maps programatically
  if(recmapstr.size())
    numReceptorTypes = gmaker.createMapFromString(recmapstr, rmap);
  else if(recmapfile.size() > 0)
    numReceptorTypes = gmaker.createAtomTypeMap(recmapfile, rmap);
  else
    numReceptorTypes = gmaker.createDefaultRecMap(rmap);


  if(ligmapstr.size())
    numLigandTypes = gmaker.createMapFromString(ligmapstr, lmap);
  else if (ligmapfile.size() > 0)
    numLigandTypes = gmaker.createAtomTypeMap(ligmapfile, lmap);
  else
    numLigandTypes = gmaker.createDefaultLigMap(lmap);

  CHECK_GT(batch_size, 0) << "Positive batch size required";
  //keep track of atoms and transformations for each example in batch
  batch_transform.resize(batch_size * param.maxgroupsize());

  //if specified, preload all gninatype information
  string reccache = param.recmolcache();
  string ligcache = param.ligmolcache();

  if(reccache.size() > 0) {
    load_cache(reccache, rmap, 0, recmolcache);
  }

  if(ligcache.size() > 0) {
    load_cache(ligcache, rmap, numReceptorTypes, ligmolcache);
  }

  //setup shape of layer
  top_shape.clear();
  unsigned number_examples = batch_size;
  if(duplicate) number_examples = batch_size*numposes;
  top_shape.push_back(number_examples);
  
  numchannels = numReceptorTypes+numLigandTypes;
  if(!duplicate && numposes > 1) numchannels = numReceptorTypes+numposes*numLigandTypes;
  top_shape.push_back(numchannels);
  top_shape.push_back(dim);
  top_shape.push_back(dim);
  top_shape.push_back(dim);

  example_size = (numchannels)*numgridpoints;

  // Reshape prefetch_data and top[0] according to the batch_size.
  top[0]->Reshape(top_shape);

  // Reshape label, affinity, rmsds
  vector<int> label_shape(1, number_examples); // [batch_size]

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

  if(param.affinity_reweight_stdcut() > 0) {
    unsigned indx = top.size()-1;
    if(ligpeturb) indx--; //this is getting cumbersome
    top[indx]->Reshape(label_shape);
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
typename BaseMolGridDataLayer<Dtype, GridMakerT>::quaternion BaseMolGridDataLayer<Dtype, GridMakerT>::axial_quaternion()
{
  unsigned rot = current_rotation;
  qt ret;
  //first rotate to a face
  switch(rot%6) {
    case 0:
      ret = quaternion(1,0,0,0); //identity
      break;
    case 1: //rotate z 90
      ret = quaternion(sqrt(0.5),0,0,sqrt(0.5));
      break;
    case 2: //rotate z 180
      ret = quaternion(0,0,0,1);
      break;
    case 3: //rotate z 270 (-90)
      ret = quaternion(sqrt(0.5),0,0,-sqrt(0.5));
      break;
    case 4: //rotate y 90
      ret = quaternion(sqrt(0.5),0,sqrt(0.5),0);
      break;
    case 5: //rotate y -90
      ret = quaternion(sqrt(0.5),0,-sqrt(0.5),0);
      break;
  }

  //now four rotations around x axis
  rot /= 6;
  switch(rot%4) {
    case 0:
      break; //identity
    case 1: //90 degrees
      ret *= quaternion(sqrt(0.5),sqrt(0.5),0,0);
      break;
    case 2: //180
      ret *= quaternion(0,1,0,0);
      break;
    case 3: //270
      ret *= quaternion(sqrt(0.5),-sqrt(0.5),0,0);
      break;
  }
  return ret;
}

//add atom information to minfo, return true if atom actually added
template <typename Dtype, class GridMakerT>
bool BaseMolGridDataLayer<Dtype, GridMakerT>::add_to_minfo(const string& file, const vector<int>& atommap, unsigned mapoffset, smt t, float x, float y, float z,  mol_info& minfo)
{
  int index = atommap[t];
  if(index >= 0)
  {
    float4 ainfo;
    ainfo.x = x;
    ainfo.y = y;
    ainfo.z  = z;
    if(fixedradius <= 0)
      ainfo.w = use_covalent_radius ? covalent_radius(t) : xs_radius(t);
    else
      ainfo.w = fixedradius;

    float3 gradient(0,0,0);
    minfo.atoms.push_back(ainfo);
    minfo.whichGrid.push_back(index+mapoffset);
    minfo.gradient.push_back(gradient);
  }
  else
  {
    static bool madewarning = false;
    if(!madewarning) {
      LOG(WARNING) << "WARNING: Unknown atom type " << t << " in " << file << ".  This atom will be discarded.  Future warnings will be suppressed\n";
      madewarning = true;
    }
    return 0;
  }
  return 1;
}

//two shared caches for the whole program
//so ugly, plz replace with macro
template<>
BaseMolGridDataLayer<float, GridMaker>::MolCache BaseMolGridDataLayer<float, GridMaker>::recmolcache = BaseMolGridDataLayer<float, GridMaker>::MolCache();
template<>
BaseMolGridDataLayer<double, GridMaker>::MolCache BaseMolGridDataLayer<double, GridMaker>::recmolcache =  BaseMolGridDataLayer<double, GridMaker>::MolCache();
template<>
BaseMolGridDataLayer<float, SubcubeGridMaker>::MolCache BaseMolGridDataLayer<float, SubcubeGridMaker>::recmolcache = BaseMolGridDataLayer<float, SubcubeGridMaker>::MolCache();
template<>
BaseMolGridDataLayer<double, SubcubeGridMaker>::MolCache BaseMolGridDataLayer<double, SubcubeGridMaker>::recmolcache =  BaseMolGridDataLayer<double, SubcubeGridMaker>::MolCache();

template<>
BaseMolGridDataLayer<float, GridMaker>::MolCache BaseMolGridDataLayer<float, GridMaker>::ligmolcache = BaseMolGridDataLayer<float, GridMaker>::MolCache();
template<>
BaseMolGridDataLayer<double, GridMaker>::MolCache BaseMolGridDataLayer<double, GridMaker>::ligmolcache =  BaseMolGridDataLayer<double, GridMaker>::MolCache();
template<>
BaseMolGridDataLayer<float, SubcubeGridMaker>::MolCache BaseMolGridDataLayer<float, SubcubeGridMaker>::ligmolcache = BaseMolGridDataLayer<float, SubcubeGridMaker>::MolCache();
template<>
BaseMolGridDataLayer<double, SubcubeGridMaker>::MolCache BaseMolGridDataLayer<double, SubcubeGridMaker>::ligmolcache =  BaseMolGridDataLayer<double, SubcubeGridMaker>::MolCache();

//load custom formatted cache file of all gninatypes into specified molcache using specified mapping and offset
template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::load_cache(const string& file, const vector<int>& atommap, unsigned mapoffset, BaseMolGridDataLayer<Dtype, GridMakerT>::MolCache& molcache)
{
  //file shoudl be organized
  //name size (1byte)
  //name (string)
  //number of atoms (4bytes)
  //atoms (3 floats and an int of the type)

  char buffer[257] = {0,};
  struct info {
    float x,y,z;
    int type;
  } atom;

  string fullpath = file;
  if(file.size() > 0 && file[0] != '/')
    fullpath = root_folder + file; //prepend dataroot if not absolute
  ifstream in(fullpath.c_str());
  CHECK(in) << "Could not read " << fullpath;

  LOG(INFO) << "Loading from " << fullpath << " with cache at size " << molcache.size() << "\n";
  while(in && in.peek() != EOF)
  {
    char sz = 0;
    int natoms = 0;
    in.read(&sz, sizeof(char));
    in.read(buffer,sizeof(char)*sz);
    buffer[(int)sz] = 0; //null terminate
    string fname(buffer);

    in.read((char*)&natoms, sizeof(int));

    if(molcache.count(fname)) {
      static int warncnt = 0;

      if(warncnt == 0) {
        LOG(WARNING) << "File " << fname << " duplicated in provided cache " << file << ".  Future warnings are supressed.";
        warncnt++;
      }
    }

    mol_info& minfo = molcache[fname];
    minfo.atoms.clear();
    minfo.whichGrid.clear();
    minfo.gradient.clear();
    int cnt = 0;
    vec center(0,0,0);

    for(unsigned i = 0; i < natoms; i++)
    {
      in.read((char*)&atom, sizeof(atom));
      smt t = (smt)atom.type;

      if(add_to_minfo(fname, atommap, mapoffset, t, atom.x, atom.y, atom.z, minfo)) {
        cnt++;
        center += vec(atom.x,atom.y,atom.z);
      }
    }

    if(cnt == 0) {
      LOG(WARNING) << "WARNING: No atoms in " << file <<"\n";
      continue;
    }

    center /= cnt;
    minfo.center = center;

    if(this->layer_param_.molgrid_data_param().fix_center_to_origin()) {
      minfo.center = vec(0,0,0);
    }
  }

  LOG(INFO) << "Done loading from " << fullpath << " with cache at size " << molcache.size() << std::endl;

}

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::set_mol_info(const string& file, const vector<int>& atommap,
    unsigned mapoffset, mol_info& minfo)
{
  //read mol info from file
  //OpenBabel is SLOW, especially for the receptor, so we cache the result
  //if this gets too annoying, can add support for spawning a thread for openbabel
  //but since this gets amortized across many hits to the same example, not a high priority
  //if clear is true, set, otherwise add
  using namespace OpenBabel;

  minfo.atoms.clear();
  minfo.whichGrid.clear();
  minfo.gradient.clear();
  
  int cnt = 0;
  vec center(0,0,0);

  //also, implemented a custom gninatypes files to precalc this info
  if(boost::algorithm::ends_with(file,".gninatypes"))
  {
    struct info {
      float x,y,z;
      int type;
    } atom;

    ifstream in(file.c_str());
    CHECK(in) << "Could not read " << file;

    while(in.read((char*)&atom, sizeof(atom)))
    {
      smt t = (smt)atom.type;

      if(add_to_minfo(file, atommap, mapoffset, t, atom.x, atom.y, atom.z, minfo)) {
        cnt++;
        center += vec(atom.x,atom.y,atom.z);
      }
    }
  }
  else if(!boost::algorithm::ends_with(file,"none")) //reserved word
  {
    //read mol from file and set mol info (atom coords and grid positions)
    //types are mapped using atommap values plus offset
    OpenBabel::OBConversion conv;
    OBMol mol;
    CHECK(conv.ReadFile(&mol, file)) << "Could not read " << file;

    if(this->layer_param_.molgrid_data_param().addh()) {
      mol.AddHydrogens();
    }

    minfo.atoms.reserve(mol.NumHvyAtoms());
    minfo.whichGrid.reserve(mol.NumHvyAtoms());
    minfo.gradient.reserve(mol.NumHvyAtoms());

    FOR_ATOMS_OF_MOL(a, mol)
    {
      smt t = obatom_to_smina_type(*a);
      if(add_to_minfo(file, atommap, mapoffset, t, a->x(), a->y(), a->z(), minfo)) {
        cnt++;
        center += vec(a->x(), a->y(), a->z());
      }
    }
  }

  if(cnt == 0) {
    std::cerr << "WARNING: No atoms in " << file <<"\n";
  }
  else {
    center /= cnt;
  }
  minfo.center = center;

  if(this->layer_param_.molgrid_data_param().fix_center_to_origin()) {
    minfo.center = vec(0,0,0);
  }

}

template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::set_grid_ex(Dtype *data, 
    const BaseMolGridDataLayer<Dtype, GridMakerT>::example& ex,
    const string& root_folder, BaseMolGridDataLayer<Dtype, GridMakerT>::mol_transform& transform, 
    int pose, output_transform& peturb, bool gpu)
{
  //set grid values for example
  //cache atom info
  //pose specifies which ligand pose to use (relevant if num_poses > 1)
  //if it is negative, use them all (but with distinct channels
  //data should be positioned at the start of the example
  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
  bool docache = param.cache_structs();
  bool doall = false;
  if(pose < 0) {
      doall = true;
      pose = 0;
  }

  CHECK_LT(pose, ex.ligands.size()) << "Incorrect pose index";
  const char* ligand = ex.ligands[pose];

  if(docache)
  {
    if(recmolcache.count(ex.receptor) == 0)
    {
      set_mol_info(root_folder+ex.receptor, rmap, 0, recmolcache[ex.receptor]);
    }
    if(ligmolcache.count(ligand) == 0)
    {
      set_mol_info(root_folder+ligand, lmap, numReceptorTypes, ligmolcache[ligand]);
    }

    if(doall) {
      //make sure every ligand is in the cache, then aggregate
      mol_info lig(ligmolcache[ligand]);
      for(unsigned p = 1, np = ex.ligands.size(); p < np; p++) {
        ligand = ex.ligands[p];
        if(ligmolcache.count(ligand) == 0)
        {
          set_mol_info(root_folder+ligand, lmap, numReceptorTypes, ligmolcache[ligand]);
        }
        lig.append(ligmolcache[ligand],numLigandTypes*p);
      }
      set_grid_minfo(data, recmolcache[ex.receptor], lig, transform, peturb, gpu);
    } else {
        set_grid_minfo(data, recmolcache[ex.receptor], ligmolcache[ligand], transform, peturb, gpu);
    }
  }
  else
  {
    mol_info rec;
    mol_info lig;
    set_mol_info(root_folder+ex.receptor, rmap, 0, rec);
    set_mol_info(root_folder+ligand, lmap, numReceptorTypes, lig);
    if(doall) {
        for(unsigned p = 1, np = ex.ligands.size(); p < np; p++) {
          mol_info tmplig;
          set_mol_info(root_folder+ligand, lmap, numReceptorTypes+numLigandTypes*p, tmplig);
          lig.append(tmplig);
        }
    }    
    set_grid_minfo(data, rec, lig, transform, peturb, gpu);
  }
}


template <typename Dtype, class GridMakerT>
void BaseMolGridDataLayer<Dtype, GridMakerT>::set_grid_minfo(Dtype *data, 
    const BaseMolGridDataLayer<Dtype, GridMakerT>::mol_info& recatoms,
    const BaseMolGridDataLayer<Dtype, GridMakerT>::mol_info& ligatoms, 
    BaseMolGridDataLayer<Dtype, GridMakerT>::mol_transform& transform,
    output_transform& peturb, bool gpu)
{
  //set grid values from mol info
  //first clear transform from the previous batch
  rng_t* rng = caffe_rng();
  transform = mol_transform();
  mol_transform ligtrans;

  //include receptor and ligand atoms
  transform.mol.append(recatoms);
  //set center to ligand center
  transform.mol.center = ligatoms.center;

  if(ligpeturb) {
    if(ligpeturb_rotate)
    {
      ligtrans.set_random_quaternion(rng);
    }
    else
    {
      ligtrans.Q = quaternion(1,0,0,0); //identity
    }
    ligtrans.add_random_displacement(rng, ligpeturb_translate);
    transform.mol.transform_and_append(ligatoms, ligtrans);

    //store the inverse transformation
    peturb.x = -ligtrans.center[0];
    peturb.y = -ligtrans.center[1];
    peturb.z = -ligtrans.center[2];

    qt qinv = conj(ligtrans.Q)/norm(ligtrans.Q); //not Cayley, not euclidean norm - already squared
    peturb.set_from_quaternion(qinv);

    //set the center to the translated value
    transform.mol.center = ligatoms.center + ligtrans.center;
  } else if(ignore_ligand) {
    //do nothing - ligand is only used to set center
  } else {
    transform.mol.append(ligatoms);
  }

  //figure out transformation
  transform.Q = quaternion(1,0,0,0);

  if(current_rotation == 0 && !randrotate)
    transform.Q = quaternion(1,0,0,0); //check real part to avoid mult

  if (randrotate)
  {
    transform.set_random_quaternion(rng);
  }

  transform.center[0] = transform.mol.center[0];
  transform.center[1] = transform.mol.center[1];
  transform.center[2] = transform.mol.center[2];
  if (randtranslate)
  {
    double radius = ligatoms.radius();
    //don't let ligand atoms translate out of sphere inscribed in box
    if(ignore_ligand) radius = 0;
    double maxtrans = max(dimension/2.0 - radius,0.0);
    updateTranslations(transform.add_random_displacement(rng, min(randtranslate,maxtrans)));
  }

  if(current_rotation > 0) {
    transform.Q *= axial_quaternion();
  }

  //TODO move this into gridmaker.setAtoms, have it just take the mol_transform as input - separate receptor transform as well
  gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);

  if(transform.mol.atoms.size() == 0) {
     std::cerr << "ERROR: No atoms in molecule.  I can't deal with this.\n";
     exit(-1); //presumably you never actually want this and it results in a cuda error
  } 
  if(jitter > 0) {
    //add small random displacement (in-place) to atoms
    for(unsigned i = 0, n = transform.mol.atoms.size(); i < n; i++) {
      float4& atom = transform.mol.atoms[i];
      float xdiff = jitter*(unit_sample(rng)*2.0-1.0);
      atom.x += xdiff;
      float ydiff = jitter*(unit_sample(rng)*2.0-1.0);
      atom.y += ydiff;
      float zdiff = jitter*(unit_sample(rng)*2.0-1.0);
      atom.z += zdiff;
    }
  }

  //compute grid from atom info arrays
  if(gpu)
  {
    unsigned natoms = transform.mol.atoms.size();
    allocateGPUMem(natoms);
    CUDA_CHECK(cudaMemcpy(gpu_gridatoms, &transform.mol.atoms[0], natoms*sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gpu_gridwhich, &transform.mol.whichGrid[0], natoms*sizeof(short), cudaMemcpyHostToDevice));

    gmaker.template setAtomsGPU<Dtype>(natoms, gpu_gridatoms, gpu_gridwhich, transform.Q, numchannels, data);
  }
  else
  {
    gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q.boost(), data, numchannels);
  }
}

template <typename Dtype>
void GroupedMolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data, 
    const typename BaseMolGridDataLayer<Dtype, GridMaker>::mol_info& recatoms,
    const typename BaseMolGridDataLayer<Dtype, GridMaker>::mol_info& ligatoms, 
    typename BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform& transform,
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
    unsigned transform_idx = example_idx % batch_size;
    transform.Q = this->batch_transform[transform_idx].Q;
    transform.center = translations[transform_idx];
    typename BaseMolGridDataLayer<Dtype, GridMaker>::mol_transform ligtrans;

    //include receptor and ligand atoms
    transform.mol.append(recatoms);
    //set center to ligand center
    transform.mol.center = ligatoms.center;

    if(this->ligpeturb) {
      ligtrans.center[0] = -this->perturbations[example_idx % batch_size].x;
      ligtrans.center[1] = -this->perturbations[example_idx % batch_size].y;
      ligtrans.center[2] = -this->perturbations[example_idx % batch_size].z;
      //FIXME: is this quaternion already normalized? if so this is wrong
      ligtrans.Q = conj(this->perturbations[example_idx % batch_size].get_quaternion());
      transform.mol.transform_and_append(ligatoms, ligtrans);

      //store the inverse transformation
      peturb = this->perturbations[example_idx % batch_size];

      //set the center to the translated value
      transform.mol.center = ligatoms.center + ligtrans.center;
    } else if(this->ignore_ligand) {
      //do nothing - ligand is only used to set center
    } else {
      transform.mol.append(ligatoms);
    }

    transform.center += transform.mol.center;

    //TODO move this into gridmaker.setAtoms, have it just take the mol_transform as input - separate receptor transform as well
    this->gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);

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

  if ((maxgroupsize-1) && !maxchunksize)
    maxchunksize = maxgroupsize;

  bool hasweights = (this->layer_param_.molgrid_data_param().affinity_reweight_stdcut() > 0);
  bool duplicate = this->layer_param_.molgrid_data_param().duplicate_poses();

  Dtype *top_data = NULL;
  if(gpu)
    top_data = top[0]->mutable_gpu_data();
  else
    top_data = top[0]->mutable_cpu_data();

  perturbations.clear();

  unsigned batch_size;
  if (subgrid_dim || (maxgroupsize-1)) {
    CHECK_EQ(numposes, 1); //TODO: combine multipose with groups?
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
    if (data2)
      dataswitch = batch_size*data_ratio/(data_ratio+1);

    for (int idx = 0; idx < batch_size + (maxchunksize-1)*batch_size; ++idx)
    {
      int batch_idx = idx % batch_size;
      example ex;
      string *root;
      if (batch_idx < dataswitch)
      {
        data->next(ex);
        root = &root_folder;
      }
      else
      {
        data2->next(ex);
        root = &root_folder2;
      }

      updateLabels(ex.label, ex.affinity, ex.rmsd, ex.affinity_weight);
      int step = idx / batch_size;
      int offset = ((batch_size * step) + batch_idx) * example_size;

      //if label == -1 then this is a padding example for grouped data; just
      //memset data to 0
      if (ex.label == -1) {
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
          set_grid_ex(top_data+offset, ex, *root, batch_transform[batch_idx], numposes > 1 ? -1 : 0, peturb, gpu);
          perturbations.push_back(peturb);
        }
        else {
      	  for(unsigned p = 0; p < numposes; p++) {
            int p_offset = batch_idx*(example_size*numposes)+example_size*p;
            set_grid_ex(top_data+p_offset, ex, *root, batch_transform[batch_idx], p, peturb, gpu);
            perturbations.push_back(peturb);
            //NOTE: num_rotations not actually implemented!
          }
        }
        //NOTE: batch_transform contains transformation of last pose only - don't use unless numposes == 1
      }

    }

  }

  copyToBlobs(top, hasaffinity, hasrmsd, hasweights, gpu);
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
        mol_transform& transform = batch_transform[item_id];
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
    mol_transform& transform = batch_transform[item_id];
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
void GroupedMolGridDataLayer<Dtype>::setBlobShape(const vector<Blob<Dtype>*>& top, 
    bool hasrmsd, bool hasaffinity) {
  int batch_size = this->batch_transform.size() / maxgroupsize;
  this->top_shape.clear();
  this->top_shape.push_back(maxchunksize);
  this->top_shape.push_back(batch_size);
  this->top_shape.push_back(this->numReceptorTypes+this->numLigandTypes);
  this->top_shape.push_back(this->dim);
  this->top_shape.push_back(this->dim);
  this->top_shape.push_back(this->dim);

  this->example_size = (this->numReceptorTypes+this->numLigandTypes)*this->numgridpoints;
  // Reshape prefetch_data and top[0] according to the batch_size.
  top[0]->Reshape(this->top_shape);

  // Reshape label, affinity, rmsds
  vector<int> label_shape;
  label_shape.push_back(maxchunksize);
  label_shape.push_back(batch_size);

  top[1]->Reshape(label_shape);
  if (hasaffinity)
  {
    top[2]->Reshape(label_shape);
    if (hasrmsd)
      top[3]->Reshape(label_shape);
  else if(hasrmsd)
    top[2]->Reshape(label_shape);
  }

  //RNN layer requires a TxN "sequence continuation" blob
  seqcont_shape.clear();
  seqcont_shape.push_back(maxchunksize);
  seqcont_shape.push_back(batch_size);
  int idx = this->ExactNumTopBlobs() - this->ligpeturb - 1;
  top[idx]->Reshape(seqcont_shape);

  if(this->ligpeturb) {
    vector<int> peturbshape(2);
    peturbshape[0] = batch_size;
    peturbshape[1] = 6; //trans+orient
    top.back()->Reshape(peturbshape);
  }
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
