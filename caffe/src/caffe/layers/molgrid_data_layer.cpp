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
RealMolGridDataLayer<Dtype, GridMakerT>::~RealMolGridDataLayer<Dtype, GridMakerT>() {
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
RealMolGridDataLayer<Dtype, GridMakerT>::example::example(RealMolGridDataLayer<Dtype, GridMakerT>::string_cache& cache, string line, bool hasaffinity, bool hasrmsd)
  : label(0), affinity(0.0), rmsd(0.0)
{
  stringstream stream(line);
  string tmp;
  //first the label
  stream >> label;
  if(hasaffinity)
   stream >> affinity;
  if(hasrmsd)
   stream >> rmsd;
  //receptor
  stream >> tmp;
  CHECK(tmp.length() > 0) << "Empty receptor, missing affinity/rmsd?";
  receptor = cache.get(tmp);
  //ligand
  tmp.clear();
  stream >> tmp;
  CHECK(tmp.length() > 0) << "Empty ligand, missing affinity/rmsd?";

  ligand = cache.get(tmp);
}


//for in-memory inputs, set the desired label (for gradient computation)
//note that zero affinity/rmsd means to ignore these
template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::setLabels(Dtype pose, Dtype affinity, Dtype rmsd)
{
  labels.clear();
  affinities.clear();
  rmsds.clear();
  labels.push_back(pose);
  affinities.push_back(affinity);
  rmsds.push_back(rmsd);
}


//the following really shouldn't be recalculated each evaluation (not including gradients)
template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::getReceptorAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
      atoms.push_back(mol.atoms[i]);
}

template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::getLigandAtoms(int batch_idx, vector<float4>& atoms)
{
  atoms.resize(0);
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
      atoms.push_back(mol.atoms[i]);
}

template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::getReceptorChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] < numReceptorTypes)
      whichGrid.push_back(mol.whichGrid[i]);
}

template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::getLigandChannels(int batch_idx, vector<short>& whichGrid)
{
  whichGrid.resize(0);
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
      whichGrid.push_back(mol.whichGrid[i]);
}

template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::getReceptorGradient(int batch_idx, vector<float3>& gradient)
{
  gradient.resize(0);
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
void RealMolGridDataLayer<Dtype, GridMakerT>::getReceptorTransformationGradient(int batch_idx, vec& force, vec& torque)
{
  force = vec(0,0,0);
  torque = vec(0,0,0);

  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  mol_info& mol = batch_transform[batch_idx].mol;
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
void RealMolGridDataLayer<Dtype, GridMakerT>::getMappedReceptorGradient(int batch_idx, unordered_map<string, float3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
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
void RealMolGridDataLayer<Dtype, GridMakerT>::getLigandGradient(int batch_idx, vector<float3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
  gradient.resize(0);
  mol_info& mol = batch_transform[batch_idx].mol;
  for (unsigned i = 0, n = mol.atoms.size(); i < n; ++i)
    if (mol.whichGrid[i] >= numReceptorTypes)
    {
      gradient.push_back(mol.gradient[i]);
    }
}

template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::getMappedLigandGradient(int batch_idx, unordered_map<string, float3>& gradient)
{
  CHECK(compute_atom_gradients) << "Gradients requested but not computed";
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
void RealMolGridDataLayer<Dtype, GridMakerT>::remove_missing_and_setup(vector<typename RealMolGridDataLayer<Dtype, GridMakerT>::balanced_example_provider>& examples)
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
      LOG(INFO) << "Dropping receptor " << tmp.receptor << " with no decoys.";
    }
  }

  swap(examples,tmp);
}

//specialized version for balanced data that remove receptors without any actives or decoys
//annoyingly, have to specialize Dtype and GridMakerT - workaround?
template <>
template <>
void RealMolGridDataLayer<float, GridMaker>::receptor_stratified_example_provider<typename RealMolGridDataLayer<float, GridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template<>
template<>
void RealMolGridDataLayer<double, GridMaker>::receptor_stratified_example_provider<typename RealMolGridDataLayer<double, GridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template <>
template <>
void RealMolGridDataLayer<float, RNNGridMaker>::receptor_stratified_example_provider<typename RealMolGridDataLayer<float, RNNGridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template<>
template<>
void RealMolGridDataLayer<double, RNNGridMaker>::receptor_stratified_example_provider<typename RealMolGridDataLayer<double, RNNGridMaker>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}


//ensure gpu memory is of sufficient size
template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::allocateGPUMem(unsigned sz)
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
typename RealMolGridDataLayer<Dtype, GridMakerT>::example_provider* RealMolGridDataLayer<Dtype, GridMakerT>::create_example_data(const MolGridDataParameter& parm)
{
  bool balanced  = parm.balanced();
  bool strat_receptor  = parm.stratify_receptor();
  bool strat_aff = parm.stratify_affinity_max() != parm.stratify_affinity_min();

  //strat_aff > strat_receptor > balanced
  if(strat_aff)
  {
    if(strat_receptor)
    {
      if(balanced) // sample 2 from each receptor
      {
        return new affinity_stratified_example_provider<receptor_stratified_example_provider<balanced_example_provider, 2> >(parm);
      }
      else //sample 1 from each receptor
      {
        return new affinity_stratified_example_provider<receptor_stratified_example_provider<uniform_example_provider, 1> >(parm);
      }
    }
    else
    {
      if(balanced)
      {
        return new affinity_stratified_example_provider<balanced_example_provider>(parm);
      }
      else //sample 1 from each receptor
      {
        return new affinity_stratified_example_provider<uniform_example_provider>(parm);
      }
    }
  }
  else if(strat_receptor)
  {
    if(balanced) // sample 2 from each receptor
    {
      return new receptor_stratified_example_provider<balanced_example_provider, 2>(parm);
    }
    else //sample 1 from each receptor
    {
      return new receptor_stratified_example_provider<uniform_example_provider, 1>(parm);
    }
  }
  else if(balanced)
  {
    return new balanced_example_provider(parm);
  }
  else
  {
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
void RealMolGridDataLayer<Dtype, GridMakerT>::populate_data(const string& root_folder, const string& source,
    RealMolGridDataLayer<Dtype, GridMakerT>::example_provider* data, bool hasaffinity, bool hasrmsd)
{
  LOG(INFO) << "Opening file " << source;
  std::ifstream infile(source.c_str());
  CHECK((bool)infile) << "Could not open " << source;
  string line;
  while (getline(infile, line))
  {
    example ex(scache, line, hasaffinity, hasrmsd);
    data->add(ex);
  }
  CHECK_GT(data->size(),0) << "No examples provided in source: " << source;

  data->setup(); //done adding

}

template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::setBlobShape(const vector<Blob<Dtype>*>& top, 
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
  //LSTM layer requires a "sequence continuation" blob
  vector<int> seqcont_shape((dimension/subgrid_dim) * (dimension/subgrid_dim) * 
      (dimension / subgrid_dim), batch_size);

  top[1]->Reshape(label_shape);

  if (hasaffinity)
  {
    top[2]->Reshape(label_shape);
    if (hasrmsd)
    {
      top[3]->Reshape(label_shape);
      if (subgrid_dim) {
        top[4]->Reshape(seqcont_shape);
      }
    }
    else if(subgrid_dim) {
      top[3]->Reshape(seqcont_shape);
    }
  }
  else if(hasrmsd)
  {
    top[2]->Reshape(label_shape);
    if (subgrid_dim) {
      top[3]->Reshape(seqcont_shape);
    }
  }
  else if (subgrid_dim) {
    top[2]->Reshape(seqcont_shape);
  }

  if(ligpeturb) {
    vector<int> peturbshape(2);
    peturbshape[0] = batch_size;
    peturbshape[1] = 6; //trans+orient
    top.back()->Reshape(peturbshape);
  }

}

//read in structure input and atom type maps
template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {

  const MolGridDataParameter& param = this->layer_param_.molgrid_data_param();
  root_folder = param.root_folder();
  num_rotations = param.rotate();
  inmem = param.inmemory();
  dimension = param.dimension();
  resolution = param.resolution();
  subgrid_dim = param.subgrid_dim();
  binary = param.binary_occupancy();
  bool spherize = param.spherical_mask();
  randtranslate = param.random_translate();
  randrotate = param.random_rotation();
  ligpeturb = param.peturb_ligand();
  ligpeturb_translate = param.peturb_ligand_translate();
  ligpeturb_rotate = param.peturb_ligand_rotate();
  radiusmultiple = param.radius_multiple();
  fixedradius = param.fixed_radius();
  bool hasaffinity = param.has_affinity();
  bool hasrmsd = param.has_rmsd();
  data_ratio = param.source_ratio();
  root_folder2 = param.root_folder2();

  if(binary) radiusmultiple = 1.0;

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
    // each line is label [affinity] [rmsd] receptor_file ligand_file
    data = create_example_data(param);
    populate_data(root_folder, source, data, hasaffinity, hasrmsd);

    if(source2.length() > 0)
    {
      CHECK_GE(data_ratio, 0) << "Must provide non-negative ratio for two data sources";
      data2 = create_example_data(param);
      populate_data(root_folder2, source2, data2, hasaffinity, hasrmsd);
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

  CHECK_GT(batch_size, 0) << "Positive batch size required";
  //keep track of atoms and transformations for each example in batch
  batch_transform.resize(batch_size);

  //initialize atom type maps
  string recmap = param.recmap();
  string ligmap = param.ligmap();

  if (recmap.size() == 0)
    numReceptorTypes = GridMaker::createDefaultRecMap(rmap);
  else
    numReceptorTypes = GridMaker::createAtomTypeMap(recmap, rmap);

  if (ligmap.size() == 0)
    numLigandTypes = GridMaker::createDefaultLigMap(lmap);
  else
    numLigandTypes = GridMaker::createAtomTypeMap(ligmap, lmap);

  setBlobShape(top, hasrmsd, hasaffinity);
}

//return quaternion representing one of 24 distinct axial rotations
template <typename Dtype, class GridMakerT>
typename RealMolGridDataLayer<Dtype, GridMakerT>::quaternion RealMolGridDataLayer<Dtype, GridMakerT>::axial_quaternion()
{
  using namespace boost::math;
  unsigned rot = current_rotation;
  quaternion ret;
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


template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::set_mol_info(const string& file, const vector<int>& atommap,
    unsigned mapoffset, mol_info& minfo)
{
  //read mol info from file
  //OpenBabel is SLOW, especially for the receptor, so we cache the result
  //if this gets too annoying, can add support for spawning a thread for openbabel
  //but since this gets amortized across many hits to the same example, not a high priority
  using namespace OpenBabel;

  vec center(0,0,0);
  minfo.atoms.clear();
  minfo.whichGrid.clear();
  minfo.gradient.clear();

  //also, implemented a custom gninatypes files to precalc this info
  if(boost::algorithm::ends_with(file,".gninatypes"))
  {
    struct info {
      float x,y,z;
      int type;
    } atom;

    ifstream in(file.c_str());
    CHECK(in) << "Could not read " << file;

    int cnt = 0;
    while(in.read((char*)&atom, sizeof(atom)))
    {
      smt t = (smt)atom.type;
      int index = atommap[t];

      if(index >= 0)
      {
        cnt++;
        float4 ainfo;
        ainfo.x = atom.x;
        ainfo.y = atom.y;
        ainfo.z = atom.z;
        if(fixedradius <= 0)
        	ainfo.w = xs_radius(t);
        else
        	ainfo.w = fixedradius;
        float3 gradient(0,0,0);

        minfo.atoms.push_back(ainfo);
        minfo.whichGrid.push_back(index+mapoffset);
        minfo.gradient.push_back(gradient);
        center += vec(atom.x,atom.y,atom.z);
      }
      else if(t > 1) //silence on hydrogens
      {
       std::cerr << "WARNING: Unknown atom type " << t << " in " << file << ".  This atom will be discarded\n";
      }
    }
    center /= cnt;
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

    int cnt = 0;
    FOR_ATOMS_OF_MOL(a, mol)
    {
      smt t = obatom_to_smina_type(*a);
      int index = atommap[t];

      if(index >= 0)
      {
        cnt++;
        float4 ainfo;
        ainfo.x = a->x();
        ainfo.y = a->y();
        ainfo.z  = a->z();
        if(fixedradius <= 0)
        	ainfo.w = xs_radius(t);
        else
        	ainfo.w = fixedradius;
        float3 gradient(0,0,0);

        minfo.atoms.push_back(ainfo);
        minfo.whichGrid.push_back(index+mapoffset);
        minfo.gradient.push_back(gradient);
        center += vec(a->x(),a->y(),a->z());
      }
      else
      {
        std::cerr << "WARNING: Unknown atom type in " << file << ".  This atom will be discarded\n";
      }	
    }
    center /= cnt;
  }

  minfo.center = center;

}

template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::set_grid_ex(Dtype *data, const RealMolGridDataLayer<Dtype, GridMakerT>::example& ex,
    const string& root_folder, RealMolGridDataLayer<Dtype, GridMakerT>::mol_transform& transform, output_transform& peturb, bool gpu, unsigned batch_size, unsigned batch_idx)
{
  //set grid values for example
  //cache atom info
  bool docache = this->layer_param_.molgrid_data_param().cache_structs();

  if(docache)
  {
    if(molcache.count(ex.receptor) == 0)
    {
      set_mol_info(root_folder+ex.receptor, rmap, 0, molcache[ex.receptor]);
    }
    if(molcache.count(ex.ligand) == 0)
    {
      set_mol_info(root_folder+ex.ligand, lmap, numReceptorTypes, molcache[ex.ligand]);
    }

    set_grid_minfo(data, molcache[ex.receptor], molcache[ex.ligand], transform, peturb, gpu, 
                   batch_idx);
  }
  else
  {
    mol_info rec;
    mol_info lig;
    set_mol_info(root_folder+ex.receptor, rmap, 0, rec);
    set_mol_info(root_folder+ex.ligand, lmap, numReceptorTypes, lig);
    set_grid_minfo(data, rec, lig, transform, peturb, gpu, batch_idx);
  }
}


template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::set_grid_minfo(Dtype *data, const RealMolGridDataLayer<Dtype, GridMakerT>::mol_info& recatoms,
  const RealMolGridDataLayer<Dtype, GridMakerT>::mol_info& ligatoms, RealMolGridDataLayer<Dtype, GridMakerT>::mol_transform& transform,
  output_transform& peturb, bool gpu, unsigned batch_size, unsigned batch_idx)
{
  //set grid values from mol info
  //first clear transform from the previous batch
  rng_t* rng = caffe_rng();
  transform = mol_transform();
  mol_transform ligtrans;

  //include receptor and ligand atoms
  transform.mol.append(recatoms);

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

    quaternion qinv = conj(ligtrans.Q)/norm(ligtrans.Q); //not Cayley, not euclidean norm - already squared
    peturb.set_from_quaternion(qinv);

  } else {
    transform.mol.append(ligatoms);
  }

  //set center to ligand center
  transform.mol.center = ligatoms.center;

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
    double maxtrans = max(dimension/2.0 - radius,0.0);
    transform.add_random_displacement(rng, min(randtranslate,maxtrans));
  }

  if(current_rotation > 0) {
    transform.Q *= axial_quaternion();
  }

  //TODO move this into gridmaker.setAtoms, have it just take the mol_transform as input
  gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);
 
  if(transform.mol.atoms.size() == 0) {
     std::cerr << "ERROR: No atoms in molecule.  I can't deal with this.\n";
     exit(-1); //presumably you never actually want this and it results in a cuda error
  } 
  //compute grid from atom info arrays
  if(gpu)
  {
    unsigned natoms = transform.mol.atoms.size();
    allocateGPUMem(natoms);
    CUDA_CHECK(cudaMemcpy(gpu_gridatoms, &transform.mol.atoms[0], natoms*sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gpu_gridwhich, &transform.mol.whichGrid[0], natoms*sizeof(short), cudaMemcpyHostToDevice));

    gmaker.template setAtomsGPU<Dtype>(natoms, gpu_gridatoms, gpu_gridwhich, transform.Q, numReceptorTypes+numLigandTypes, data, batch_idx, batch_size);
  }
  else
  {
    if (subgrid_dim) {
      unsigned ncubes = (dimension / subgrid_dim) * (dimension / subgrid_dim) * 
        (dimension / subgrid_dim);
      unsigned subgrid_dim_in_points = dim / (dimension / subgrid_dim);
      RNNGrids grids(data, boost::extents[ncubes][batch_size][numReceptorTypes+numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);
      gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids, 
                         batch_idx);
    }
    else {
      Grids grids(data, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);
      gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids, 
                         batch_idx);
    }
  }
}


//return a string representation of the atom type(s) represented by index
//in map - this isn't particularly efficient, but is only for debug purposes
template <typename Dtype, class GridMakerT>
string RealMolGridDataLayer<Dtype, GridMakerT>::getIndexName(const vector<int>& map, unsigned index) const
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
void RealMolGridDataLayer<Dtype, GridMakerT>::outputDXGrid(std::ostream& out, Grids& grid, unsigned g, double scale) const
{
  unsigned n = dim;
  if (subgrid_dim)
    n = dim / (dimension / subgrid_dim);
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
//if doing subcubes, output a separate file for each subcube (for now) to
//confirm that they look reasonable
template<typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::dumpDiffDX(const std::string& prefix,
		Blob<Dtype>* top, double scale) const
{
  Dtype* diff = top->mutable_cpu_diff();
	Grids grids(diff,
			boost::extents[numReceptorTypes + numLigandTypes][dim][dim][dim]);
    CHECK_GT(mem_lig.atoms.size(),0) << "DX dump only works with in-memory ligand";
    CHECK_EQ(randrotate, false) << "DX dump requires no rotation";
  unsigned ncubes;
  unsigned cubes_per_dim;
  unsigned subgrid_dim_in_points;
  unsigned instance_size;
	for (unsigned a = 0, na = numReceptorTypes; a < na; a++) {
		string name = getIndexName(rmap, a);
    if (subgrid_dim) {
      cubes_per_dim = dimension / subgrid_dim;
      ncubes = cubes_per_dim * cubes_per_dim * cubes_per_dim;
      subgrid_dim_in_points = dim / cubes_per_dim;
      instance_size = (numReceptorTypes + numLigandTypes) * 
        subgrid_dim_in_points * subgrid_dim_in_points * subgrid_dim_in_points;
      for (size_t i=0; i<ncubes; ++i) {
		    string fname = prefix + "_cube" + std::to_string(i) + "_rec_" + name + ".dx";
		    ofstream out(fname.c_str());
        unsigned offset = i * instance_size; 
        Grids subgrids(diff+offset, boost::extents[numReceptorTypes+numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);
		    outputDXGrid(out, subgrids, a, scale);
      }
    }
    else {
		  string fname = prefix + "_rec_" + name + ".dx";
		  ofstream out(fname.c_str());
		  outputDXGrid(out, grids, a, scale);
    }
	}
	for (unsigned a = 0, na = numLigandTypes; a < na; a++) {
			string name = getIndexName(lmap, a);
      if (subgrid_dim) {
        for (size_t i=0; i<ncubes; ++i) {
		      string fname = prefix + "_cube" + std::to_string(i) + "_lig_" + name + ".dx";
		      ofstream out(fname.c_str());
          unsigned offset = i * instance_size; 
          Grids subgrids(diff+offset, boost::extents[numReceptorTypes+numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);
		      outputDXGrid(out, subgrids, a, scale);
        }
      }
      else {
			  string fname = prefix + "_lig_" + name + ".dx";
			  ofstream out(fname.c_str());
			  outputDXGrid(out, grids, numReceptorTypes+a, scale);
      }
	}

}


template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top)
{
  forward(bottom, top, false);
}


template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu)
{
  bool hasaffinity = this->layer_param_.molgrid_data_param().has_affinity();
  bool hasrmsd = this->layer_param_.molgrid_data_param().has_rmsd();

  Dtype *top_data = NULL;
  if(gpu)
    top_data = top[0]->mutable_gpu_data();
  else
    top_data = top[0]->mutable_cpu_data();

  perturbations.clear();
  unsigned batch_size = top_shape[0];
  output_transform peturb;

  //if in memory must be set programmatically
  if(inmem)
  {
    CHECK_GT(mem_rec.atoms.size(),0) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(mem_lig.atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(top_data, mem_rec, mem_lig, batch_transform[0], peturb, gpu, 
        batch_size); //TODO how do we know what batch position?
    perturbations.push_back(peturb);

    if (num_rotations > 0) {
      current_rotation = (current_rotation+1)%num_rotations;
    }

    CHECK_GT(labels.size(),0) << "Did not set labels in memory based molgrid";

  }
  else
  {
    //clear batch labels
    labels.clear();
    affinities.clear();
    rmsds.clear();

    //percent of batch from first data source
    unsigned dataswitch = batch_size;
    if (data2)
      dataswitch = batch_size*data_ratio/(data_ratio+1);

    if (subgrid_dim) {
      unsigned ncubes = (dimension / subgrid_dim) * (dimension/ subgrid_dim) * 
        (dimension / subgrid_dim);
      unsigned n_examples = ncubes * batch_size;
      for (size_t i=0; i<n_examples; ++i) {
        if (i < batch_size)
          seqcont.push_back(0);
        else
          seqcont.push_back(1);
      }
    }
    for (int batch_idx = 0; batch_idx < batch_size; ++batch_idx)
    {
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

      labels.push_back(ex.label);
      affinities.push_back(ex.affinity);
      rmsds.push_back(ex.rmsd);

      int offset = batch_idx*example_size;
      if (subgrid_dim) 
        offset = 0;
      set_grid_ex(top_data+offset, ex, *root, batch_transform[batch_idx], peturb, gpu, 
          batch_size, batch_idx);
      perturbations.push_back(peturb);
      //NOTE: num_rotations not actually implemented!
    }

  }

  if(gpu) {
    caffe_copy(labels.size(), &labels[0], top[1]->mutable_gpu_data());
    if(hasaffinity) {
      caffe_copy(affinities.size(), &affinities[0], top[2]->mutable_gpu_data());
      if(hasrmsd) {
        caffe_copy(rmsds.size(), &rmsds[0], top[3]->mutable_gpu_data());
        if (subgrid_dim) {
          caffe_copy(seqcont.size(), &seqcont[0], top[4]->mutable_gpu_data());
        }
      }
      else if (subgrid_dim) {
        caffe_copy(seqcont.size(), &seqcont[0], top[3]->mutable_gpu_data());
      }
    } 
    else if(hasrmsd) {
      caffe_copy(rmsds.size(), &rmsds[0], top[2]->mutable_gpu_data());
        if (subgrid_dim) {
          caffe_copy(seqcont.size(), &seqcont[0], top[3]->mutable_gpu_data());
        }
    }
    else if (subgrid_dim) {
      caffe_copy(seqcont.size(), &seqcont[0], top[2]->mutable_gpu_data());
    }

    if(ligpeturb) {
      //trusting struct layout is normal
      caffe_copy(perturbations.size()*6, (Dtype*)&perturbations[0], top.back()->mutable_gpu_data());
    }

  }
  else {
    caffe_copy(labels.size(), &labels[0], top[1]->mutable_cpu_data());
    if(hasaffinity) {
      caffe_copy(affinities.size(), &affinities[0], top[2]->mutable_cpu_data());
      if(hasrmsd) {
        caffe_copy(rmsds.size(), &rmsds[0], top[3]->mutable_cpu_data());
        if (subgrid_dim) {
          caffe_copy(seqcont.size(), &seqcont[0], top[4]->mutable_cpu_data());
        }
      }
      else if (subgrid_dim) {
        caffe_copy(seqcont.size(), &seqcont[0], top[3]->mutable_cpu_data());
      }
    } 
    else if(hasrmsd) {
      caffe_copy(rmsds.size(), &rmsds[0], top[2]->mutable_cpu_data());
      if (subgrid_dim) {
        caffe_copy(seqcont.size(), &seqcont[0], top[3]->mutable_cpu_data());
      }
    }
    else if (subgrid_dim) {
      caffe_copy(seqcont.size(), &seqcont[0], top[2]->mutable_cpu_data());
    }
    if(ligpeturb) {
      //trusting struct layout is normal
      caffe_copy(perturbations.size()*6, (Dtype*)&perturbations[0], top.back()->mutable_cpu_data());
    }
  }

}


template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  backward(top, bottom, false);
}


/* backpropagates gradients onto atoms (note there is not actual bottom)
 * Only performed when compute_atom_gradients is true.
 */
template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom, bool gpu)
{
  //propagate gradient grid onto atom positions
  if(compute_atom_gradients || true) {
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

        if (subgrid_dim) {
          unsigned ncubes = (dimension / subgrid_dim) * (dimension / subgrid_dim) * 
            (dimension / subgrid_dim);
          unsigned subgrid_dim_in_points = dim / (dimension / subgrid_dim);
          Grids grids(diff, boost::extents[ncubes][batch_size][numReceptorTypes+numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);
          gmaker.setAtomGradientsCPU(transform.mol.atoms, transform.mol.whichGrid, 
                  transform.Q, grids, transform.mol.gradient, item_id);
        }
        else {
          Grids grids(diff+offset, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);
	      // boost::timer::cpu_timer time;
          gmaker.setAtomGradientsCPU(transform.mol.atoms, transform.mol.whichGrid, 
                  transform.Q, grids, transform.mol.gradient, item_id);
	      // std::cout << "CPU grid time " << time.elapsed().wall/1000000000.0 << "\n";
        }
      }
    }
  }
}

template <typename Dtype, class GridMakerT>
void RealMolGridDataLayer<Dtype, GridMakerT>::Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{

  /*
  float top_sum = 0.0;

  std::cout << "MOLGRID TOP:";
  for(int i = 0; i < top[0]->count(); i++)
  {
          top_sum += top[0]->cpu_diff()[i];
  }
  std::cout << "MOLGRID TOP: " << top_sum << '\n';
  */

  Dtype *diff = top[0]->mutable_cpu_diff(); //TODO: implement gpu

  //propagate gradient grid onto atom positions
  unsigned batch_size = top_shape[0];
  for (int item_id = 0; item_id < batch_size; ++item_id) {

    int offset = item_id*example_size;
    mol_transform& transform = batch_transform[item_id];
    gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);

    if (subgrid_dim) {
      unsigned ncubes = (dimension / subgrid_dim) * (dimension / subgrid_dim) * 
        (dimension / subgrid_dim);
      unsigned subgrid_dim_in_points = dim / (dimension / subgrid_dim);
      Grids grids(diff, boost::extents[ncubes][batch_size][numReceptorTypes+numLigandTypes][subgrid_dim_in_points][subgrid_dim_in_points][subgrid_dim_in_points]);

      gmaker.setAtomRelevanceCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids,
          transform.mol.gradient, item_id);
    }
    else {
      Grids grids(diff+offset, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);

      gmaker.setAtomRelevanceCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids,
          transform.mol.gradient, item_id);
    }
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
  if (mgrid_param.subgrid_dim()) 
    return shared_ptr<Layer<Dtype> >(new RealMolGridDataLayer<Dtype, RNNGridMaker>(param));
  else
    return shared_ptr<Layer<Dtype> >(new RealMolGridDataLayer<Dtype, GridMaker>(param));
}

INSTANTIATE_CLASS(MolGridDataLayer);
REGISTER_LAYER_CREATOR(MolGridData, GetMolGridDataLayer);

}  // namespace caffe
