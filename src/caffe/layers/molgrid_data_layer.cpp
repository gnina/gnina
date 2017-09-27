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

template <typename Dtype>
MolGridDataLayer<Dtype>::~MolGridDataLayer<Dtype>() {
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

template <typename Dtype>
MolGridDataLayer<Dtype>::example::example(MolGridDataLayer<Dtype>::string_cache& cache, string line, bool hasaffinity, bool hasrmsd)
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


//modify examples to remove any without both actives an inactives
//factored this into its own function due to the need to fully specialize setup below
template<typename Dtype>
void MolGridDataLayer<Dtype>::remove_missing_and_setup(vector<typename MolGridDataLayer<Dtype>::balanced_example_provider>& examples)
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
//annoyingly, have to specialize Dtype
template<>
template<>
void MolGridDataLayer<float>::receptor_stratified_example_provider<typename MolGridDataLayer<float>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}

template<>
template<>
void MolGridDataLayer<double>::receptor_stratified_example_provider<typename MolGridDataLayer<double>::balanced_example_provider, 2>::setup()
{
  currenti = 0; currentk = 0;
  remove_missing_and_setup(examples);
  //also shuffle receptors
  if(randomize) shuffle(examples.begin(), examples.end(), caffe::caffe_rng());
}



//ensure gpu memory is of sufficient size
template <typename Dtype>
void MolGridDataLayer<Dtype>::allocateGPUMem(unsigned sz)
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
template <typename Dtype>
typename MolGridDataLayer<Dtype>::example_provider* MolGridDataLayer<Dtype>::create_example_data(const MolGridDataParameter& parm)
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
template <typename Dtype>
void MolGridDataLayer<Dtype>::populate_data(const string& root_folder, const string& source,
    MolGridDataLayer<Dtype>::example_provider* data, bool hasaffinity, bool hasrmsd)
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

//read in structure input and atom type maps
template <typename Dtype>
void MolGridDataLayer<Dtype>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
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
  radiusmultiple = param.radius_multiple();
  fixedradius = param.fixed_radius();
  bool hasaffinity = param.has_affinity();
  bool hasrmsd = param.has_rmsd();
  data_ratio = param.source_ratio();
  root_folder2 = param.root_folder2();

  if(binary) radiusmultiple = 1.0;

  gmaker.initialize(resolution, dimension, radiusmultiple, binary);

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
    {
      top[3]->Reshape(label_shape);
    }
  }
  else if(hasrmsd)
  {
    top[2]->Reshape(label_shape);
  }
}

//return quaternion representing one of 24 distinct axial rotations
template <typename Dtype>
typename MolGridDataLayer<Dtype>::quaternion MolGridDataLayer<Dtype>::axial_quaternion()
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


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_mol_info(const string& file, const vector<int>& atommap,
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

template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_ex(Dtype *data, const MolGridDataLayer<Dtype>::example& ex,
    const string& root_folder, MolGridDataLayer<Dtype>::mol_transform& transform, bool gpu)
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

    set_grid_minfo(data, molcache[ex.receptor], molcache[ex.ligand], transform, gpu);
  }
  else
  {
    mol_info rec;
    mol_info lig;
    set_mol_info(root_folder+ex.receptor, rmap, 0, rec);
    set_mol_info(root_folder+ex.ligand, lmap, numReceptorTypes, lig);
    set_grid_minfo(data, rec, lig, transform, gpu);
  }
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data, const MolGridDataLayer<Dtype>::mol_info& recatoms,
  const MolGridDataLayer<Dtype>::mol_info& ligatoms, MolGridDataLayer<Dtype>::mol_transform& transform, bool gpu)
{
  //set grid values from mol info
  //first clear transform from the previous batch
  transform = mol_transform();

  //include receptor and ligand atoms
  transform.mol.append(recatoms);
  transform.mol.append(ligatoms);

  //set center to ligand center
  transform.mol.center = ligatoms.center;

  //figure out transformation
  transform.Q = quaternion(1,0,0,0);
  if(current_rotation == 0 && !randrotate)
    transform.Q = quaternion(0,0,0,0); //check real part to avoid mult
  rng_t* rng = caffe_rng();
  if (randrotate)
  {
    double d =  ((*rng)() - rng->min()) / double(rng->max());
    double r1 = ((*rng)() - rng->min()) / double(rng->max());
    double r2 = ((*rng)() - rng->min()) / double(rng->max());
    double r3 = ((*rng)() - rng->min()) / double(rng->max());
    transform.Q = quaternion(1, r1 / d, r2 / d, r3 / d);
  }

  transform.center[0] = transform.mol.center[0];
  transform.center[1] = transform.mol.center[1];
  transform.center[2] = transform.mol.center[2];
  if (randtranslate)
  {
    double offx = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    double offy = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    double offz = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    transform.center[0] += offx * randtranslate;
    transform.center[1] += offy * randtranslate;
    transform.center[2] += offz * randtranslate;
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

    gmaker.setAtomsGPU<Dtype>(natoms, gpu_gridatoms, gpu_gridwhich, transform.Q, numReceptorTypes+numLigandTypes, data);
  }
  else
  {
    Grids grids(data, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);
    gmaker.setAtomsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids);
  }
}


//return a string representation of the atom type(s) represented by index
//in map - this isn't particularly efficient, but is only for debug purposes
template <typename Dtype>
string MolGridDataLayer<Dtype>::getIndexName(const vector<int>& map, unsigned index) const
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
template<typename Dtype>
void MolGridDataLayer<Dtype>::outputDXGrid(std::ostream& out, Grids& grid, unsigned g, double scale) const
{
  unsigned n = dim;
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

  out.precision(10);
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
template<typename Dtype>
void MolGridDataLayer<Dtype>::dumpDiffDX(const std::string& prefix,
		Blob<Dtype>* top, double scale) const
{
	Grids grids(top->mutable_cpu_diff(),
			boost::extents[numReceptorTypes + numLigandTypes][dim][dim][dim]);
    CHECK_GT(mem_lig.atoms.size(),0) << "DX dump only works with in-memory ligand";
    CHECK_EQ(randrotate, false) << "DX dump requires no rotation";
	for (unsigned a = 0, na = numReceptorTypes; a < na; a++) {
		string name = getIndexName(rmap, a);
		string fname = prefix + "_rec_" + name + ".dx";
		ofstream out(fname.c_str());
		outputDXGrid(out, grids, a, scale);
	}
	for (unsigned a = 0, na = numLigandTypes; a < na; a++) {
			string name = getIndexName(lmap, a);
			string fname = prefix + "_lig_" + name + ".dx";
			ofstream out(fname.c_str());
			outputDXGrid(out, grids, numReceptorTypes+a, scale);
	}

}


template <typename Dtype>
void MolGridDataLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top)
{
  forward(bottom, top, false);
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::forward(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top, bool gpu)
{
  bool hasaffinity = this->layer_param_.molgrid_data_param().has_affinity();
  bool hasrmsd = this->layer_param_.molgrid_data_param().has_rmsd();

  Dtype *top_data = NULL;
  if(gpu)
    top_data = top[0]->mutable_gpu_data();
  else
    top_data = top[0]->mutable_cpu_data();

  //clear batch labels
  labels.clear();
  affinities.clear();
  rmsds.clear();
  unsigned batch_size = top_shape[0];

  //if in memory must be set programmatically
  if(inmem)
  {
    CHECK_GT(mem_rec.atoms.size(),0) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(mem_lig.atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(top_data, mem_rec, mem_lig, batch_transform[0], gpu); //TODO how do we know what batch position?
    if (num_rotations > 0) {
      current_rotation = (current_rotation+1)%num_rotations;
    }
  }
  else
  {
    //percent of batch from first data source
    unsigned dataswitch = batch_size;
    if (data2)
      dataswitch = batch_size*data_ratio/(data_ratio+1);

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
      set_grid_ex(top_data+offset, ex, *root, batch_transform[batch_idx], gpu);
      //NOTE: num_rotations not actually implemented!
    }


    if(gpu) {
      caffe_copy(labels.size(), &labels[0], top[1]->mutable_gpu_data());
      if(hasaffinity) {
    	  caffe_copy(affinities.size(), &affinities[0], top[2]->mutable_gpu_data());
    	  if(hasrmsd) {
    		  caffe_copy(rmsds.size(), &rmsds[0], top[3]->mutable_gpu_data());
    	  }
      } else if(hasrmsd) {
    	  caffe_copy(rmsds.size(), &rmsds[0], top[2]->mutable_gpu_data());
      }

    }
    else {
      caffe_copy(labels.size(), &labels[0], top[1]->mutable_cpu_data());
      if(hasaffinity) {
    	  caffe_copy(affinities.size(), &affinities[0], top[2]->mutable_cpu_data());
    	  if(hasrmsd) {
    		  caffe_copy(rmsds.size(), &rmsds[0], top[3]->mutable_cpu_data());
    	  }
      } else if(hasrmsd) {
    	  caffe_copy(rmsds.size(), &rmsds[0], top[2]->mutable_cpu_data());
      }
    }
  }
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
  backward(top, bottom, false);
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::backward(const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom,
    bool gpu)
{
  Dtype *diff = NULL;
  if(gpu)
    diff = top[0]->mutable_cpu_diff(); //TODO
  else
    diff = top[0]->mutable_cpu_diff();

  //propagate gradient grid onto atom positions
  unsigned batch_size = top_shape[0];
  for (int item_id = 0; item_id < batch_size; ++item_id) {

    int offset = item_id*example_size;
    Grids grids(diff+offset, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);

    mol_transform& transform = batch_transform[item_id];
    gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);
    gmaker.setAtomGradientsCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids,
        transform.mol.gradient);
  }
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
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
    Grids grids(diff+offset, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);

    mol_transform& transform = batch_transform[item_id];
    gmaker.setCenter(transform.center[0], transform.center[1], transform.center[2]);
    gmaker.setAtomRelevanceCPU(transform.mol.atoms, transform.mol.whichGrid, transform.Q, grids,
        transform.mol.gradient);
  }

  //float bottom_sum = 0.0;
  //for(int i = 0; i < bottom[0]->count(); i++)
  //{
  //        bottom_sum += bottom[0]->cpu_diff()[i];
  //}
  //std::cout << "MOLGRID BOTTOM: " << bottom_sum << '\n';


}



INSTANTIATE_CLASS(MolGridDataLayer);
REGISTER_LAYER_CLASS(MolGridData);

}  // namespace caffe
