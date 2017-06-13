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
#include "caffe/util/rng.hpp"

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
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::paired_examples::add(const MolGridDataLayer::example& ex)
{
	//have we seen this receptor before?
	if(recmap.count(ex.receptor) == 0) {
		//if not, assign index and resize vectors
		recmap[ex.receptor] = receptors.size();
		receptors.push_back(ex.receptor); //honestly, don't really need a vector as opposed to a counter..
		actives.resize(receptors.size());
		decoys.resize(receptors.size());
	}

	unsigned rindex = recmap[ex.receptor];

	//add to appropriate sub-vector
	if(ex.label != 0) {

		//if active, add to indices
		indices.push_back(make_pair(rindex, actives[rindex].size()));
		actives[rindex].push_back(ex);
	}
	else {
		decoys[rindex].first = 0;
		decoys[rindex].second.push_back(ex);
	}
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::paired_examples::shuffle_pairs() //randomize - only necessary at start
{
    shuffle(indices.begin(), indices.end(), caffe_rng());
    //shuffle decoys
    for(unsigned i = 0, n = decoys.size(); i < n; i++) {
    	shuffle(decoys[i].second.begin(), decoys[i].second.end(), caffe_rng());
    }
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::paired_examples::next(example& active, example& decoy)
{
	assert(indices.size() > 0);
	if(curr_index >= indices.size()) {
		//need to wrap around and re-shuffle
		curr_index = 0;
		shuffle(indices.begin(), indices.end(), caffe_rng());
		//decoys get shuffled on demand
	}

	unsigned rindex = indices[curr_index].first;

	while(decoys[rindex].second.size() == 0) {
		//skip targets with empty decoys
		curr_index++;
		if(curr_index == indices.size())
			break;
		rindex = indices[curr_index].first;
	}

	//allow one wrap
	if(curr_index == indices.size()) {
		curr_index = 0;
		rindex = indices[curr_index].first;
		while(decoys[rindex].second.size() == 0) {
			curr_index++; 
			CHECK_LT(curr_index, indices.size()) << "No decoy examples for pairing.";
			rindex = indices[curr_index].first;
		}
	}
	unsigned aindex = indices[curr_index].second;
	assert(rindex < actives.size());
	assert(aindex < actives[rindex].size());
	active = actives[rindex][aindex];

	//now get decoy, reshuffling if necessary
	assert(rindex < decoys.size());
	unsigned dindex = decoys[rindex].first;
	vector<example>& decvec = decoys[rindex].second;
	if(dindex >= decvec.size()) {
		dindex = 0;
		shuffle(decvec.begin(), decvec.end(), caffe_rng());
	}
	decoy = decvec[dindex];
	//increment indices
	decoys[rindex].first++;
	curr_index++;

}

template <typename Dtype>
MolGridDataLayer<Dtype>::example::example(string line, bool hasaffinity, bool hasrmsd)
  : label(0), affinity(0.0), rmsd(0.0)
{
  stringstream stream(line);

  //first the label
  stream >> label;
  if(hasaffinity)
   stream >> affinity;
  if(hasrmsd)
   stream >> rmsd;
  //receptor
  stream >> receptor;
  //ligand
  stream >> ligand;

  CHECK(receptor.length() > 0) << "Empty receptor, missing affinity/rmsd?";
  CHECK(ligand.length() > 0) << "Empty ligand, missing affinity/rmsd?";
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::examples::add(const example& ex)
{
  if (store_all)
  {
    all.push_back(ex);
  }
  if (store_actives_decoys)
  {
    if (ex.label)
      actives.push_back(ex);
    else
      decoys.push_back(ex);
  }
  if (store_pairs)
  {
    pairs.add(ex);
  }
  count++;
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::examples::shuffle_()
{
  shuffle(all.begin(), all.end(), caffe_rng());
  shuffle(actives.begin(), actives.end(), caffe_rng());
  shuffle(decoys.begin(), decoys.end(), caffe_rng());
  pairs.shuffle_pairs();
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::examples::next(example& ex)
{
  ex = all[all_index];
  all_index++;
  if (all_index >= all.size())
  {
    all_index = 0;
    if (shuffle_on_wrap)
      shuffle(all.begin(), all.end(), caffe_rng());
  }
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::examples::next_active(example& ex)
{
  ex = actives[actives_index];
  actives_index++;
  if (actives_index >= actives.size())
  {
    actives_index = 0;
    if (shuffle_on_wrap)
      shuffle(actives.begin(), actives.end(), caffe_rng());
  }
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::examples::next_decoy(example& ex)
{
  ex = decoys[decoys_index];
  decoys_index++;
  if (decoys_index >= decoys.size())
  {
    decoys_index = 0;
    if (shuffle_on_wrap)
      shuffle(decoys.begin(), decoys.end(), caffe_rng());
  }
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


//read in structure input and atom type maps
template <typename Dtype>
void MolGridDataLayer<Dtype>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {

  string root_folder = this->layer_param_.molgrid_data_param().root_folder();
  balanced  = this->layer_param_.molgrid_data_param().balanced();
  paired  = this->layer_param_.molgrid_data_param().paired();
  num_rotations = this->layer_param_.molgrid_data_param().rotate();
  inmem = this->layer_param_.molgrid_data_param().inmemory();
  dimension = this->layer_param_.molgrid_data_param().dimension();
  resolution = this->layer_param_.molgrid_data_param().resolution();
  binary = this->layer_param_.molgrid_data_param().binary_occupancy();
  randtranslate = this->layer_param_.molgrid_data_param().random_translate();
  randrotate = this->layer_param_.molgrid_data_param().random_rotation();
  radiusmultiple = this->layer_param_.molgrid_data_param().radius_multiple();
  fixedradius = this->layer_param_.molgrid_data_param().fixed_radius();
  bool hasaffinity = this->layer_param_.molgrid_data_param().has_affinity();
  bool hasrmsd = this->layer_param_.molgrid_data_param().has_rmsd();
  data_ratio = this->layer_param_.molgrid_data_param().source_ratio();
  string root_folder2 = this->layer_param_.molgrid_data_param().root_folder2();
  bool shuffle = this->layer_param_.molgrid_data_param().shuffle();

  //make sure root folder(s) have trailing slash
  if (root_folder.length() > 0 && root_folder[root_folder.length()-1] != '/')
    root_folder = root_folder + "/";
  if (root_folder2.length() > 0 && root_folder2[root_folder2.length()-1] != '/')
    root_folder2 = root_folder2 + "/";

  if(binary) radiusmultiple = 1.0;

  gmaker.initialize(resolution, dimension, radiusmultiple, binary);

  dim = round(dimension/resolution)+1; //number of grid points on a side
  numgridpoints = dim*dim*dim;
  if(numgridpoints % 512 != 0)
    LOG(INFO) << "Total number of grid points (" << numgridpoints << ") is not evenly divisible by 512.";

  //shape must come from parameters
  const int batch_size = this->layer_param_.molgrid_data_param().batch_size();
  CHECK_GT(batch_size, 0) << "Positive batch size required";

  //keep track of atoms and transformations for each example in batch
  batch_transform.resize(batch_size);

  CHECK_LE(inmem + paired + balanced, 1) << "Only one of inmemory, paired, and balanced can be set";

  if(!inmem)
  {
    const string& source = this->layer_param_.molgrid_data_param().source();
    const string& source2 = this->layer_param_.molgrid_data_param().source2();

    CHECK_GT(source.length(), 0) << "No data source file provided";

    // Read source file(s) with labels and structures,
    // each line is label [affinity] [rmsd] receptor_file ligand_file
    LOG(INFO) << "Opening file " << source;
    std::ifstream infile(source.c_str());
    CHECK((bool)infile) << "Could not open " << source;

    bool all = !(balanced || paired);
    data = examples(all, balanced, paired);
    data.root_folder = root_folder;
    data.shuffle_on_wrap = shuffle;

    string line;
    while (getline(infile, line))
    {
      example ex(line, hasaffinity, hasrmsd);
      data.add(ex);
    }

    if (source2.length() > 0)
    {
      CHECK_GE(data_ratio, 0) << "Must provide non-negative ratio for two data sources";

      LOG(INFO) << "Opening file " << source2;
      std::ifstream infile(source2.c_str());
      CHECK((bool)infile) << "Could not open " << source2;

      data2 = examples(all, balanced, paired);
      data2.root_folder = root_folder2;
      data2.shuffle_on_wrap = shuffle;

      while (getline(infile, line))
      {
        example ex(line, hasaffinity, hasrmsd);
        data2.add(ex);
      }

      LOG(INFO) << "Combining at ratio " << data_ratio;
      two_data_sources = true;
    }
    else
    {
      CHECK_EQ(data_ratio, -1) << "Need two data sources to use ratio";
      two_data_sources = false;
    }

    if (shuffle)
    {
      // randomly shuffle data
      LOG(INFO) << "Shuffling data";
      data.shuffle_();
      if (two_data_sources)
        data2.shuffle_();
    }

    LOG(INFO) << "Total examples: " << data.count + data2.count;

    // Check if we would need to randomly skip a few data points
    if (this->layer_param_.molgrid_data_param().rand_skip())
    {
      unsigned int skip = caffe_rng_rand() %
          this->layer_param_.molgrid_data_param().rand_skip();

      LOG(INFO) << "Skipping first " << skip << " data points.";
      CHECK_GT(data.count, skip) << "Not enough points to skip in " << source;

      data.all_index = skip % data.all.size();
      data.actives_index = skip % data.actives.size();
      data.decoys_index = skip % data.decoys.size();

      if (two_data_sources)
      {
        CHECK_GT(data2.count, skip) << "Not enough points to skip in " << source2;

        data2.all_index = skip % data2.all.size();
        data2.actives_index = skip % data2.actives.size();
        data2.decoys_index = skip % data2.decoys.size();
      }
    }

    if (balanced)
    {
      CHECK_GT(batch_size, 1) << "Batch size must be > 1 to balance actives/decoys";

      CHECK_GT(data.actives.size(), 0) << "Need non-zero number of actives to balance " << source;
      CHECK_GT(data.decoys.size(), 0) << "Need non-zero number of decoys to balance " << source;

      if (two_data_sources)
      {
        CHECK_GT(data.actives.size(), 0) << "Need non-zero number of actives to balance " << source2;
        CHECK_GT(data.decoys.size(), 0) << "Need non-zero number of decoys to balance " << source2;
      }
    }
  }

  //initialize atom type maps
  string recmap = this->layer_param_.molgrid_data_param().recmap();
  string ligmap = this->layer_param_.molgrid_data_param().ligmap();

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
  else
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

  //percent of batch from first data source
  float pct_data = 1.0;
  if (two_data_sources)
    pct_data = data_ratio/(data_ratio+1);

  //if in memory must be set programmatically
  if(inmem) {
    CHECK_GT(mem_rec.atoms.size(),0) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(mem_lig.atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(top_data, mem_rec, mem_lig, batch_transform[0], gpu); //TODO how do we know what batch position?
    if (num_rotations > 0) {
      current_rotation = (current_rotation+1)%num_rotations;
    }
  }
  else {
	if(paired) {
		CHECK(!two_data_sources) << "Paired input not implemented for two data sources"; //TODO
		CHECK_EQ(batch_size % 2, 0) << "Paired input requires even batch size in MolGridDataLayer";
		example active;
		example decoy;
		unsigned npairs = batch_size/2;
		for(unsigned i = 0; i < npairs; i++) {
	        int offset = labels.size()*example_size;
			data.pairs.next(active, decoy);

			//active
			labels.push_back(active.label);
			affinities.push_back(active.affinity);
			rmsds.push_back(active.rmsd);
			set_grid_ex(top_data+offset, active, data.root_folder, batch_transform[2*i], gpu);

			//then decoy
			offset += example_size;
			labels.push_back(decoy.label);
			affinities.push_back(decoy.affinity);
			rmsds.push_back(decoy.rmsd);
			set_grid_ex(top_data+offset, decoy, data.root_folder, batch_transform[2*i+1], gpu);
		}

	}
    else if (balanced) //load equally from actives and decoys
    {
      for (int batch_idx = 0; batch_idx < batch_size; ++batch_idx)
      {
        example ex;
        string root_folder;
        if (batch_idx < batch_size/2) //actives
        {
          if (batch_idx < pct_data*batch_size/2)
          {
            data.next_active(ex);
            root_folder = data.root_folder;
          }
          else
          {
            data2.next_active(ex);
            root_folder = data2.root_folder;
          }
        }
        else //decoys
        {
          if (batch_idx < (pct_data+1)*batch_size/2)
          {
            data.next_decoy(ex);
            root_folder = data.root_folder;
          }
          else
          {
            data2.next_decoy(ex);
            root_folder = data2.root_folder;
          }
        }

        labels.push_back(ex.label);
        affinities.push_back(ex.affinity);
        rmsds.push_back(ex.rmsd);

        int offset = batch_idx*example_size;
        set_grid_ex(top_data+offset, ex, root_folder, batch_transform[batch_idx], gpu);

        //this is less than ideal, gets reset for both actives and decoys
        //TODO increment current_rotation when data wraps around
      }
    } else //load from all
    {
      for (int batch_idx = 0; batch_idx < batch_size; ++batch_idx)
      {
        example ex;
        string root_folder;
        if (batch_idx < pct_data*batch_size)
        {
          data.next(ex);
          root_folder = data.root_folder;
        }
        else
        {
          data2.next(ex);
          root_folder= data2.root_folder;
        }

        labels.push_back(ex.label);
        affinities.push_back(ex.affinity);
        rmsds.push_back(ex.rmsd);

        int offset = batch_idx*example_size;
        set_grid_ex(top_data+offset, ex, root_folder, batch_transform[batch_idx], gpu);

        //TODO increment current_rotation when data wraps around
      }
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

  float top_sum = 0.0;

  //std::cout << "MOLGRID TOP:";
  for(int i = 0; i < top[0]->count(); i++)
  {
          //std::cout << top[0]->cpu_diff()[i] << "|";
          top_sum += top[0]->cpu_diff()[i];
  }
  //std::cout << '\n';
  std::cout << "MOLGRID TOP: " << top_sum << '\n';

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
