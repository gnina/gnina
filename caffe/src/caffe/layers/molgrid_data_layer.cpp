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

  root_folder = this->layer_param_.molgrid_data_param().root_folder();
  balanced  = this->layer_param_.molgrid_data_param().balanced();
  paired  = this->layer_param_.molgrid_data_param().paired();
  num_rotations = this->layer_param_.molgrid_data_param().rotate();
  all_pos_ = actives_pos_ = decoys_pos_ = 0;
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

  if(root_folder.length() > 0 && root_folder[root_folder.length()-1] != '/')
    root_folder = root_folder + "/"; //make sure we have trailing slash

  int rmsdindex = -1;
  if(hasaffinity) {
	  rmsdindex = 3;
  }
  else if(hasrmsd) {
	  rmsdindex = 2;
  }

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

  if(!inmem)
  {
    // Read the file with labels and structure, each line is
    // label receptor_file ligand_file
    const string& source = this->layer_param_.molgrid_data_param().source();

    LOG(INFO) << "Opening file " << source;
    std::ifstream infile(source.c_str());
    CHECK((bool)infile) << "Could not load " << source;


    string line, recname, ligname;

    while (getline(infile, line)) {
      stringstream ex(line);
      int label = 0;
      double affinity = 0.0;
      double rmsd = 0.0;
      ligname = string();
      recname = string();
      //first the label
      ex >> label;
      if(hasaffinity)
    	  ex >> affinity;
      if(hasrmsd)
    	  ex >> rmsd;
      //receptor
      ex >> recname;
      //ligand
      ex >> ligname;

      CHECK(recname.length() > 0) << "Empty receptor, missing affinity/rmsd?";
      CHECK(ligname.length() > 0) << "Empty ligand, missing affinity/rmsd?";

      if(root_folder.length() > 0) recname = root_folder + recname;
      if(root_folder.length() > 0) ligname = root_folder + ligname;

      example lineex(label, affinity, rmsd, recname, ligname);
      all_.push_back(lineex);

      if(label) actives_.push_back(lineex);
      else decoys_.push_back(lineex);

      pairs_.add(lineex);
    }

    if (this->layer_param_.molgrid_data_param().shuffle()) {
      // randomly shuffle data
      LOG(INFO) << "Shuffling data";
      Shuffle();
    }
    LOG(INFO) << "A total of " << all_.size() << " examples.";

    // Check if we would need to randomly skip a few data points
    if (this->layer_param_.molgrid_data_param().rand_skip()) {
      unsigned int skip = caffe_rng_rand() %
          this->layer_param_.molgrid_data_param().rand_skip();
      LOG(INFO) << "Skipping first " << skip << " data points.";
      CHECK_GT(all_.size(), skip) << "Not enough points to skip";
      all_pos_ = skip;
      actives_pos_ = skip % actives_.size();
      decoys_pos_ = skip % decoys_.size();
    }

    if(balanced) {
      CHECK_GT(batch_size, 1) << "Batch size must be > 1 with balanced option.";
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

  // label, affinity, rmsds
  vector<int> label_shape(1, batch_size);
  top[1]->Reshape(label_shape);

  if(hasaffinity)
	  top[2]->Reshape(label_shape);
  if(hasrmsd)
	  top[rmsdindex]->Reshape(label_shape);

}

template <typename Dtype>
void MolGridDataLayer<Dtype>::Shuffle() {

    shuffle(actives_.begin(), actives_.end(), caffe_rng());
    shuffle(decoys_.begin(), decoys_.end(), caffe_rng());
    shuffle(all_.begin(), all_.end(), caffe_rng());
    pairs_.shuffle_pairs();
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
    }
    center /= cnt;
  }
  else
  {
  //read mol from file and set mol info (atom coords and grid positions)
  //types are mapped using atommap values plus offset
    OpenBabel::OBConversion conv;
    OBMol mol;
    CHECK(conv.ReadFile(&mol, root_folder + file)) << "Could not read " << file;

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
    }
    center /= cnt;
  }

  minfo.center = center;

}

template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_ex(Dtype *data, const MolGridDataLayer<Dtype>::example& ex,
    MolGridDataLayer<Dtype>::mol_transform& transform, bool gpu)
{
  //output grid values for provided example
  //cache atom info
  if(molcache.count(ex.receptor) == 0)
  {
    set_mol_info(ex.receptor, rmap, 0, molcache[ex.receptor]);
  }
  if(molcache.count(ex.ligand) == 0)
  {
    set_mol_info(ex.ligand, lmap, numReceptorTypes, molcache[ex.ligand]);
  }

  set_grid_minfo(data, molcache[ex.receptor], molcache[ex.ligand], transform, gpu);

}


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data, const MolGridDataLayer<Dtype>::mol_info& recatoms,
  const MolGridDataLayer<Dtype>::mol_info& ligatoms, MolGridDataLayer<Dtype>::mol_transform& transform, bool gpu)
{
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
  Dtype *data = NULL;
  if(gpu)
    data = top[0]->mutable_gpu_data();
  else
    data = top[0]->mutable_cpu_data();

  labels.resize(0);
  affinities.resize(0);
  rmsds.resize(0);
  unsigned batch_size = top_shape[0];
  //if in memory must be set programmatically
  if(inmem) {
    CHECK_GT(mem_rec.atoms.size(),0) << "Receptor not set in MolGridDataLayer";
    CHECK_GT(mem_lig.atoms.size(),0) << "Ligand not set in MolGridDataLayer";
    //memory is now available
    set_grid_minfo(data, mem_rec, mem_lig, batch_transform[0], gpu); //TODO how do we know what batch position?
    if (num_rotations > 0) {
      current_rotation = (current_rotation+1)%num_rotations;
    }
  }
  else {
	if(paired) {
		CHECK_EQ(batch_size % 2, 0) << "Paired input requires even batch size in MolGridDataLayer";
		example active;
		example decoy;
		unsigned npairs = batch_size/2;
		for(unsigned i = 0; i < npairs; i++) {
	        int offset = labels.size()*example_size;
			pairs_.next(active, decoy);

			//active
			labels.push_back(active.label);
			affinities.push_back(active.affinity);
			rmsds.push_back(active.rmsd);
			set_grid_ex(data+offset, active, batch_transform[2*i], gpu);

			//then decoy
			offset += example_size;
			labels.push_back(decoy.label);
			affinities.push_back(decoy.affinity);
			rmsds.push_back(decoy.rmsd);
			set_grid_ex(data+offset, decoy, batch_transform[2*i+1], gpu);
		}

	}
	else if(balanced) { //load equally from actives/decoys
      unsigned nactives = batch_size/2;

      CHECK_GT(actives_.size(), 0) << "Need non-zero number of actives for balanced input in MolGridDataLayer";
      CHECK_GT(decoys_.size(), 0) << "Need non-zero number of decoys for balanced input in MolGridDataLayer";

      int item_id = 0;
      unsigned asz = actives_.size();
      for (item_id = 0; item_id < nactives; ++item_id) {
        int offset = item_id*example_size;
        labels.push_back(1.0);
        affinities.push_back(actives_[actives_pos_].affinity);
        rmsds.push_back(actives_[actives_pos_].rmsd);

        set_grid_ex(data+offset, actives_[actives_pos_], batch_transform[item_id], gpu);

        actives_pos_++;
        if(actives_pos_ >= asz) {
          DLOG(INFO) << "Restarting actives data  from start.";
          actives_pos_ = 0;
          if (this->layer_param_.molgrid_data_param().shuffle()) {
            shuffle(actives_.begin(), actives_.end(), caffe_rng());
          }
          //this is less than ideal, gets reset for both actives and decoys
          if (num_rotations > 0) {
            current_rotation = (current_rotation+1)%num_rotations;
          }
        }
      }
      unsigned dsz = decoys_.size();
      for (; item_id < batch_size; ++item_id) {
        int offset = item_id*example_size;
        labels.push_back(0.0);
        affinities.push_back(decoys_[decoys_pos_].affinity);
        rmsds.push_back(decoys_[decoys_pos_].rmsd);
        set_grid_ex(data+offset, decoys_[decoys_pos_], batch_transform[item_id], gpu);

        decoys_pos_++;
        if(decoys_pos_ >= dsz) {
          DLOG(INFO) << "Restarting decoys data  from start.";
          decoys_pos_ = 0;
          if (this->layer_param_.molgrid_data_param().shuffle()) {
            shuffle(decoys_.begin(), decoys_.end(), caffe_rng());
          }
          if (num_rotations > 0) {
            current_rotation = (current_rotation+1)%num_rotations;
          }
        }
      }

    } else {
      //load from all
      unsigned sz = all_.size();
      for (int item_id = 0; item_id < batch_size; ++item_id) {
        int offset = item_id*example_size;
        labels.push_back(all_[all_pos_].label);
        affinities.push_back(all_[all_pos_].affinity);
        rmsds.push_back(all_[all_pos_].rmsd);
        set_grid_ex(data+offset, all_[all_pos_], batch_transform[item_id], gpu);

        all_pos_++;
        if(all_pos_ >= sz) {
          DLOG(INFO) << "Restarting data  from start.";
          all_pos_ = 0;
          if (this->layer_param_.molgrid_data_param().shuffle()) {
            shuffle(all_.begin(), all_.end(), caffe_rng());
          }
          if (num_rotations > 0) {
            current_rotation = (current_rotation+1)%num_rotations;
          }
        }
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
  //grid the examples

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
