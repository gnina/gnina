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

  dim = round(dimension/resolution)+1; //number of grid points on a size
  numgridpoints = dim*dim*dim;
  if(numgridpoints % 512 != 0)
    LOG(INFO) << "Total number of grid points (" << numgridpoints << ") is not evenly divisible by 512.";

  //shape must come from parameters
  const int batch_size = this->layer_param_.molgrid_data_param().batch_size();
  CHECK_GT(batch_size, 0) << "Positive batch size required";

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
  string recmap =  this->layer_param_.molgrid_data_param().recmap();
  string ligmap =  this->layer_param_.molgrid_data_param().ligmap();

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
        ainfo.z  = atom.z;
        if(fixedradius <= 0)
        	ainfo.w = xs_radius(t);
        else
        	ainfo.w = fixedradius;

        minfo.atoms.push_back(ainfo);
        minfo.whichGrid.push_back(index+mapoffset);
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

        minfo.atoms.push_back(ainfo);
        minfo.whichGrid.push_back(index+mapoffset);
        center += vec(a->x(),a->y(),a->z());
      }
    }
    center /= cnt;
  }

  minfo.center = center;

}

template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_ex(Dtype *data, const MolGridDataLayer<Dtype>::example& ex, bool gpu)
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

  set_grid_minfo(data, molcache[ex.receptor], molcache[ex.ligand], gpu);

}


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid_minfo(Dtype *data, const MolGridDataLayer<Dtype>::mol_info& recatoms, const MolGridDataLayer<Dtype>::mol_info& ligatoms, bool gpu)
{
  mol_info gridatoms; //includes receptor and ligand
  gridatoms.append(recatoms);
  gridatoms.append(ligatoms);
  //set center to ligand center
  gridatoms.center = ligatoms.center;

  //figure out transformation
  quaternion Q(1,0,0,0);
  if(current_rotation == 0 && !randrotate)
    Q = quaternion(0,0,0,0); //check real part to avoid mult
  rng_t* rng = caffe_rng();
  if (randrotate)
  {
    double d =  ((*rng)() - rng->min()) / double(rng->max());
    double r1 = ((*rng)() - rng->min()) / double(rng->max());
    double r2 = ((*rng)() - rng->min()) / double(rng->max());
    double r3 = ((*rng)() - rng->min()) / double(rng->max());
    Q = quaternion(1, r1 / d, r2 / d, r3 / d);
  }

  vec center(gridatoms.center);
  if (randtranslate)
  {
    double offx = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    double offy = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    double offz = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    center[0] += offx * randtranslate;
    center[1] += offy * randtranslate;
    center[2] += offz * randtranslate;
  }

  if(current_rotation > 0) {
    Q *= axial_quaternion();
  }

  gmaker.setCenter(center[0], center[1], center[2]);
  
  //compute grid from atom info arrays
  if(gpu)
  {
    unsigned natoms = gridatoms.atoms.size();
    allocateGPUMem(natoms);
    CUDA_CHECK(cudaMemcpy(gpu_gridatoms, &gridatoms.atoms[0], natoms*sizeof(float4),cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(gpu_gridwhich, &gridatoms.whichGrid[0], natoms*sizeof(short),cudaMemcpyHostToDevice));

    gmaker.setAtomsGPU<Dtype>(natoms, gpu_gridatoms, gpu_gridwhich, Q, numReceptorTypes+numLigandTypes, data);
  }
  else
  {
    Grids grids(data, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);
    gmaker.setAtomsCPU(gridatoms.atoms, gridatoms.whichGrid,Q, grids);
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
    set_grid_minfo(data, mem_rec, mem_lig, gpu);
    if (num_rotations > 0) {
      current_rotation = (current_rotation+1)%num_rotations;
    }
  }
  else {
    if(balanced) { //load equally from actives/decoys
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

        set_grid_ex(data+offset, actives_[actives_pos_], gpu);

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
        set_grid_ex(data+offset, decoys_[decoys_pos_], gpu);

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
        set_grid_ex(data+offset, all_[all_pos_], gpu);

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


INSTANTIATE_CLASS(MolGridDataLayer);
REGISTER_LAYER_CLASS(MolGridData);

}  // namespace caffe
