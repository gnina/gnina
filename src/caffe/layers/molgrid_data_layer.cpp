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
}



//create a mapping from atom type ids to a unique id given a file specifying
//what types we care about (anything missing is ignored); if multiple types are
//on the same line, they are merged, if the file isn't specified, use default mapping
//return total number of types
//map is indexed by smina_atom_type, maps to -1 if type should be ignored
static unsigned createAtomTypeMap(const string& fname, vector<int>& map)
{
	using namespace std;
	using namespace boost::algorithm;
  map.assign(smina_atom_type::NumTypes, -1);
  CHECK(fname.size() > 0) <<  "Map file not specified\n";

	unsigned cnt = 0;
	ifstream in(fname.c_str());

	CHECK(in) << "Could not open " << fname << "\n";

	string line;
	while (getline(in, line))
	{
		vector<string> types;
		split(types, line, is_any_of("\t \n"));
		for (unsigned i = 0, n = types.size(); i < n; i++)
		{
			const string& name = types[i];
			smt t = string_to_smina_type(name);
			if (t < smina_atom_type::NumTypes) //valid
			{
				map[t] = cnt;
			}
			else if (name.size() > 0) //this ignores consecutive delimiters
			{
				cerr << "Invalid atom type " << name << "\n";
				exit(-1);
			}
		}
		if (types.size() > 0)
			cnt++;
	}
	return cnt;
}

static unsigned createDefaultMap(const char *names[], vector<int>& map)
{
  map.assign(smina_atom_type::NumTypes, -1);
  const char **nameptr = names;
  unsigned cnt = 0;
  while (*nameptr != NULL)
  {
    string name(*nameptr);
    //note that if we every start using merged atom types by default
    //this code will have to be updated
    smt t = string_to_smina_type(name);
    CHECK_LT(t, smina_atom_type::NumTypes) << "Invalid atom type " << name << "\n";
    map[t] = cnt;
    cnt++;
    nameptr++;
  }
  return cnt;
}

//initialize default receptor/ligand maps
//these were determined by evaluating how common various atoms types are
static unsigned createDefaultRecMap(vector<int>& map)
{
  const char *names[] =
  { "AliphaticCarbonXSHydrophobe",
      "AliphaticCarbonXSNonHydrophobe",
      "AromaticCarbonXSHydrophobe",
      "AromaticCarbonXSNonHydrophobe",
      "Calcium",
      "Iron",
      "Magnesium",
      "Nitrogen",
      "NitrogenXSAcceptor",
      "NitrogenXSDonor",
      "NitrogenXSDonorAcceptor",
      "OxygenXSAcceptor",
      "OxygenXSDonorAcceptor",
      "Phosphorus",
      "Sulfur",
      "Zinc", NULL };

  return createDefaultMap(names, map);
}

static unsigned createDefaultLigMap(vector<int>& map)
{
  const char *names[] =
  { "AliphaticCarbonXSHydrophobe",
      "AliphaticCarbonXSNonHydrophobe",
      "AromaticCarbonXSHydrophobe",
      "AromaticCarbonXSNonHydrophobe",
      "Bromine",
      "Chlorine",
      "Fluorine",
      "Nitrogen",
      "NitrogenXSAcceptor",
      "NitrogenXSDonor",
      "NitrogenXSDonorAcceptor",
      "Oxygen",
      "OxygenXSAcceptor",
      "OxygenXSDonorAcceptor",
      "Phosphorus",
      "Sulfur",
      "SulfurAcceptor",
      "Iodine",
      NULL };
  return createDefaultMap(names, map);
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

  if(root_folder.length() > 0 && root_folder[root_folder.length()-1] != '/')
    root_folder = root_folder + "/"; //make sure we have trailing slash

  if(binary) radiusmultiple = 1.0;

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
      //first the label
      ex >> label;
      //receptor
      ex >> recname;
      //ligand
      ex >> ligname;

      example lineex(label, recname, ligname);
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
    numReceptorTypes = createDefaultRecMap(rmap);
  else
    numReceptorTypes = createAtomTypeMap(recmap, rmap);

  if (ligmap.size() == 0)
    numLigandTypes = createDefaultLigMap(lmap);
  else
    numLigandTypes = createAtomTypeMap(ligmap, lmap);


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

  // label
  vector<int> label_shape(1, batch_size);
  top[1]->Reshape(label_shape);

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
void MolGridDataLayer<Dtype>::memoryIsSet()
{
  boost::unique_lock<boost::mutex> lock(mem_mutex);
  unsigned batch_size = top_shape[0];
  unsigned add = 1;

  if(num_rotations > 0) {
    CHECK_LE(batch_size, num_rotations);
    CHECK_EQ(num_rotations % batch_size, 0);
    add = num_rotations/batch_size;
  }
  data_avail += add;

  mem_cond.notify_one();
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::set_mol_info(const string& file, const vector<int>& atommap,
    unsigned mapoffset, mol_info& minfo)
{
  using namespace OpenBabel;
  //read mol from file and set mol info (atom coords and grid positions)
  //types are mapped using atommap values plus offset
  OpenBabel::OBConversion conv;
  OBMol mol;
  CHECK(conv.ReadFile(&mol, root_folder + file)) << "Could not read " << file;

  minfo.atoms.clear(); minfo.atoms.reserve(mol.NumHvyAtoms());
  minfo.whichGrid.clear(); minfo.whichGrid.reserve(mol.NumHvyAtoms());
  vec center(0,0,0);

  int cnt = 0;
  FOR_ATOMS_OF_MOL(a, mol)
  {
    smt t = obatom_to_smina_type(*a);
    int index = atommap[t];
    if(index >= 0)
    {
      cnt++;
      atom_info ainfo;
      ainfo.coord[0] = a->x();
      ainfo.coord[1] = a->y();
      ainfo.coord[2]  = a->z();
      ainfo.radius = xs_radius(t);

      minfo.atoms.push_back(ainfo);
      minfo.whichGrid.push_back(index+mapoffset);
      center += ainfo.coord;
    }
  }
  center /= cnt;

  minfo.setCenter(dimension, center);

}

template <typename Dtype>
void MolGridDataLayer<Dtype>::zeroGrids(Grids& grids)
{
  std::fill(grids.data(), grids.data() + grids.num_elements(), 0.0);
}



//return the range of grid points spanned from c-r to c+r within dim
template <typename Dtype>
pair<unsigned, unsigned>
MolGridDataLayer<Dtype>::getrange(const pair<float, float>& d, double c,
    double r)
{
  pair<unsigned, unsigned> ret(0, 0);
  double low = c - r - d.first;
  if (low > 0)
  {
    ret.first = floor(low / resolution);
  }

  double high = c + r - d.first;
  if (high > 0) //otherwise zero
  {
    ret.second = std::min(dim, (unsigned) ceil(high / resolution));
  }
  return ret;
}


//return the occupancy for atom a at point x,y,z
template <typename Dtype>
Dtype MolGridDataLayer<Dtype>::calcPoint(const vec& coords, double ar, const vec& pt)
{
  Dtype rsq = (pt - coords).norm_sqr();
  if (binary)
  {
    //is point within radius?
    if (rsq < ar * ar)
      return 1.0;
    else
      return 0.0;
  }
  else
  {
    //for non binary we want a gaussian were 2 std occurs at the radius
    //after which which switch to a quadratic
    //the quadratic is to fit to have both the same value and first order
    //derivative at the cross over point and a value and derivative of zero
    //at 1.5*radius
    double dist = sqrt(rsq);
    if (dist >= ar * radiusmultiple)
    {
      return 0.0;
    }
    else if (dist <= ar)
    {
      //return gaussian
      Dtype h = 0.5 * ar;
      Dtype ex = -dist * dist / (2 * h * h);
      return exp(ex);
    }
    else //return quadratic
    {
      Dtype h = 0.5 * ar;
      Dtype eval = 1.0 / (M_E * M_E); //e^(-2)
      Dtype q = dist * dist * eval / (h * h) - 6.0 * eval * dist / h
          + 9.0 * eval;
      return q;
    }
  }
  return 0.0;
}

//set the relevant grid points for a
//return false if atom not in grid
template <typename Dtype>
void MolGridDataLayer<Dtype>::set_atom(const mol_info& mol, const atom_info& atom, int whichgrid, const quaternion& Q, Grids& grids)
{
  vec coords;
  if (Q.real() != 0)
  { //apply rotation
    vec tpt = atom.coord - mol.center;
    quaternion p(0, tpt[0], tpt[1], tpt[2]);
    p = Q * p * (conj(Q) / norm(Q));

    vec newpt(p.R_component_2(), p.R_component_3(), p.R_component_4());
    newpt += mol.center;

    coords = newpt;
  }
  else
  {
    coords = atom.coord;
  }

  double r = atom.radius * radiusmultiple;
  vector<pair<unsigned, unsigned> > ranges;
  for (unsigned i = 0; i < 3; i++)
  {
    ranges.push_back(getrange(mol.dims[i], coords[i], r));
  }
  //for every grid point possibly overlapped by this atom
  for (unsigned i = ranges[0].first, iend = ranges[0].second; i < iend; i++)
  {
    for (unsigned j = ranges[1].first, jend = ranges[1].second; j < jend;
        j++)
    {
      for (unsigned k = ranges[2].first, kend = ranges[2].second;
          k < kend; k++)
      {
        Dtype x = mol.dims[0].first + i * resolution;
        Dtype y = mol.dims[1].first + j * resolution;
        Dtype z = mol.dims[2].first + k * resolution;
        Dtype val = calcPoint(coords, atom.radius, vec(x, y, z));

        if (binary)
        {
          if (val != 0)
            grids[whichgrid][i][j][k] = 1.0; //don't add, just 1 or 0
        }
        else
          grids[whichgrid][i][j][k] += val;

      }
    }
  }
}

//set the relevant grid points for passed info
template <typename Dtype>
void MolGridDataLayer<Dtype>::set_atoms(const mol_info& mol, const quaternion& Q, Grids& grids)
{
  zeroGrids(grids);
  assert(mol.atoms.size() == mol.whichGrid.size());
  for (unsigned i = 0, n = mol.atoms.size(); i < n; i++)
  {
    int pos = mol.whichGrid[i];
    if (pos >= 0)
      set_atom(mol, mol.atoms[i], pos, Q, grids);
  }
}


template <typename Dtype>
void MolGridDataLayer<Dtype>::set_grid(Dtype *data, MolGridDataLayer<Dtype>::example ex, bool gpu)
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
  mol_info gridatoms; //includes receptor and ligand
  gridatoms.append(molcache[ex.receptor]);

  const mol_info linfo = molcache[ex.ligand];
  gridatoms.append(linfo);
  //set center to ligand center
  gridatoms.center = linfo.center;

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

  if (randtranslate)
  {
    double offx = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    double offy = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    double offz = ((*rng)() - rng->min()) / double(rng->max()/2.0)-1.0;
    gridatoms.center[0] += offx * randtranslate;
    gridatoms.center[1] += offy * randtranslate;
    gridatoms.center[2] += offz * randtranslate;
  }

  if(current_rotation > 0) {
    Q *= axial_quaternion();
  }

  //compute grid from atom info arrays
  if(gpu)
  {
    abort();
  }
  else
  {
    Grids grids(data, boost::extents[numReceptorTypes+numLigandTypes][dim][dim][dim]);
    set_atoms(gridatoms, Q, grids);
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
  labels.resize(0);
  unsigned batch_size = top_shape[0];
  //if in memory must be set programmatically
  if(inmem) {
    boost::unique_lock<boost::mutex> lock(mem_mutex);
    while(data_avail == 0)
    {
      mem_cond.wait(lock);
    }
    data_avail--;

    //memory is now available
    //TODO - something
    abort();
  }
  else {
    if(balanced) { //load equally from actives/decoys
      unsigned nactives = batch_size/2;

      int item_id = 0;
      unsigned asz = actives_.size();
      for (item_id = 0; item_id < nactives; ++item_id) {
        int offset = item_id*example_size;
        labels.push_back(1.0);

        set_grid(top[0]->mutable_cpu_data()+offset, actives_[actives_pos_], gpu);

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

        set_grid(top[0]->mutable_cpu_data()+offset, decoys_[decoys_pos_], gpu);

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

        set_grid(top[0]->mutable_cpu_data()+offset, all_[all_pos_], gpu);

        all_pos_++;
        if(all_pos_ >= sz) {
          DLOG(INFO) << "Restarting data  from start.";
          all_pos_ = 0;
          if (this->layer_param_.ndim_data_param().shuffle()) {
            shuffle(all_.begin(), all_.end(), caffe_rng());
          }
          if (num_rotations > 0) {
            current_rotation = (current_rotation+1)%num_rotations;
          }
        }
      }
    }

    if(gpu)
      caffe_copy(labels.size(), &labels[0], top[1]->mutable_gpu_data());
    else
      caffe_copy(labels.size(), &labels[0], top[1]->mutable_cpu_data());
  }
  //grid the examples

}


INSTANTIATE_CLASS(MolGridDataLayer);
REGISTER_LAYER_CLASS(MolGridData);

}  // namespace caffe
