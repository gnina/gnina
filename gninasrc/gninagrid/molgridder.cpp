/*
 * molgridder.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: dkoes
 */

#include "molgridder.h"
#include <libmolgrid/cartesian_grid.h>
#include <boost/timer/timer.hpp>
#include <boost/lexical_cast.hpp>


using namespace std;
using namespace libmolgrid;

MolGridder::MolGridder(const gridoptions& opt) :
    dimension(opt.dim), resolution(opt.res),
    random_rotate(opt.randrotate), random_translate(opt.randtranslate),
    gmaker(opt.res, opt.dim), gpu(opt.gpu) {

  ex.sets.resize(2); //ligand and receptor

  //setup typers
  rectyper = std::make_shared<FileMappedGninaTyper>(defaultGninaReceptorTyper);
  ligtyper = std::make_shared<FileMappedGninaTyper>(defaultGninaLigandTyper);

  if (opt.recmap.size() > 0) {
    rectyper = std::make_shared<FileMappedGninaTyper>(opt.recmap);
  }
  if (opt.ligmap.size() > 0) {
    ligtyper = std::make_shared<FileMappedGninaTyper>(opt.ligmap);
  }

  if (opt.examplegrid.size() > 0) set_from_example(opt.examplegrid);
  if (opt.usergrids.size()) set_usergrids(opt.usergrids);

  gmaker.set_resolution(resolution);
  gmaker.set_dimension(dimension);
  gmaker.set_binary(opt.binary);

  float3 dims = gmaker.get_grid_dims();
  grid = MGrid4f(rectyper->num_types()+ligtyper->num_types()+usergrids.size(), dims.x, dims.y, dims.z);
  tee log(true);
  FlexInfo finfo(log); //dummy
  mols.create_init_model(opt.receptorfile, "", finfo, log);
  mols.setInputFile(opt.ligandfile);
  N = 1 + round(dimension / resolution);

  ex.sets[0].set_num_types(rectyper->num_types());
  ex.sets[1].set_num_types(ligtyper->num_types());
  setReceptor(mols.getInitModel());

  if(opt.separate) setGrid(gpu);

}

//set receptor coordinate set from model
void MolGridder::setReceptor(const model& m) {
  const atomv& atoms = m.get_fixed_atoms();
  vector<float3> coords; coords.reserve(atoms.size());
  vector<int> t; t.reserve(atoms.size());
  vector<float> r; t.reserve(atoms.size());

  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    const atom& a = atoms[i];
    auto t_r = rectyper->get_int_type(a.sm);
    if (t_r.first >= 0) {
      coords.push_back(gfloat3(a.coords));
      t.push_back(t_r.first);
      r.push_back(t_r.second);
    }
  }

  ex.sets[0] = CoordinateSet(coords, t, r, rectyper->num_types());
}

void MolGridder::setLigand(const model& m) {
  const atomv& atoms = m.get_movable_atoms();
  assert(atoms.size() == m.coordinates().size());

  vector<float3> coords; coords.reserve(atoms.size());
  vector<int> t; t.reserve(atoms.size());
  vector<float> r; t.reserve(atoms.size());

  for (unsigned i = 0, n = atoms.size(); i < n; i++) {
    const atom& a = atoms[i];
    auto t_r = ligtyper->get_int_type(a.sm);
    if (t_r.first >= 0) {
      coords.push_back(gfloat3(m.coordinates()[i]));
      t.push_back(t_r.first);
      r.push_back(t_r.second);
    }
  }

  ex.sets[1] = CoordinateSet(coords, t, r, ligtyper->num_types());
}

//convert ex to a grid
void MolGridder::setGrid(bool use_gpu) {
  if(!center_set) {
    //get center from lig
    center = ex.sets[1].center();
  }

  if(random_translate > 0 || random_rotate) {
    //update transform
    current_transform = Transform(center, random_translate, random_rotate);
  } else {
    current_transform.set_rotation_center(center);
  }

  if(usergrids.size() > 0) { //not particularly optimized
    //copy into first so many channels
    unsigned n = usergrids.size();

    for(unsigned i = 0; i < n; i++) {
      grid[i].copyFrom(usergrids[i].cpu());
    }

    size_t offset = grid[0].size()*n;
    unsigned channels = rectyper->num_types()+ligtyper->num_types();
    if(use_gpu) {
      Grid4fCUDA g(grid.gpu().data()+offset, channels, N, N, N);
      gmaker.forward(ex, current_transform, g);
    } else {
      Grid4f g(grid.cpu().data()+offset, channels, N, N, N);
      gmaker.forward(ex, current_transform, g);
    }
  }
  else { //no user grids
    if(use_gpu) {
      gmaker.forward(ex, current_transform, grid.gpu());
    } else {
      gmaker.forward(ex, current_transform, grid.cpu());
    }
  }
}

void MolGridder::cpuSetGridCheck() {
  //assume current grid is right
  MGrid4f saved = grid.clone();

  //for recompute, apply same transformation
  bool saverot = random_rotate;
  float savetrans = random_translate;
  setGrid(false);

  random_rotate = saverot;
  random_translate = savetrans;

  //now compare
  if(saved.size() != grid.size()) {
    cerr << "Different sized grids in compare\n";
    exit(-1);
  }
  float* a = saved.data();
  float* b = grid.data();
  for(unsigned i = 0, n = grid.size(); i < n; i++) {
    float diff = a[i]-b[i];
    if(fabs(diff) > 0.0001) {
      cerr << "Values differ " << a[i] << " != " << b[i] << " at index " << i << "\n";
      exit(1);
    }
  }
}


void MolGridder::set_center(gfloat3 c) {
  center_set = true;
  center = c;
}

//set grid parameters from an example dx
void MolGridder::set_from_example(const string& examplefile) {
  CartesianMGrid cgrid = read_dx<float>(examplefile);
  center = cgrid.center();
  resolution = cgrid.resolution();
  dimension = resolution * (cgrid.grid().dimension(0) - 1);
  center_set = true;
  N = 1 + round(dimension / resolution);
}

//read in specified grid files and set usergrids and grid parameters appropriately
void MolGridder::set_usergrids(const vector<string>& userfiles) {

  if (userfiles.size() > 0) {
    //use grid specified by user
    if (random_rotate || random_translate) {
      cerr << "Random rotation/translation is not supported with user grids.\n";
      exit(1);
    }
  }
  usergrids.resize(0);
  usergrids.reserve(userfiles.size());
  center_set = true; //defiend by user grids
  for (unsigned i = 0, ng = userfiles.size(); i < ng; i++) {
    //load dx file
    CartesianMGrid cgrid = read_dx<float>(userfiles[i]);
    unsigned npts = cgrid.grid().dimension(0);
    usergrids.push_back(cgrid.grid());
    //check resolution/dimensions
    if (i == 0) {
      resolution = cgrid.resolution();
      center = cgrid.center();
      dimension = resolution * (npts - 1); //fencepost
    } else {
      if (cgrid.resolution() != resolution) {
        cerr << "Inconsistent resolutions in grids: " << cgrid.resolution()
            << " vs "
            << resolution << "\n";
        exit(1);
      }
      double dim = resolution * (npts - 1);
      if (dim != dimension) {
        cerr << "Inconsistent dimensions in grids: " << dim << " vs "
            << dimension << "\n";
        exit(1);
      }

      if (center != cgrid.center()) {
        cerr << "Inconsistent centers in grids\n";
        exit(1);
      }
    }
  }
  N = 1 + round(dimension / resolution);
}


bool MolGridder::readMolecule(bool timeit) {

  model m;
  if (!mols.readMoleculeIntoModel(m)) return false;
  setLigand(m);

  boost::timer::cpu_timer t;
  setGrid(gpu);

  if (timeit) {
    cout << "Grid Time: " << t.elapsed().wall << "\n";
    cpuSetGridCheck();
  }
  return true;
}

//return true if grid only contains zeroes
static bool gridIsEmpty(const Grid3f& grid) {
  for (const float *ptr = grid.data(), *end = grid.data() + grid.size();
      ptr != end; ptr++) {
    if (*ptr != 0.0) return false;
  }
  return true;
}

//output map for each grid
void MolGridder::outputMAP(const std::string& base) {

  for(unsigned a = 0, na = usergrids.size(); a < na; a++) {
    string fname = base + "_usergrid_" + boost::lexical_cast<string>(a) + ".dx";
    ofstream out(fname.c_str());
    libmolgrid::write_dx<float>(out,grid[a], center, resolution);
  }

  unsigned roff = usergrids.size();
  std::vector<std::string> recnames = rectyper->get_type_names();
  for (unsigned a = 0, na = recnames.size(); a < na; a++) {
    //this is for debugging, so avoid outputting empty grids
    if (!gridIsEmpty(grid[roff+a])) {
      string name = recnames[a];
      string fname = base + "_rec_" + name + ".map";
      ofstream out(fname.c_str());
      libmolgrid::write_map<float>(out, grid[roff+a], center, resolution);
    }
  }
  std::vector<std::string> lignames = ligtyper->get_type_names();
  roff += recnames.size();
  for (unsigned a = 0, na = lignames.size(); a < na; a++) {
    if (!gridIsEmpty(grid[roff+a])) {
      string name = lignames[a];
      string fname = base + "_lig_" + name + ".map";
      ofstream out(fname.c_str());
      libmolgrid::write_map<float>(out,grid[roff+a], center, resolution);
    }
  }
}

//output an dx map for each grid
void MolGridder::outputDX(const std::string& base) {

  for(unsigned a = 0, na = usergrids.size(); a < na; a++) {
    string fname = base + "_lig_" + boost::lexical_cast<string>(a) + ".dx";
    ofstream out(fname.c_str());
    libmolgrid::write_dx<float>(out,grid[a], center, resolution);
  }

  std::vector<std::string> recnames = rectyper->get_type_names();
  unsigned roff = usergrids.size();
  for (unsigned a = 0, na = recnames.size(); a < na; a++) {
    //this is for debugging, so avoid outputting empty grids
    if (!gridIsEmpty(grid[roff+a])) {
      string name = recnames[a];
      string fname = base + "_rec_" + name + ".dx";
      ofstream out(fname.c_str());
      libmolgrid::write_dx<float>(out, grid[roff+a], center, resolution);
    }
  }

  std::vector<std::string> lignames = ligtyper->get_type_names();
  roff += recnames.size();
  for (unsigned a = 0, na = lignames.size(); a < na; a++) {
    if (!gridIsEmpty(grid[roff+a])) {
      string name = lignames[a];
      string fname = base + "_lig_" + name + ".dx";
      ofstream out(fname.c_str());
      libmolgrid::write_dx<float>(out,grid[roff+a], center, resolution);
    }
  }
}

//output binary form of raw data in 3D multi-channel form
void MolGridder::outputBIN(const std::string& base, bool outputrec, bool outputlig) {
  unsigned chan = 0;
  if (outputrec) chan += rectyper->num_types() + usergrids.size();
  if (outputlig) chan += ligtyper->num_types();
  string outname = base + "." + boost::lexical_cast<string>(N) + "." + boost::lexical_cast<string>(chan)+".binmap";

  ofstream binout(outname.c_str());
  if(!binout) {
    throw file_error(outname, false);
  }
  unsigned roff = 0;
  if(outputrec) {
    for(unsigned i = 0, n = usergrids.size(); i < n; i++) {
      write_bin(binout, grid[i]);
    }
    roff += usergrids.size();
    for (unsigned a = 0, na = rectyper->num_types(); a < na; a++) {
      write_bin(binout, grid[roff+a]);
    }
  }
  roff += rectyper->num_types();
  if(outputlig) {
    for (unsigned a = 0, na = ligtyper->num_types(); a < na; a++) {
      write_bin(binout, grid[roff+a]);
    }
  }
}
