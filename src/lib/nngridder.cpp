/*
 * gninagrid.cpp
 *  GPL v2 (due to OpenBabel)/ BSD License (if you remove OpenBabel dependency).
 *  Created on: Nov 6, 2015
 *      Author: dkoes
 *
 * Class for converting smina model into grids of atom type occupancies.
 */

#include "nngridder.h"
#include <cmath>

using namespace boost;

static void createDefaultMap(const char *names[], vector<int>& map)
{
	map.assign(smina_atom_type::NumTypes, -1);
	const char **nameptr = names;
	unsigned cnt = 0;
	while(*nameptr != NULL) {
		string name(*nameptr);
		//note that if we every start using merged atom types by default
		//this code will have to be updated
		smt t = string_to_smina_type(name);
		if(t < smina_atom_type::NumTypes) //valid
		{
			map[t] = cnt;
			cnt++;
		}
		else //should never happen
		{
			cerr << "Invalid atom type " << name << "\n";
			exit(-1);
		}
		nameptr++;
	}
}

//initialize default receptor/ligand maps
//these were determined by evaluating how common various atoms types are
void NNGridder::createDefaultRecMap(vector<int>& map)
{
	const char *names[] = {"AliphaticCarbonXSHydrophobe",
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
			"Zinc", NULL};

	createDefaultMap(names, map);
}

void NNGridder::createDefaultLigMap(vector<int>& map) {
	const char *names[] = {"AliphaticCarbonXSHydrophobe",
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
			NULL};
	createDefaultMap(names, map);
}


//return the occupancy for atom a at point x,y,z
float NNGridder::calcPoint(const atom& a, const vec& pt)
{
	double rsq = (pt - a.coords).norm_sqr();
	double ar = xs_radius(a.sm);
	if(binary)
	{
		//is point within radius?
		if(rsq < ar*ar)
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
		if(dist >= ar*1.5)
		{
			return 0.0;
		}
		else if(dist <= ar)
		{
			//return gaussian
			double h = 0.5*ar;
			double ex = -dist*dist/(2*h*h);
			return exp(ex);
		}
		else //return quadratic
		{
			double h = 0.5*ar;
			double eval = 1.0/(M_E*M_E); //e^(-2)
			double q = dist*dist*eval/(h*h) - 6.0 *eval*dist/h + 9.0 * eval;
			return q;
		}
	}
	return 0.0;
}


//return the range of grid points spanned from c-r to c+r within dim
pair<unsigned, unsigned> NNGridder::getrange(const grid_dim& dim, double c, double r)
{
	pair<unsigned, unsigned> ret(0,0);
	double low = c-r-dim.begin;
	if(low > 0)
	{
		ret.first = floor(low/resolution);
	}

	double high = c+r-dim.begin;
	if(high > 0) //otherwise zero
	{
		ret.second = min(dim.n,(sz)ceil(high/resolution));
	}
	return ret;
}


//set the relevant grid points for a
//figure out what volume of the grid is relevant for this atom and for each
//grid point in this volume, convert it into world coordinates and call calcpoint
//to get its value
bool NNGridder::setAtom(const atom& origa, boost::multi_array<float, 3>& grid)
{
  atom a = origa;
  if(Q.real() != 0) { //apply rotation
    vec tpt = a.coords - trans;
    quaternion p(0,tpt[0], tpt[1], tpt[2]);
    p = Q*p*(conj(Q)/norm(Q));

    vec newpt(p.R_component_2(),p.R_component_3(),p.R_component_4());
    newpt += trans;

    a.coords = newpt;
  }

	double r = xs_radius(a.sm)*radiusmultiple;
	vector< pair<unsigned,unsigned> > ranges;
	for(unsigned i = 0; i < 3; i++)
	{
		ranges.push_back(getrange(dims[i], a.coords[i], r));
	}
	bool isset = false; //is this atom actaully within grid?
	//for every grid point possibly overlapped by this atom
	for(unsigned i = ranges[0].first, iend = ranges[0].second; i < iend; i++)
	{
		for(unsigned j = ranges[1].first, jend = ranges[1].second; j < jend; j++)
		{
			for(unsigned k = ranges[2].first, kend = ranges[2].second; k < kend; k++)
			{
				double x = dims[0].begin+i*resolution;
				double y = dims[1].begin+j*resolution;
				double z= dims[2].begin+k*resolution;
				double val = calcPoint(a, vec(x,y,z));

				if(binary)
				{
					if(val != 0)
						grid[i][j][k] = 1.0; //don't add, just 1 or 0
				}
				else
					grid[i][j][k] += val;

				isset = true;
			}
		}
	}
	return isset;
}


//output a grid the file in map format (for debug)
void NNGridder::outputMAPGrid(ostream& out, boost::multi_array<float, 3>& grid)
{
	unsigned max = dims[0].n + 1;
	out.precision(5);
	out << "GRID_PARAMETER_FILE\nGRID_DATA_FILE\nMACROMOLECULE\n";
	out << "SPACING " << resolution << "\n";
	out << "NELEMENTS " << max - 1 << " " << max - 1 << " " << max - 1 << "\n";
	out << "CENTER";
	for (unsigned i = 0; i < 3; i++)
	{
		double c = (dims[i].end + dims[i].begin) / 2.0;
		out << " " << c;
	}
	out << "\n";

	//now coordinates - z,y,x
	for (unsigned k = 0; k < max; k++)
	{
		for (unsigned j = 0; j < max; j++)
		{
			for (unsigned i = 0; i < max; i++)
			{
				out << grid[i][j][k] << "\n";
			}
		}
	}
}

//return a string representation of the atom type(s) represented by index
//in map - this isn't particularly efficient, but is only for debug purposes
string NNGridder::getIndexName(const vector<int>& map, unsigned index) const
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

	if(ret.str().length() > 32) //there are limits on file name lengths
		return altret.str();
	else
		return ret.str();
}

//create a mapping from atom type ids to a unique id given a file specifying
//what types we care about (anything missing is ignored); if multiple types are
//on the same line, they are merged, if the file isn't specified, use default mapping
//return total number of types
//map is indexed by smina_atom_type, maps to -1 if type should be ignored
static int createAtomTypeMap(const string& fname, vector<int>& map)
{
	map.assign(smina_atom_type::NumTypes, -1);

	if(fname.size() == 0)
	{
		//default mapping
		cerr << "Map file not specified\n";
		exit(-1);
	}
	else
	{
		int cnt = 0;
		ifstream in(fname.c_str());

		if(!in)
		{
			cerr << "Could not open " << fname << "\n";
			exit(-1);
		}
		string line;
		while(getline(in, line))
		{
			vector<string> types;
			split(types, line, is_any_of("\t \n"));
			for(unsigned i = 0, n = types.size(); i < n; i++)
			{
				const string& name = types[i];
				smt t = string_to_smina_type(name);
				if(t < smina_atom_type::NumTypes) //valid
				{
					map[t] = cnt;
				}
				else if(name.size() > 0) //this ignores consecutive delimiters
				{
					cerr << "Invalid atom type " << name << "\n";
					exit(-1);
				}
			}
			if(types.size() > 0)
				cnt++;
		}
		return cnt;
	}
}


//read a molecule (return false if unsuccessful)
//set the ligand grid appropriately
bool NNMolsGridder::readMolecule()
{
	model m;
	if (!mols.readMoleculeIntoModel(m))
		return false;

	//clear ligand array
	for (unsigned i = 0, n = ligandGrids.size(); i < n; i++)
	{
		fill_n(ligandGrids[i].data(), ligandGrids[i].num_elements(), 0.0);
	}

	//fill in heavy atoms
	const atomv& atoms = m.get_movable_atoms();
	assert(atoms.size() == m.coordinates().size());
	for (unsigned i = 0, n = atoms.size(); i < n; i++)
	{
		atom a = atoms[i];
		a.coords = m.coordinates()[i]; //have to explicitly set coords
		int pos = lmap[a.sm];
		assert(pos < ligandGrids.size());
		if (pos >= 0)
			setAtom(a, ligandGrids[pos]);
	}
	return true;
}

//return string detailing the configuration (size.channels)
string NNGridder::getParamString(bool outputrec, bool outputlig) const
{
	unsigned n = dims[0].n + 1;
	unsigned chan = 0;
	if(outputrec)
		chan += receptorGrids.size();
	if(outputlig)
		chan += ligandGrids.size();
	return lexical_cast<string>(n) + "." + lexical_cast<string>(chan);
}

//return true if grid only contains zeroes
static bool gridIsEmpty(const multi_array<float, 3>& grid)
{
	for(const float *ptr = grid.data(), *end = grid.data() + grid.num_elements(); ptr != end; ptr++)
	{
		if(*ptr != 0.0) return false;
	}
	return true;
}

//output an AD4 map for each grid
void NNGridder::outputMAP(const string& base)
{
	for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
	{
		//this is for debugging, so avoid outputting empty grids
		if(!gridIsEmpty(receptorGrids[a]))
		{
			string name = getIndexName(rmap, a);
			string fname = base + "_rec_" + name + ".map";
			ofstream out(fname.c_str());
			outputMAPGrid(out, receptorGrids[a]);
		}
	}
	for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
	{
		if(!gridIsEmpty(ligandGrids[a]))
		{
			string name = getIndexName(lmap, a);
			string fname = base + "_lig_" + name + ".map";
			ofstream out(fname.c_str());
			outputMAPGrid(out, ligandGrids[a]);
		}
	}

}

//output binary form of raw data in 3D multi-channel form (types are last)
void NNGridder::outputBIN(ostream& out, bool outputrec, bool outputlig)
{
	unsigned n = dims[0].n + 1;
  if(outputrec)
  {
    for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
    {
      for (unsigned i = 0; i < n; i++)
      {
        for (unsigned j = 0; j < n; j++)
        {
          for (unsigned k = 0; k < n; k++)
          {
            //when you see this many loops you known you're going to generate a lot of data..

                out.write((char*) &receptorGrids[a][i][j][k],
                    sizeof(float));
            }
          }
        }
    }
	}
  if(outputlig)
  {
    for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
    {
      for (unsigned i = 0; i < n; i++)
      {
        for (unsigned j = 0; j < n; j++)
        {
          for (unsigned k = 0; k < n; k++)
          {
            out.write((char*) &ligandGrids[a][i][j][k], sizeof(float));
          }
        }
      }
    }
  }
}

void NNGridder::outputMem(vector<float>& out)
{
	unsigned n = dims[0].n + 1;
	unsigned gsize = n*n*n;
	out.resize(gsize*receptorGrids.size()+gsize*ligandGrids.size());

	float *ptr = &out[0];
  for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
  {
  	memcpy(ptr, receptorGrids[a].origin(), gsize*sizeof(float));
  	ptr += gsize;
  }
  for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
  {
  	memcpy(ptr, ligandGrids[a].origin(), gsize*sizeof(float));
  	ptr += gsize;
  }
}


//set all the elements of a grid to zero
void NNGridder::zeroGrids(vector<boost::multi_array<float, 3> >& grid)
{
	for(unsigned i = 0, n = grid.size(); i < n; i++)
	{
		boost::multi_array<float, 3>& g = grid[i];
		std::fill(g.data(),g.data()+g.num_elements(), 0.0);
	}
}


void NNGridder::setMapsAndGrids(const gridoptions& opt)
{
	if(opt.recmap.size() == 0)
		NNGridder::createDefaultRecMap(rmap);
	else
		createAtomTypeMap(opt.recmap,rmap);

	if(opt.ligmap.size() == 0)
		NNGridder::createDefaultLigMap(lmap);
	else
		createAtomTypeMap(opt.ligmap,lmap);

	//setup grids,
	trans = vec(opt.x,opt.y,opt.z);
	int numpts = round(opt.dim / opt.res);
	double half = opt.dim / 2.0;
	dims[0].begin = opt.x - half;
	dims[0].end = opt.x + half;
	dims[0].n = numpts;

	dims[1].begin = opt.y - half;
	dims[1].end = opt.y + half;
	dims[1].n = numpts;

	dims[2].begin = opt.z - half;
	dims[2].end = opt.z + half;
	dims[2].n = numpts;
	unsigned n = numpts + 1; //fencepost

	receptorGrids.reserve(smina_atom_type::NumTypes);
	ligandGrids.reserve(smina_atom_type::NumTypes);

	for (unsigned at = 0; at < smina_atom_type::NumTypes; at++)
	{
		if (rmap[at] >= 0) //valid type for receptor
		{
			unsigned i = rmap[at];
			if (receptorGrids.size() <= i)
				receptorGrids.resize(i + 1);
			if(receptorGrids[i].num_elements() == 0)
			{
				receptorGrids[i].resize(extents[n][n][n]);
				fill_n(receptorGrids[i].data(), receptorGrids[i].num_elements(), 0.0);
			}
		}

		if (lmap[at] >= 0)
		{
			unsigned i = lmap[at];
			if (ligandGrids.size() <= i)
				ligandGrids.resize(i + 1);
			if(ligandGrids[i].num_elements() == 0)
			{
				ligandGrids[i].resize(extents[n][n][n]);
				fill_n(ligandGrids[i].data(), ligandGrids[i].num_elements(), 0.0);
			}
		}
	}

	//check for empty mappings
	for(unsigned i = 0, nr = receptorGrids.size(); i < nr; i++)
	{
	  if(receptorGrids[i].num_elements() == 0)
	  {
	    cerr << "Empty slot in receptor types: " << i << ", possible duplicate?\n";
	    receptorGrids[i].resize(extents[n][n][n]);
      fill_n(receptorGrids[i].data(), receptorGrids[i].num_elements(), 0.0);
	  }
	}
  for(unsigned i = 0, nl = ligandGrids.size(); i < nl; i++)
  {
    if(ligandGrids[i].num_elements() == 0)
    {
      cerr << "Empty slot in ligand types: " << i << ", possible duplicate?\n";
      ligandGrids[i].resize(extents[n][n][n]);
      fill_n(ligandGrids[i].data(), ligandGrids[i].num_elements(), 0.0);
    }
  }

}

NNMolsGridder::NNMolsGridder(const gridoptions& opt, quaternion q)
{
	binary = opt.binary;
	resolution = opt.res;
	radiusmultiple = 1.5;
	Q = q;

	if(binary) radiusmultiple = 1.0;

	setMapsAndGrids(opt);

	//open receptor
	tee log(true);
	FlexInfo finfo(log); //dummy
	mols.create_init_model(opt.receptorfile, "", finfo, log);

	//initialize receptor
	const model& m = mols.getInitModel();

	const atomv& atoms = m.get_fixed_atoms();
	for (unsigned i = 0, n = atoms.size(); i < n; i++)
	{
		const atom& a = atoms[i];
		int pos = rmap[a.sm];
		if (pos >= 0)
			setAtom(a, receptorGrids[pos]);
	}

	//set ligand file
	mols.setInputFile(opt.ligandfile);
}


void NNModelGridder::initialize(const gridoptions& opt)
{
	binary = opt.binary;
	resolution = opt.res;
	radiusmultiple = 1.5;
	Q = quaternion(1,0,0,0);

	if(binary) radiusmultiple = 1.0;

	setMapsAndGrids(opt);

	//setRecptor needs to be called to init receptor grids
}

void NNModelGridder::setReceptor(const model& m)
{
	zeroGrids(receptorGrids);
	const atomv& atoms = m.get_fixed_atoms();
	for(unsigned i = 0, n = atoms.size(); i < n; i++) {
		const atom& a = atoms[i];
		int pos = rmap[a.sm];
		if(pos >= 0)
			setAtom(a, receptorGrids[pos]);
	}
}

void NNModelGridder::setLigand(const model& m)
{
	zeroGrids(ligandGrids);
	//fill in heavy atoms
	const atomv& atoms = m.get_movable_atoms();
	assert(atoms.size() == m.coordinates().size());
	for (unsigned i = 0, n = atoms.size(); i < n; i++)
	{
		atom a = atoms[i];
		a.coords = m.coordinates()[i]; //have to explicitly set coords
		int pos = lmap[a.sm];
		assert(pos < ligandGrids.size());
		if (pos >= 0)
			setAtom(a, ligandGrids[pos]);
	}
}
