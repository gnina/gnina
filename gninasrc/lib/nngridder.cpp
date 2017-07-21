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
#include <boost/timer/timer.hpp>
#include <cuda_runtime.h>
#include "gpu_math.h"


#define CUDA_CHECK(condition) \
  /* Code block avoids redefinition of cudaError_t error */ \
  do { \
    cudaError_t error = condition; \
    if(error != cudaSuccess) { cerr << " " << cudaGetErrorString(error) << ": " << __FILE__ << ":" << __LINE__ << "\n"; exit(1); } \
  } while (0)

using namespace boost;



//output a grid the file in map format (for debug)
void NNGridder::outputMAPGrid(ostream& out, Grid& grid)
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

//output a grid the file in dx format (for debug)
void NNGridder::outputDXGrid(ostream& out, Grid& grid)
{
  unsigned n = dims[0].n+1;
  out.precision(5);
  setprecision(5);
  out << fixed;
  out << "object 1 class gridpositions counts " << n << " " << n << " " << " " << n << "\n";
  out << "origin";
  for (unsigned i = 0; i < 3; i++)
  {
    out << " " << dims[i].begin;
  }
  out << "\n";
  out << "delta " << resolution << " 0 0\ndelta 0 " << resolution << " 0\ndelta 0 0 " << resolution << "\n";
  out << "object 2 class gridconnections counts " << n << " " << n << " " << " " << n << "\n";
  out << "object 3 class array type double rank 0 items [ " << n*n*n << "] data follows\n";
  //now coordinates - x,y,z
  unsigned total = 0;
  for (unsigned i = 0; i < n; i++)
  {
    for (unsigned j = 0; j < n; j++)
    {
      for (unsigned k = 0; k < n; k++)
      {
        out << grid[i][j][k];
        total++;
        if(total % 3 == 0) out << "\n";
        else out << " ";
      }
    }
  }
}

//read a dx file into grid
bool NNGridder::readDXGrid(istream& in, vec& center, double& res, Grid& grid)
{
	string line;
	vector<string> tokens;

	res = 0;
	getline(in, line);
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	if(tokens.size() != 8) return false;
	unsigned n = lexical_cast<unsigned>(tokens[7]);
	if(lexical_cast<unsigned>(tokens[6]) != n) return false;
	if(lexical_cast<unsigned>(tokens[5]) != n) return false;

	//the center
	getline(in, line);
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	if(tokens.size() != 4) return false;
	double x = lexical_cast<double>(tokens[1]);
	double y = lexical_cast<double>(tokens[2]);
	double z = lexical_cast<double>(tokens[3]);

	//the transformation matrix, which has the resolution
	getline(in, line);
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	if(tokens.size() != 4) return false;
	res = lexical_cast<float>(tokens[1]);

	getline(in, line);
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	if(tokens.size() != 4) return false;
	if(res != lexical_cast<float>(tokens[2])) return false;

	getline(in, line);
	split(tokens, line, is_any_of(" \t"), token_compress_on);
	if(tokens.size() != 4) return false;
	if(res != lexical_cast<float>(tokens[3])) return false;

	//figure out center
	double half = res*n/2.0;
	center[0] = x+half;
	center[1] = y+half;
	center[2] = z+half;

	//grid connections
	getline(in, line);
	//object 3
	getline(in, line);

	//data begins
	grid.resize(extents[n][n][n]);
	fill_n(grid.data(), grid.num_elements(), 0.0);

	unsigned total = 0;
	for (unsigned i = 0; i < n; i++)
	{
		for (unsigned j = 0; j < n; j++)
		{
			for (unsigned k = 0; k < n; k++)
			{
				in >> grid[i][j][k];
				total++;
			}
		}
	}
	if(total != n*n*n) return false;

	return true;
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

	if (ret.str().length() > 32) //there are limits on file name lengths
		return altret.str();
	else
		return ret.str();
}



//return string detailing the configuration (size.channels)
string NNGridder::getParamString(bool outputrec, bool outputlig) const
		{
	unsigned n = dims[0].n + 1;
	unsigned chan = 0;
	if (outputrec)
		chan += receptorGrids.size()+userGrids.size();
	if (outputlig)
		chan += ligandGrids.size();
	return lexical_cast<string>(n) + "." + lexical_cast<string>(chan);
}

//return true if grid only contains zeroes
static bool gridIsEmpty(const NNGridder::Grid& grid)
{
	for (const float *ptr = grid.data(), *end = grid.data()
			+ grid.num_elements(); ptr != end; ptr++)
	{
		if (*ptr != 0.0)
			return false;
	}
	return true;
}

//output an AD4 map for each grid
void NNGridder::outputMAP(const string& base)
{
	for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
	{
		//this is for debugging, so avoid outputting empty grids
		if (!gridIsEmpty(receptorGrids[a]))
		{
			string name = getIndexName(rmap, a);
			string fname = base + "_rec_" + name + ".map";
			ofstream out(fname.c_str());
			outputMAPGrid(out, receptorGrids[a]);
		}
	}
	for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
	{
		if (!gridIsEmpty(ligandGrids[a]))
		{
			string name = getIndexName(lmap, a);
			string fname = base + "_lig_" + name + ".map";
			ofstream out(fname.c_str());
			outputMAPGrid(out, ligandGrids[a]);
		}
	}

}

//output an AD4 map for each grid
void NNGridder::outputDX(const string& base)
{
  for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
  {
    //this is for debugging, so avoid outputting empty grids
    if (!gridIsEmpty(receptorGrids[a]))
    {
      string name = getIndexName(rmap, a);
      string fname = base + "_rec_" + name + ".dx";
      ofstream out(fname.c_str());
      outputDXGrid(out, receptorGrids[a]);
    }
  }
  for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
  {
    if (!gridIsEmpty(ligandGrids[a]))
    {
      string name = getIndexName(lmap, a);
      string fname = base + "_lig_" + name + ".dx";
      ofstream out(fname.c_str());
      outputDXGrid(out, ligandGrids[a]);
    }
  }

}

//output binary form of raw data in 3D multi-channel form (types are last)
void NNGridder::outputBIN(ostream& out, bool outputrec, bool outputlig)
{
	unsigned n = dims[0].n + 1;

	if(outputrec)
	{
		//receptor
		for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
		{
			assert(receptorGrids[a].shape()[0] == n);

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

		//user grids
		for (unsigned a = 0, na = userGrids.size(); a < na; a++)
		{
			assert(userGrids[a].shape()[0] == n);
			for (unsigned i = 0; i < n; i++)
			{
				for (unsigned j = 0; j < n; j++)
				{
					for (unsigned k = 0; k < n; k++)
					{
						out.write((char*) &userGrids[a][i][j][k],
								sizeof(float));
					}
				}
			}
		}
	}

	if(outputlig)
	{
		//ligand
		for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
		{
			assert(ligandGrids[a].shape()[0] == n);
			for (unsigned i = 0; i < n; i++)
			{
				for (unsigned j = 0; j < n; j++)
				{
					for (unsigned k = 0; k < n; k++)
					{
						out.write((char*) &ligandGrids[a][i][j][k],
								sizeof(float));
					}
				}
			}
		}
	}

}

void NNGridder::outputMem(vector<float>& out)
{
	unsigned n = dims[0].n + 1;
	unsigned gsize = n * n * n;
	out.resize(gsize * receptorGrids.size() + gsize * ligandGrids.size());

	float *ptr = &out[0];
	for (unsigned a = 0, na = receptorGrids.size(); a < na; a++)
	{
		memcpy(ptr, receptorGrids[a].origin(), gsize * sizeof(float));
		ptr += gsize;
	}
	for (unsigned a = 0, na = ligandGrids.size(); a < na; a++)
	{
		memcpy(ptr, ligandGrids[a].origin(), gsize * sizeof(float));
		ptr += gsize;
	}
}

//copy gpu grid to passed cpu grid
//TODO: rearchitect to use flat grids on the cpu?
void NNGridder::cudaCopyGrids(vector<Grid>& grid, float* gpu_grid)
{
	for(unsigned i = 0, n = grid.size(); i < n; i++)
	{
		float *cpu = grid[i].data();
		unsigned sz = grid[i].num_elements();
		CUDA_CHECK(cudaMemcpyAsync(cpu, gpu_grid, sz*sizeof(float),cudaMemcpyDeviceToHost));
		gpu_grid += sz;
	}
}

bool NNGridder::compareGrids(Grid& g1, Grid& g2, const char *name, int index)
{
	if(g1.size() != g2.size())
	{
		cerr << "Invalid initialize grid size " << name << index << "\n";
		return false;
	}
	for(unsigned i = 0, I = g1.size(); i < I; i++)
	{
		if(g1[i].size() != g2[i].size())
		{
			cerr << "Invalid secondary grid size " << name << index << "\n";
			return false;
		}
		for(unsigned j = 0, J = g1[i].size(); j < J; j++)
		{
			if(g1[i][j].size() != g2[i][j].size())
			{
				cerr << "Invalid tertiary grid size " << name << index << "\n";
				return false;
			}
			for(unsigned k = 0, K = g1[i][j].size(); k < K; k++)
			{
				float diff = g1[i][j][k] - g2[i][j][k];
				if(fabs(diff) > 0.0001)
				{
					cerr << "Values differ " << g1[i][j][k] <<  " != " << g2[i][j][k] << " " << name << index << " " << i <<","<<j<<","<<k <<"\n";
					return false;
				}
			}
		}
	}
	return true;
}



void NNGridder::setCenter(double x, double y, double z)
{
	if(userGrids.size() > 0)
	{
		cerr << "Attempting to re-set grid center with user specified grids.\n";
		abort();
	}
	gmaker.setCenter(x,y,z);

	trans = vec(x, y, z);
	int numpts = round(dimension / resolution);
	double half = dimension / 2.0;
	dims[0].begin = x - half;
	dims[0].end = x + half;
	dims[0].n = numpts;

	dims[1].begin = y - half;
	dims[1].end = y + half;
	dims[1].n = numpts;

	dims[2].begin = z - half;
	dims[2].end = z + half;
	dims[2].n = numpts;
}

void NNGridder::setMapsAndGrids(const gridoptions& opt)
{
	if (opt.recmap.size() == 0)
		GridMaker::createDefaultRecMap(rmap);
	else
		GridMaker::createAtomTypeMap(opt.recmap, rmap);

	if (opt.ligmap.size() == 0)
		GridMaker::createDefaultLigMap(lmap);
	else
		GridMaker::createAtomTypeMap(opt.ligmap, lmap);

	//setup grids,
	dimension = opt.dim;
	resolution = opt.res;

	int numpts = round(dimension / resolution);
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
			if (receptorGrids[i].num_elements() == 0)
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
			if (ligandGrids[i].num_elements() == 0)
			{
				ligandGrids[i].resize(extents[n][n][n]);
				fill_n(ligandGrids[i].data(), ligandGrids[i].num_elements(), 0.0);
			}
		}
	}

	//check for empty mappings
	for (unsigned i = 0, nr = receptorGrids.size(); i < nr; i++)
	{
		if (receptorGrids[i].num_elements() == 0)
		{
			cerr << "Empty slot in receptor types: " << i
					<< ", possible duplicate?\n";
			receptorGrids[i].resize(extents[n][n][n]);
			fill_n(receptorGrids[i].data(), receptorGrids[i].num_elements(), 0.0);
		}
	}
	for (unsigned i = 0, nl = ligandGrids.size(); i < nl; i++)
	{
		if (ligandGrids[i].num_elements() == 0)
		{
			cerr << "Empty slot in ligand types: " << i
					<< ", possible duplicate?\n";
			ligandGrids[i].resize(extents[n][n][n]);
			fill_n(ligandGrids[i].data(), ligandGrids[i].num_elements(), 0.0);
		}
	}

	vec savedcenter;
	userGrids.resize(0); userGrids.reserve(opt.usergrids.size());
	for(unsigned i = 0, ng = opt.usergrids.size(); i < ng; i++)
	{
		//load dx file
		ifstream gridfile(opt.usergrids[i].c_str());
		if(!gridfile)
		{
			cerr << "Could not open grid file " << opt.usergrids[i] << "\n";
			abort();
		}
		vec gridcenter;
		double gridres = 0;
		Grid g;
		if(!readDXGrid(gridfile,gridcenter, gridres, g))
		{
			cerr << "I couldn't understand the provided dx file " << opt.usergrids[i] << ". I apologize for not being more informative and possibly being too picky about my file formats.\n";
			abort();
		}

		unsigned npts = g.shape()[0];

		//check resolution/dimensions
		if(i == 0)
		{
			resolution = gridres;
			savedcenter = gridcenter;
			dimension = resolution*(npts-1); //fencepost
			setCenter(savedcenter[0],savedcenter[1],savedcenter[2]);
		}
		else
		{
			if(gridres != resolution)
			{
				cerr << "Inconsistent resolutions in grids: " << gridres << " vs " << resolution << "\n";
				abort();
			}
			double dim = resolution*(npts-1);
			if(dim != dimension)
			{
				cerr << "Inconsistent dimensions in grids: " << dim << " vs " << dimension << "\n";
				abort();
			}

			for(unsigned c = 0; c < 3; c++)
			{
				if(savedcenter[c] != gridcenter[c])
				{
					cerr << "Inconsistent center in grids at position " << c << ": " << savedcenter[c] << " vs " << gridcenter[c] << "\n";
					abort();
				}
			}
		}

		//add grid
		userGrids.push_back(g);

	}

}



NNMolsGridder::NNMolsGridder(const gridoptions& opt)
{
	initialize(opt);
	//open receptor
	tee log(true);
	FlexInfo finfo(log); //dummy
	mols.create_init_model(opt.receptorfile, "", finfo, log);

	//set ligand file
	mols.setInputFile(opt.ligandfile);
}

void NNGridder::initialize(const gridoptions& opt)
{
	binary = opt.binary;
	resolution = opt.res;
	radiusmultiple = 1.5;
	randtranslate = opt.randtranslate;
	randrotate = opt.randrotate;
	gpu = opt.gpu;
	Q = quaternion(0, 0, 0, 0);

	gmaker.initialize(resolution, opt.dim, radiusmultiple, binary);

	if (binary)
		radiusmultiple = 1.0;

	setMapsAndGrids(opt);

	if(gpu)
	{
		//allocate gpu memory for grids
		unsigned nrgrids = receptorGrids.size();
		unsigned nlgrids = ligandGrids.size();
		assert(nrgrids > 0);
		assert(nlgrids > 0);
		unsigned n = receptorGrids[0].num_elements();
		CUDA_CHECK(cudaMalloc(&gpu_receptorGrids, nrgrids*n*sizeof(float)));
		CUDA_CHECK(cudaMalloc(&gpu_ligandGrids, nlgrids*n*sizeof(float)));
	}
}

//allocate (if neccessary) and copy recCoords,radii, and whichgrid
//assumes cpu version is set
void NNGridder::setRecGPU()
{
	if(gpu_receptorAInfo == NULL) {
		assert(sizeof(vec) == sizeof(float3));
		CUDA_CHECK(cudaMalloc(&gpu_receptorAInfo, recAInfo.size()*sizeof(float4)));
		CUDA_CHECK(cudaMemcpy(gpu_receptorAInfo, &recAInfo[0], recAInfo.size()*sizeof(float4),cudaMemcpyHostToDevice));
	}

	if(gpu_recWhichGrid == NULL) {
		CUDA_CHECK(cudaMalloc(&gpu_recWhichGrid, recWhichGrid.size()*sizeof(float3)));
		CUDA_CHECK(cudaMemcpy(gpu_recWhichGrid, &recWhichGrid[0], recWhichGrid.size()*sizeof(float3),cudaMemcpyHostToDevice));
	}
}

void NNGridder::setLigGPU()
{
	if(gpu_ligandAInfo == NULL) {
		CUDA_CHECK(cudaMalloc(&gpu_ligandAInfo, ligWhichGrid.size()*sizeof(float4)));
		//ligand coordinates
	}

	if(gpu_ligWhichGrid == NULL) {
		CUDA_CHECK(cudaMalloc(&gpu_ligWhichGrid, ligWhichGrid.size()*sizeof(short)));
		CUDA_CHECK(cudaMemcpy(gpu_ligWhichGrid, &ligWhichGrid[0], ligWhichGrid.size()*sizeof(short),cudaMemcpyHostToDevice));
	}
}

void NNGridder::setModel(const model& m, bool reinitlig, bool reinitrec)
{
	//compute center from ligand
	const atomv& atoms = m.get_movable_atoms();
	assert(atoms.size() == m.coordinates().size());

	if(userGrids.size() > 0)
	{
		//use grid specified by user
		if(randrotate || randtranslate)
		{
			cerr << "Random rotation/translation is not supported with user grids.\n";
			abort();
		}
		//the grid is already set - do nothing
	}
	else //center around ligand
	{
		vec center(0,0,0);
		for (unsigned i = 0, n = atoms.size(); i < n; i++)
		{
			center += m.coordinates()[i];
		}
		center /= atoms.size();

		//apply random modifications
		if (randrotate)
		{
			double d = rand() / double(RAND_MAX);
			double r1 = rand() / double(RAND_MAX);
			double r2 = rand() / double(RAND_MAX);
			double r3 = rand() / double(RAND_MAX);
			Q = NNGridder::quaternion(1, r1 / d, r2 / d, r3 / d);
		}

		if (randtranslate)
		{
			double offx = rand() / double(RAND_MAX / 2.0) - 1.0;
			double offy = rand() / double(RAND_MAX / 2.0) - 1.0;
			double offz = rand() / double(RAND_MAX / 2.0) - 1.0;
			center[0] += offx * randtranslate;
			center[1] += offy * randtranslate;
			center[2] += offz * randtranslate;
		}

		setCenter(center[0], center[1], center[2]);
	}

	//set constant arrays if needed
	if(recAInfo.size() == 0 || reinitrec)
	{
		const atomv& atoms = m.get_fixed_atoms();
		recAInfo.resize(0); recAInfo.reserve(atoms.size());
		recWhichGrid.resize(0); recWhichGrid.reserve(atoms.size());

		for (unsigned i = 0, n = atoms.size(); i < n; i++)
		{
			const atom& a = atoms[i];
			if(rmap[a.sm] >= 0) {
				float4 ai = {a.coords[0], a.coords[1], a.coords[2], xs_radius(a.sm)};
				recWhichGrid.push_back(rmap[a.sm]);
				recAInfo.push_back(ai);
			}
		}

		if(gpu) setRecGPU();

	}


	if(ligRadii.size() == 0 || reinitlig)
	{
		const atomv& atoms = m.get_movable_atoms();
		assert(atoms.size() == m.coordinates().size());
		ligRadii.resize(atoms.size());
		ligWhichGrid.resize(atoms.size());
		//we can't omit stupid atoms since they are included in the coordinates
		for (unsigned i = 0, n = atoms.size(); i < n; i++)
		{
			atom a = atoms[i];
			ligWhichGrid[i] = lmap[a.sm];
			ligRadii[i] = xs_radius(a.sm);
		}

		if(gpu) setLigGPU();
	}

	const vecv& coords = m.coordinates();
	vector<float4> ainfo; ainfo.reserve(coords.size());
	for(unsigned i = 0, na = coords.size(); i < na; i++)
	{
		float4 ai = {coords[i][0],coords[i][1],coords[i][2],ligRadii[i]};
		ainfo.push_back(ai);
	}

	if(gpu)
	{
		unsigned nlatoms = m.coordinates().size();
		CUDA_CHECK(cudaMemcpy(gpu_ligandAInfo, &ainfo[0], nlatoms*sizeof(float4),cudaMemcpyHostToDevice));
		//cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024*1024*4);

		gmaker.setAtomsGPU<float>(recAInfo.size(),gpu_receptorAInfo, gpu_recWhichGrid, Q, receptorGrids.size(), gpu_receptorGrids);
		cudaCopyGrids(receptorGrids, gpu_receptorGrids);

		gmaker.setAtomsGPU<float>(nlatoms, gpu_ligandAInfo, gpu_ligWhichGrid, Q, ligandGrids.size(), gpu_ligandGrids);
		cudaCopyGrids(ligandGrids, gpu_ligandGrids);

		CUDA_CHECK(cudaDeviceSynchronize());
	}
	else
	{
		//put in to gridmaker format

		gmaker.setAtomsCPU(recAInfo, recWhichGrid, Q, receptorGrids);
		gmaker.setAtomsCPU(ainfo, ligWhichGrid,  Q, ligandGrids);
	}
}

bool NNGridder::cpuSetModelCheck(const model& m, bool reinitlig, bool reinitrec)
{
	vector<Grid> savedRGrid = receptorGrids;
	vector<Grid> savedLGrid = ligandGrids;

	bool savedgpu = gpu;
	gpu = false;
	bool saverandrotate = randrotate;
	randrotate = false; //keep current Q, don't have a convenient way to chekc rand translate
	setModel(m, reinitlig, reinitrec);
	randrotate = saverandrotate;

	gpu = savedgpu;

	for(unsigned i = 0, n = receptorGrids.size(); i < n; i++)
	{
		compareGrids(receptorGrids[i], savedRGrid[i], "receptor", i);
	}
	for(unsigned i = 0, n = ligandGrids.size(); i < n; i++)
	{
		compareGrids(ligandGrids[i], savedLGrid[i], "ligand", i);
	}
}


//read a molecule (return false if unsuccessful)
//set the ligand grid appropriately
bool NNMolsGridder::readMolecule(bool timeit)
{
	model m;
	if (!mols.readMoleculeIntoModel(m))
		return false;

	timer::cpu_timer t;
	setModel(m, true);

	if(timeit)
	{
		cout << "Grid Time: " << t.elapsed().wall << "\n";
		//DEBUG CODE BELOW
		if(gpu) cpuSetModelCheck(m, true);

	}
	return true;
}
