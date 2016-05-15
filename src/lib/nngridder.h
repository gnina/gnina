#ifndef NNGRIDDER_H
#define NNGRIDDER_H

#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/quaternion.hpp>
#include <cuda.h>
#include <vector_types.h>
#include "atom_type.h"
#include "box.h"
#include "gridoptions.h"
#include "molgetter.h"

using namespace std;

/* Maintains atom grid information.  Stores a model of the receptor/ligand with
 * MolGetter, but also numerical grids for every protein/ligand atom type.
 */
class NNGridder
{
public:
  typedef boost::math::quaternion<double> quaternion;
protected:
	grid_dims dims; //this is a cube
	quaternion Q;
	vec trans;
	double resolution;
	double dimension;
	double radiusmultiple; //extra to consider past vdw radius
	double randtranslate;
	bool binary; //produce binary occupancies
	bool randrotate;
	bool gpu; //use gpu

	vector<boost::multi_array<float, 3> > receptorGrids;
	vector<boost::multi_array<float, 3> > ligandGrids;
	vector<int> rmap; //map atom types to position in grid vectors
	vector<int> lmap;

	vector<vec> recCoords; //these don't change
	vector<float> recRadii; //nor do these, note for complex type maps the same
	vector<short> recWhichGrid; // the atom type based grid index

	vector<float> ligRadii;
	vector<short> ligWhichGrid; //only change if ligand changes


	//gpu data structures, these all point to device mem
	float *gpu_receptorGrids;
	float *gpu_ligandGrids;

	float3 *gpu_receptorCoords;
	float *gpu_recRadii;
	short *gpu_recWhichGrid;

	float3 *gpu_ligandCoords;
	float *gpu_ligRadii;
	short *gpu_ligWhichGrid;

	void setRecGPU();
	void setLigGPU();

	pair<unsigned, unsigned> getrange(const grid_dim& dim, double c, double r);

	//return the density value for atom a at the provided point
	float calcPoint(const vec& coords, double ar, const vec& pt);

	//set the relevant grid points for a
	//return false if atom not in grid
	bool setAtom(const vec& coords, double radius, boost::multi_array<float, 3>& grid);

	//set the relevant grid points for passed info
	void setAtoms(const vector<vec>& coords, const vector<short>& gridindex, const vector<float>& radii, vector<boost::multi_array<float, 3> >& grids);

	//GPU accelerated version
	void setAtomsGPU(unsigned natoms, float3 *coords, short *gridindex, float *radii, unsigned ngrids, float *grids);

	//output a grid the file in map format (for debug)
	void outputMAPGrid(ostream& out, boost::multi_array<float, 3>& grid);

	//return a string representation of the atom type(s) represented by index
	//in map - this isn't particularly efficient, but is only for debug purposes
	string getIndexName(const vector<int>& map, unsigned index) const;

	//setup ligand/receptor maps
	//setup grid dimensions and zero-init
	void setMapsAndGrids(const gridoptions& opt);

	//set the center of the grid, must reset receptor/ligand
	void setCenter(double x, double y, double z);

	static void zeroGrids(vector<boost::multi_array<float, 3> >& grid);
	static void cudaCopyGrids(vector<boost::multi_array<float, 3> >& grid, float* gpu_grid);

	//for debugging
	static bool compareGrids(boost::multi_array<float, 3>& g1, boost::multi_array<float, 3>& g2, const char *name, int index);


public:

	NNGridder(): resolution(0.5), dimension(24), radiusmultiple(1.5),
			randtranslate(0), binary(false), randrotate(false), gpu(false),
			gpu_receptorGrids(NULL), gpu_ligandGrids(NULL),
			gpu_receptorCoords(NULL), gpu_recRadii(NULL), gpu_recWhichGrid(NULL),
			gpu_ligandCoords(NULL), gpu_ligRadii(NULL), gpu_ligWhichGrid(NULL)
			{}

	void initialize(const gridoptions& opt);

	//set grids (receptor and ligand)
	//reinits should be set to true if have different molecule than previously scene
	void setModel(const model& m, bool reinitlig=false, bool reinitrec=false);

	//return string detailing the configuration (size.channels)
	string getParamString(bool outputrec, bool outputlig) const;

	//output an AD4 map for each grid
	void outputMAP(const string& base);

	//output binary form of raw data in 3D multi-channel form (types are last)
	void outputBIN(ostream& out, bool outputrec = true, bool outputlig = true);

	//set vector to full set of grids
	void outputMem(vector<float>& out);

	//initialize default receptor/ligand maps
	static void createDefaultRecMap(vector<int>& map);
	static void createDefaultLigMap(vector<int>& map);

	unsigned nchannels() const { return receptorGrids.size() + ligandGrids.size(); }

	//for debugging, run non-gpu code and compre to values in current grids
	bool cpuSetModelCheck(const model& m, bool reinitlig=false, bool reinitrec=false);
};

/* This gridder uses a MolGetter to read molecules */
class NNMolsGridder : public NNGridder
{
public:
  typedef boost::math::quaternion<double> quaternion;
private:
	MolGetter mols; //this stores the models

public:

	NNMolsGridder(const gridoptions& opt);

	//read a molecule (return false if unsuccessful)
	//set the ligand grid appropriately
	bool readMolecule(bool timeit);

};


#define CUDA_CHECK(condition) \
  /* Code block avoids redefinition of cudaError_t error */ \
  do { \
    cudaError_t error = condition; \
    if(error != cudaSuccess) { cerr << " " << cudaGetErrorString(error) << ": " << __FILE__ << ":" << __LINE__ << "\n"; exit(1); } \
  } while (0)


#endif
