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

#include "atom_type.h"
#include "box.h"
#include "molgetter.h"
#include "options.h"

using namespace std;

/* Maintains atom grid information.  Stores a model of the receptor/ligand with
 * MolGetter, but also numerical grids for every protein/ligand atom type.
 */
class NNGridder
{
	MolGetter mols; //this stores the models
	grid_dims dims; //this is  cuge
	double resolution;
	double radiusmultiple; //extra to consider past vdw radius
	bool binary; //produce binary occupancies

	vector<boost::multi_array<float, 3> > receptorGrids;
	vector<boost::multi_array<float, 3> > ligandGrids;
	vector<int> rmap; //map atom types to position in grid vectors
	vector<int> lmap;

	pair<unsigned, unsigned> getrange(const grid_dim& dim, double c, double r);

	//return the density value for atom a at the provided point
	float calcPoint(const atom& a, const vec& pt);

	//set the relevant grid points for a
	void setAtom(const atom& a, boost::multi_array<float, 3>& grid);

	//output a grid the file in map format (for debug)
	void outputMAPGrid(ostream& out, boost::multi_array<float, 3>& grid);

	//return a string representation of the atom type(s) represented by index
	//in map - this isn't particularly efficient, but is only for debug purposes
	string getIndexName(const vector<int>& map, unsigned index) const;

public:
	NNGridder(const cmdoptions& opt, const vector<int>& recmap, const vector<int>& ligmap);

	//read a molecule (return false if unsuccessful)
	//set the ligand grid appropriately
	bool readMolecule();

	//return string detailing the configuration (size.channels)
	string getParamString() const;

	//output an AD4 map for each grid
	void outputMAP(const string& base);

	//output binary form of raw data in 3D multi-channel form (types are last)
	void outputBIN(ostream& out);
};

#endif
