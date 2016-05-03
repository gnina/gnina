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
	//return false if atom not in grid
	bool setAtom(const atom& a, boost::multi_array<float, 3>& grid);

	//output a grid the file in map format (for debug)
	void outputMAPGrid(ostream& out, boost::multi_array<float, 3>& grid);

	//return a string representation of the atom type(s) represented by index
	//in map - this isn't particularly efficient, but is only for debug purposes
	string getIndexName(const vector<int>& map, unsigned index) const;

	//setup ligand/receptor maps
	//setup grid dimensions and zero-init
	void setMapsAndGrids(const gridoptions& opt);

	static void zeroGrids(vector<boost::multi_array<float, 3> >& grid);
public:

	NNGridder(): resolution(0.5), radiusmultiple(1.5), binary(false) {}

	//return string detailing the configuration (size.channels)
	string getParamString(bool outputrec, bool outputlig) const;

	//output an AD4 map for each grid
	void outputMAP(const string& base);

	//output binary form of raw data in 3D multi-channel form (types are last)
	void outputBIN(ostream& out, bool outputrec = true, bool outputlig = true);

	//initialize default receptor/ligand maps
	static void createDefaultRecMap(vector<int>& map);
	static void createDefaultLigMap(vector<int>& map);

	unsigned nchannels() const { return receptorGrids.size() + ligandGrids.size(); }
};

/* This gridder uses a MolGetter to read molecules */
class NNMolsGridder : public NNGridder
{
public:
  typedef boost::math::quaternion<double> quaternion;
private:
	MolGetter mols; //this stores the models

public:

	NNMolsGridder(const gridoptions& opt, quaternion q = quaternion(1,0,0,0));

	//read a molecule (return false if unsuccessful)
	//set the ligand grid appropriately
	bool readMolecule();

};

/* This gridder extracts ligand from model rather than reading from file */
class NNModelGridder : public NNGridder
{
public:
  typedef boost::math::quaternion<double> quaternion;
private:

public:

	NNModelGridder() {}

	void initialize(const gridoptions& opt);

	void setReceptor(const model& m);
	void setLigand(const model& m);

};

#endif
