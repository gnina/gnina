/*
 * SminaConverter.cpp
 *
 *  Created on: Jun 4, 2014
 *      Author: dkoes
 *
 *  Convert internal molecular data (ie OBMol) into smina parse tree.
 */

#include "SminaConverter.h"
#include "parsing.h"
#include "PDBQTUtilities.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <map>

namespace SminaConverter {

using namespace OpenBabel;
using namespace std;


MCMolConverter::MCMolConverter(OpenBabel::OBMol& m): mol(m)
{
	//precompute fragments and tree
	int nc = mol.NumConformers();
	mol.AddHydrogens();

	mol.SetAutomaticFormalCharge(false);
	DeleteHydrogens(mol); //leaves just polars

	//we kind of assume a connected molecule
	unsigned best_root_atom = FindFragments(mol, rigid_fragments);
	torsdof=rigid_fragments.size()-1;

	unsigned int root_piece = 0;
    for (unsigned j = 0; j < rigid_fragments.size(); j++)
    {
      if (IsIn((rigid_fragments[j]), best_root_atom)) {root_piece=j; break;} //this is the root rigid molecule fragment
    }

    ConstructTree(tree, rigid_fragments, root_piece, mol, true);

    assert(nc == mol.NumConformers()); //didn't lose any in analysis, did we?
}


//output data for this conformer
void MCMolConverter::convertConformer(unsigned conf, std::ostream& out)
{
	parsing_struct p;
	context c;

	mol.SetConformer(conf);
	std::map<unsigned int, obbranch> tmptree(tree); //tree gets modified by outputtree
    OutputTree(mol, c, p, tmptree, torsdof);

	boost::iostreams::filtering_stream<boost::iostreams::output> strm;
	strm.push(boost::iostreams::gzip_compressor());
	strm.push(out);

	boost::archive::binary_oarchive serialout(strm,boost::archive::no_header|boost::archive::no_tracking);

	serialout << torsdof;
	serialout << p;
	serialout << c;
}

//sets up data structures used by both text and binary
//we link with smina to ensure compatibility
//rootatom, an obatom index (starting at 1) can be specified, if not
//the "best" root is chosen
unsigned convertParsing(OBMol& mol, parsing_struct& p, context& c, int rootatom /* = 0*/)
{
	mol.AddHydrogens();

	mol.SetAutomaticFormalCharge(false);
	DeleteHydrogens(mol); //leaves just polars

	vector<vector<int> > rigid_fragments; //the vector of all the rigid molecule fragments, using atom indexes
	map<unsigned int, obbranch> tree;

	//we kind of assume a connected molecule
	unsigned best_root_atom = FindFragments(mol, rigid_fragments);
	unsigned torsdof=rigid_fragments.size()-1;

	if(rootatom > 0)
	{
		//user user supplied root
		best_root_atom = rootatom;
	}

	unsigned int root_piece = 0;
    for (unsigned j = 0; j < rigid_fragments.size(); j++)
    {
      if (IsIn((rigid_fragments[j]), best_root_atom)) {root_piece=j; break;} //this is the root rigid molecule fragment
    }

    ConstructTree(tree, rigid_fragments, root_piece, mol, true);

    OutputTree(mol, c, p, tree, torsdof);

    return torsdof;
}

template <class T>
static void convert(OBMol& mol, T& serialout, ostream& out, int rootatom=0)
{
	parsing_struct p;
	context c;
	unsigned torsdof = convertParsing(mol, p, c, rootatom);
	serialout << torsdof;
	serialout << p;
	serialout << c;
}

//text output
void convertText(OBMol& mol, ostream& out, int rootatom/*=0*/)
{
	boost::archive::text_oarchive serialout(out,boost::archive::no_header|boost::archive::no_tracking);
	convert(mol, serialout, out, rootatom);
}

//binary output
void convertBinary(OBMol& mol,  ostream& out, int rootatom/*=0*/)
{
	//by definition, smina format is gzipped
	boost::iostreams::filtering_stream<boost::iostreams::output> strm;
	strm.push(boost::iostreams::gzip_compressor());
	strm.push(out);

	boost::archive::binary_oarchive serialout(strm,boost::archive::no_header|boost::archive::no_tracking);
	convert(mol, serialout, strm, rootatom);
}

} //namespace SminaConverter
