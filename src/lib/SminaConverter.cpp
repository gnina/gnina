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
//sets up data structures used by both text and binary
//we link with smina to nsure compatibility
unsigned convertParsing(OBMol& mol, parsing_struct& p, context& c)
{
	mol.AddHydrogens();

	mol.SetAutomaticFormalCharge(false);
	DeleteHydrogens(mol); //leaves just polars

	vector<vector<int> > rigid_fragments; //the vector of all the rigid molecule fragments, using atom indexes
	map<unsigned int, obbranch> tree;

	//we kind of assume a connected molecule
	unsigned best_root_atom = FindFragments(mol, rigid_fragments);
	unsigned torsdof=rigid_fragments.size()-1;

	unsigned int root_piece = 0;
    for (unsigned j = 0; j < rigid_fragments.size(); j++)
    {
      if (IsIn((rigid_fragments[j]), best_root_atom)) {root_piece=j; break;} //this is the root rigid molecule fragment
    }

    ConstructTree(tree, rigid_fragments, root_piece, mol, true);

    OutputTree(mol, c, p, tree, torsdof);

  /*  stringstream str;
    str << "TORSDOF  ";
    str << torsdof;
    add_pdbqt_context(c, str.str()); */
    return torsdof;
}

template <class T>
static void convert(OBMol& mol, T& serialout, ostream& out)
{
	parsing_struct p;
	context c;
	unsigned torsdof = convertParsing(mol, p, c);
	serialout << torsdof;
	serialout << p;
	serialout << c;
}

//text output
void convertText(OBMol& mol, ostream& out)
{
	boost::archive::text_oarchive serialout(out,boost::archive::no_header|boost::archive::no_tracking);
	convert(mol, serialout, out);
}

//binary output
void convertBinary(OBMol& mol,  ostream& out)
{
	//by definition, smina format is gzipped
	boost::iostreams::filtering_stream<boost::iostreams::output> strm;
	strm.push(boost::iostreams::gzip_compressor());
	strm.push(out);

	boost::archive::binary_oarchive serialout(strm,boost::archive::no_header|boost::archive::no_tracking);
	convert(mol, serialout, strm);
}

} //namespace SminaConverter
