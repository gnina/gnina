/*
 * GninaConverter.cpp
 *
 *  Created on: Jun 4, 2014
 *      Author: dkoes
 *
 *  Convert internal molecular data (ie OBMol) into gnina parse tree.
 */

#include "GninaConverter.h"
#include "parsing.h"
#include "PDBQTUtilities.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <map>
#include <openbabel/obconversion.h>

namespace GninaConverter {

using namespace OpenBabel;
using namespace std;

MCMolConverter::MCMolConverter(OpenBabel::OBMol& m)
    : mol(m) {
  //precompute fragments and tree
  int nc = mol.NumConformers();
  mol.SetConformer(0);
  mol.AddHydrogens();

  mol.SetAutomaticFormalCharge(false);
  DeleteHydrogens(mol); //leaves just polars
  if (mol.NumAtoms() == 0) return;

  //we kind of assume a connected molecule
  unsigned best_root_atom = FindFragments(mol, rigid_fragments);
  torsdof = rigid_fragments.size() - 1;

  unsigned int root_piece = 0;
  for (unsigned j = 0; j < rigid_fragments.size(); j++) {
    if (IsIn((rigid_fragments[j]), best_root_atom)) {
      root_piece = j;
      break;
    } //this is the root rigid molecule fragment
  }

  ConstructTree(tree, rigid_fragments, root_piece, mol, true);

  if (nc != mol.NumConformers())  //didn't lose any in analysis, did we?
      {
    abort(); //there was a bug in openbabel where addhydrogens would eliminate conformers
  }
}

//output data for this conformer
void MCMolConverter::convertConformer(unsigned conf, std::ostream& out) {
  parsing_struct p;
  context c;

  mol.SetConformer(conf);

  std::map<unsigned int, obbranch> tmptree(tree); //tree gets modified by outputtree
  OutputTree(mol, c, p, tmptree, torsdof);

  boost::iostreams::filtering_stream<boost::iostreams::output> strm;
  strm.push(boost::iostreams::gzip_compressor());
  strm.push(out);

  boost::archive::binary_oarchive serialout(strm,
      boost::archive::no_header | boost::archive::no_tracking);

  serialout << torsdof;
  serialout << p;
  serialout << c;
}

//sets up data structures used by both text and binary
//we link with gnina to ensure compatibility
//rootatom, an obatom index (starting at 1) can be specified, if not
//the "best" root is chosen
unsigned convertParsing(OBMol& mol, parsing_struct& p, context& c, int rootatom,
    const vector<int>& norotate, bool addH) {
  if (addH) mol.AddHydrogens();
  mol.PerceiveBondOrders();
  mol.SetAromaticPerceived();
  mol.SetAutomaticFormalCharge(false);
  DeleteHydrogens(mol); //leaves just polars

  vector<vector<int> > rigid_fragments; //the vector of all the rigid molecule fragments, using atom indexes
  map<unsigned int, obbranch> tree;

  //we kind of assume a connected molecule
  unsigned best_root_atom = FindFragments(mol, rigid_fragments, rootatom,
      norotate);
  unsigned torsdof = rigid_fragments.size() - 1;

  if (rootatom > 0) {
    //user user supplied root
    best_root_atom = rootatom;
  }

  unsigned int root_piece = 0;
  for (unsigned j = 0; j < rigid_fragments.size(); j++) {
    if (IsIn((rigid_fragments[j]), best_root_atom)) {
      root_piece = j;
      break;
    } //this is the root rigid molecule fragment
  }

  ConstructTree(tree, rigid_fragments, root_piece, mol, true);

  OutputTree(mol, c, p, tree, torsdof);

  return torsdof;
}

unsigned convertParsing(OpenBabel::OBMol& mol, parsing_struct& p, context& c,
    bool addH) {
  std::vector<int> nr;
  return convertParsing(mol, p, c, 0, nr, addH);
}

template<class T>
static void convert(OBMol& mol, T& serialout, ostream& out, int rootatom,
    const vector<int>& norotate) {
  parsing_struct p;
  context c;
  unsigned torsdof = convertParsing(mol, p, c, rootatom, norotate);
  serialout << torsdof;
  serialout << p;
  serialout << c;
}

//text output
void convertText(OBMol& mol, ostream& out, int rootatom,
    const vector<int>& norotate) {
  boost::archive::text_oarchive serialout(out,
      boost::archive::no_header | boost::archive::no_tracking);
  convert(mol, serialout, out, rootatom, norotate);
}

void convertText(OpenBabel::OBMol& mol, std::ostream& out) {
  std::vector<int> nr;
  convertText(mol, out, 0, nr);
}

//binary output
void convertBinary(OBMol& mol, ostream& out, int rootatom,
    const vector<int>& norotate) {
  //by definition, gnina format is gzipped
  boost::iostreams::filtering_stream<boost::iostreams::output> strm;
  strm.push(boost::iostreams::gzip_compressor());
  strm.push(out);

  boost::archive::binary_oarchive serialout(strm,
      boost::archive::no_header | boost::archive::no_tracking);
  convert(mol, serialout, strm, rootatom, norotate);
}

void convertBinary(OpenBabel::OBMol& mol, std::ostream& out) {
  std::vector<int> nr;
  convertBinary(mol, out, 0, nr);
}

} //namespace GninaConverter
