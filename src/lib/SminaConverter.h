/*
 * SminaOBMolConverter.h
 *
 *  Created on: Jun 4, 2014
 *      Author: dkoes
 */

#ifndef SMINACONVERTER_H_
#define SMINACONVERTER_H_

#include <openbabel/mol.h>
#include <iostream>
#include <vector>
#include "parsing.h"
#include "PDBQTUtilities.h"

/* Routines for converting a molecule to smina format.
 */
namespace SminaConverter
{
	//can optoinal specify a desired root atom and atoms to that should not have rotatable bonds
	//text output
	void convertText(OpenBabel::OBMol& mol, std::ostream& out, int rootatom, const std::vector<int>& norotate);
	void convertText(OpenBabel::OBMol& mol, std::ostream& out);
	//binary output
	void convertBinary(OpenBabel::OBMol& mol, std::ostream& out, int rootatom, const std::vector<int>& norotate);
	void convertBinary(OpenBabel::OBMol& mol, std::ostream& out);

	//convert obmol to smina parsing struct and context; return numtors
	unsigned convertParsing(OpenBabel::OBMol& mol, parsing_struct& p, context& c, int rootatom, const std::vector<int>& norotate, bool addH=true);
	unsigned convertParsing(OpenBabel::OBMol& mol, parsing_struct& p, context& c, bool addH=true);

	//class for efficiently converting multi-conformer molecule
	class MCMolConverter
	{
		OpenBabel::OBMol mol;
		std::vector<std::vector<int> > rigid_fragments; //the vector of all the rigid molecule fragments, using atom indexes
		std::map<unsigned int, obbranch> tree;
		unsigned torsdof;
	public:
		MCMolConverter(OpenBabel::OBMol& m);
		//output smina data for specified conformer
		void convertConformer(unsigned c, std::ostream& out);
	};
};

#endif /* SMINACONVERTER_H_ */
