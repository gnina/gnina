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
#include "parsing.h"

/* Routines for converting a molecule to smina format.
 */
namespace SminaConverter
{
	//text output
	void convertText(OpenBabel::OBMol& mol, std::ostream& out);
	//binary output
	void convertBinary(OpenBabel::OBMol& mol, std::ostream& out);

	//convert obmol to smina parsing struct and context; return numtors
	unsigned convertParsing(OpenBabel::OBMol& mol, parsing_struct& p, context& c);

};

#endif /* SMINACONVERTER_H_ */
