/*
 * molgetter.h
 *
 *  Created on: Jun 5, 2014
 *      Author: dkoes
 */

#ifndef MOLGETTER_H_
#define MOLGETTER_H_

#include "model.h"
#include "obmolopener.h"


//this class abstracts reading molecules from a file
//we have three means of input:
//openbabel for general molecular data (default)
//vina parse_pdbqt for pdbqt files (one ligand, obey rotational bonds)
//smina format
class MolGetter
{
	const model& initm;
	enum Type {OB, PDBQT, SMINA}; //different inputs

	Type type;
	path lpath;
	bool add_hydrogens;
	//openbabel data structs
	OpenBabel::OBConversion conv;
	obmol_opener infileopener;

	//smina data structs
	izfile infile;

	//pdbqt data
	bool pdbqtdone;
public:
	MolGetter(const model& m, bool addH): initm(m), add_hydrogens(addH), pdbqtdone(false) {}

	//setup for reading from fname
	void setInputFile(const std::string& fname);

	//initialize model to initm and add next molecule
	//return false if no molecule available;
	bool readMoleculeIntoModel(model &m);
};



#endif /* MOLGETTER_H_ */
