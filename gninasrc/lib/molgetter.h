/*
 * molgetter.h
 *
 *  Created on: Jun 5, 2014
 *      Author: dkoes
 */

#ifndef MOLGETTER_H_
#define MOLGETTER_H_

#include "covinfo.h"
#include "flexinfo.h"
#include "model.h"
#include "obmolopener.h"

// this class abstracts reading molecules from a file
// we have three means of input:
// openbabel for general molecular data (default)
// vina parse_pdbqt for pdbqt files (one ligand, obey rotational bonds)
// smina format
class MolGetter {
  model initm;
  tee *log;
  CovInfo cinfo;
  enum Type { OB, PDBQT, SMINA, GNINA, NONE }; // different inputs

  Type type;
  path lpath;
  bool add_hydrogens;   // add hydrogens before calculating atom types
  bool strip_hydrogens; // strip them after (more efficient)
  // openbabel data structs
  OpenBabel::OBConversion conv;
  obmol_opener infileopener;

  // smina data structs
  izfile infile;

  // pdbqt data
  bool pdbqtdone;

  // covalent data
  OpenBabel::OBMol covres;              // covalently bonding residue on receptor
  OpenBabel::OBAtom *covatom = nullptr; // covalently bonding atom within this residue
  vec covpos;                           // position for covalently bonding ligand atom
  bool covres_isflex = false; //true if covalently bonded residue should be flexible

  OpenBabel::OBMol covmol;                  // current ligand being docked
  OpenBabel::OBMol origcovmol; //original ligand before covalent modifications
  std::vector<std::vector<int>> match_list; // smarts matches
  unsigned matchpos = UINT_MAX;             // current position in match_list

public:
  MolGetter(bool addH = true, bool stripH = true)
      : add_hydrogens(addH), strip_hydrogens(stripH), type(NONE), pdbqtdone(false) {}

  MolGetter(const std::string &rigid_name, const std::string &flex_name, FlexInfo &finfo, CovInfo &ci, bool addH,
            bool stripH, tee &l)
      : log(&l), cinfo(ci), add_hydrogens(addH), strip_hydrogens(stripH), type(NONE), pdbqtdone(false) {
    create_init_model(rigid_name, flex_name, finfo, l);
  }

  // create the initial model from the specified receptor files
  void create_init_model(const std::string &rigid_name, const std::string &flex_name, FlexInfo &finfo, tee &log);

  // setup for reading from fname
  void setInputFile(const std::string &fname);

  // initialize model to initm and add next molecule
  // return false if no molecule available;
  bool readMoleculeIntoModel(model &m);

  // return model without ligand
  const model &getInitModel() const { return initm; }

private:
  bool createCovalentMoleculeInModel(model &m);
};

#endif /* MOLGETTER_H_ */
