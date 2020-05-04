/*
 * molgetter.h
 *
 *  Created on: Jun 5, 2014
 *      Author: dkoes
 */
#include "molgetter.h"
#include "parse_pdbqt.h"
#include "parsing.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/timer/timer.hpp>
#include "GninaConverter.h"

//create the initial model from the specified receptor files
//mostly because Matt kept complaining about it, this will automatically create
//pdbqts if necessary using open babel
void MolGetter::create_init_model(const std::string& rigid_name,
    const std::string& flex_name, FlexInfo& finfo, tee& log) {
  if (rigid_name.size() > 0) {
    //support specifying flexible residues explicitly as pdbqt, but only
    //in compatibility mode where receptor is pdbqt as well
    if (flex_name.size() > 0) {
      if (finfo.hasContent()) {
        throw usage_error(
            "Cannot mix -flex option with -flexres or -flexdist options.");
      }
      if (boost::filesystem::extension(rigid_name) != ".pdbqt") {
        throw usage_error("Cannot use -flex option with non-PDBQT receptor.");
      }
      ifile rigidin(rigid_name);
      ifile flexin(flex_name);
      initm = parse_receptor_pdbqt(rigid_name, rigidin, flex_name, flexin);
    } else
      if (!finfo.hasContent()
          && boost::filesystem::extension(rigid_name) == ".pdbqt") {
        //compatibility mode - read pdbqt directly with no openbabel shenanigans
        ifile rigidin(rigid_name);
        initm = parse_receptor_pdbqt(rigid_name, rigidin);
      } else {
        //default, openbabel mode
        using namespace OpenBabel;
        obmol_opener fileopener;
        OBConversion conv;
        conv.SetOutFormat("PDBQT");
        conv.AddOption("r", OBConversion::OUTOPTIONS); //rigid molecule, otherwise really slow and useless analysis is triggered
        conv.AddOption("c", OBConversion::OUTOPTIONS); //single combined molecule
        fileopener.openForInput(conv, rigid_name);
        OBMol rec;
        if (!conv.Read(&rec)) throw file_error(rigid_name, true);

        rec.AddHydrogens(true);
        FOR_ATOMS_OF_MOL(a, rec){
          a->GetPartialCharge();
        }
        OBMol rigid;
        std::string flexstr;

        if(rec.NumResidues() > 0){
          try{
            // Can fail with std::runtime_error if `--flex_limit` is set
            finfo.extractFlex(rec, rigid, flexstr);
          }
          catch(std::runtime_error &e){
            // --flex_limit exceeded; print error and quit
            log << e.what() << "\n";
            std::exit(-1);
          }
        }
        else{
          // No information about residues, whole receptor treated as rigid
          rigid = rec;
        }

        std::string recstr = conv.WriteString(&rigid);
        std::stringstream recstream(recstr);
        
        if (flexstr.size() > 0) //have flexible component
        {
          std::stringstream flexstream(flexstr);
          initm = parse_receptor_pdbqt(rigid_name, recstream, flex_name,
              flexstream);
        } else { //rigid only
          initm = parse_receptor_pdbqt(rigid_name, recstream);
        }

      }

  }

  if (strip_hydrogens) initm.strip_hydrogens();
}

//setup for reading from fname
void MolGetter::setInputFile(const std::string& fname) {
  if (fname.size() > 0) //zer if no_lig
      {
    lpath = path(fname);
    if (lpath.extension() == ".pdbqt") {
      //built-in pdbqt parsing that respects rotabable bonds in pdbqt
      type = PDBQT;
      pdbqtdone = false;
    } else
      if (infile.open(lpath, ".smina", true)) //smina always gzipped
          {
        type = SMINA;
      } else
        if (infile.open(lpath, ".gnina", true)) //gnina always gzipped
            {
          type = GNINA;
        } else
          if (fname.length() > 0) //openbabel
              {
            type = OB;
            //clear in case we had previous file
            infileopener.clear();
            infileopener.openForInput(conv, fname);
            VINA_CHECK(conv.SetOutFormat("PDBQT"));

          }
  }
}

//initialize model to initm and add next molecule
//return false if no molecule available;
bool MolGetter::readMoleculeIntoModel(model &m) {
  //reinit the model
  m = initm;
  switch (type) {
  case SMINA:
  case GNINA: {
    parsing_struct p;
    context c;
    unsigned torsdof = 0;
    if (!infile) return false;
    try {
      boost::archive::binary_iarchive serialin(infile,
          boost::archive::no_header | boost::archive::no_tracking);
      serialin >> torsdof;
      serialin >> p;
      serialin >> c;

      non_rigid_parsed nr;
      postprocess_ligand(nr, p, c, torsdof);
      VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());

      pdbqt_initializer tmp;
      tmp.initialize_from_nrp(nr, c, true);
      tmp.initialize(nr.mobility_matrix());

      if (c.sdftext.valid()) {
        //set name
        m.set_name(c.sdftext.name);
      }

      if (strip_hydrogens) tmp.m.strip_hydrogens();
      m.append(tmp.m);

      return true;
    } catch (boost::archive::archive_exception& e) {
      return false;
    }
  }
    break;
  case PDBQT: {
    if (pdbqtdone) return false; //can only read one
    model lig = parse_ligand_pdbqt(lpath);
    if (strip_hydrogens) lig.strip_hydrogens();
    m.append(lig);
    pdbqtdone = true;
    return true;
  }
    break;
  case OB: {
    OpenBabel::OBMol mol;
    while (conv.Read(&mol)) //will return after first success
    {
      std::string name = mol.GetTitle();
      mol.StripSalts();
      m.set_name(name);
      try {
        parsing_struct p;
        context c;
        unsigned torsdof = GninaConverter::convertParsing(mol, p, c,
            add_hydrogens);
        non_rigid_parsed nr;
        postprocess_ligand(nr, p, c, torsdof);
        VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());

        pdbqt_initializer tmp;
        tmp.initialize_from_nrp(nr, c, true);
        tmp.initialize(nr.mobility_matrix());
        if (strip_hydrogens) tmp.m.strip_hydrogens();

        m.append(tmp.m);
        return true;
      } catch (parse_error& e) {
        std::cerr << "\n\nParse error with molecule " << mol.GetTitle()
            << " in file \"" << e.file.string() << "\": " << e.reason << '\n';
        continue;
      }
    }

    return false; //no valid molecules read
  }
  case NONE:
    return true; //nolig
    break;
  }
#ifndef __NVCC__
  return false; //shouldn't get here
#endif
}
