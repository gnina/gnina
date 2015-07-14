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
#include "SminaConverter.h"

//setup for reading from fname
void MolGetter::setInputFile(const std::string& fname)
{
	if (fname.size() > 0) //zer if no_lig
	{
		lpath = path(fname);
		if (lpath.extension() == ".pdbqt")
		{
			//built-in pdbqt parsing that respects rotabable bonds in pdbqt
			type = PDBQT;
			pdbqtdone = false;
		}
		else if (infile.open(lpath, ".smina", true)) //smina always gzipped
		{
			type = SMINA;
		}
		else //openbabel
		{
			type = OB;
			//clear in case we had previous vile
			infileopener.clear();
			infileopener.openForInput(conv, fname);
			VINA_CHECK(conv.SetOutFormat("PDBQT"));
		}
	}
}

//initialize model to initm and add next molecule
//return false if no molecule available;
bool MolGetter::readMoleculeIntoModel(model &m)
{
	//reinit the model
	m = initm;
	switch (type)
	{
	case SMINA:
		{
		parsing_struct p;
		context c;
		unsigned torsdof = 0;
		if (!infile)
			return false;
		try
		{
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

			if(c.sdftext.valid())
			{
				//set name
				m.set_name(c.sdftext.name);
			}

			m.append(tmp.m);
			return true;
		}
		catch (boost::archive::archive_exception& e)
		{
			return false;
		}
	}
		break;
	case PDBQT:
		{
		if (pdbqtdone)
			return false; //can only read one
		m.append(parse_ligand_pdbqt(lpath));
		pdbqtdone = true;
		return true;
	}
		break;
	case OB:
		{
		OpenBabel::OBMol mol;
		while (conv.Read(&mol)) //will return after first success
		{
			std::string name = mol.GetTitle();
			mol.StripSalts();
			m.set_name(name);
			try
			{
				parsing_struct p;
				context c;
				unsigned torsdof = SminaConverter::convertParsing(mol, p, c, add_hydrogens);
				non_rigid_parsed nr;
				postprocess_ligand(nr, p, c, torsdof);
				VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());

				pdbqt_initializer tmp;
				tmp.initialize_from_nrp(nr, c, true);
				tmp.initialize(nr.mobility_matrix());

				m.append(tmp.m);
				return true;
			}
			catch (parse_error& e)
			{
				std::cerr << "\n\nParse error with molecule "
						<< mol.GetTitle() << " in file \""
						<< e.file.string() << "\": " << e.reason
						<< '\n';
				continue;
			}
		}
		return false; //no valid molecules read
	}
		break;
	}
	return false; //shouldn't get here
}
