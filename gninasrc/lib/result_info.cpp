/*
 * result_info.cpp
 *
 *  Created on: Jun 10, 2014
 *      Author: dkoes
 */

#include "result_info.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/generic.h>
#include <cstring>
#include <boost/lexical_cast.hpp>

void result_info::setMolecule(const model& m) {
  std::stringstream str;
  if (m.write_sdf(str)) {
    sdfvalid = true; //can do native sdf output, //TODO - fix flex residue output - todone? seems like it works
    molstr = str.str();
  } else {
    m.write_ligand(str);
    molstr = str.str();
  }

  if (m.num_flex() > 0) { //save flex residue info
    std::stringstream fstr;
    m.write_flex(fstr);
    flexstr = fstr.str();
  }
  name = m.get_name();
}

//write a table (w/header) of per atom values to out
void result_info::writeAtomValues(std::ostream& out,
    const weighted_terms *wt) const {
  const terms *t = wt->unweighted_terms();
  out << "atomid el pos";
  std::vector<std::string> names = t->get_names(true);
  for (unsigned j = 0, m = names.size(); j < m; j++) {
    out << " " << names[j];
  }
  out << "\n";
  out << atominfo;
}

//computes per-atom term values and formats them into the atominfo string
void result_info::setAtomValues(const model& m, const weighted_terms *wt) {
  std::vector<flv> values;
  const terms *t = wt->unweighted_terms();
  t->evale_robust(m, values);
  std::stringstream str;
  vecv coords = m.get_ligand_coords();
  assert(values.size() == coords.size());
  for (unsigned i = 0, n = values.size(); i < n; i++) {
    assert(values[i].size() == t->size());
    str << m.ligand_atom_str(i) << " ";
    coords[i].print(str);
    for (unsigned j = 0, m = values[i].size(); j < m; j++) {
      str << " " << values[i][j] * wt->weight(j);
    }
    str << "\n";
  }
  str << "END\n";
  atominfo = str.str();
}

//helper function for setting molecular data that will wrap data in remark
//if output format is pdb; copy by value because we might change things
static void setMolData(OpenBabel::OBFormat *format, OpenBabel::OBMol& mol,
    std::string attr, std::string value) {
  using namespace OpenBabel;
  OBPairData* sddata = new OBPairData();
  if (strcmp(format->GetID(), "pdb") == 0 || strcmp(format->GetID(), "ent") == 0
      || strcmp(format->GetID(), "pdbqt") == 0) {
    //note that pdbqt currently ignores these.. hopefully this will be fixed
    //put attribute name into value so attr can be REMARK
    value = " " + attr + " " + value;
    attr = "REMARK";
  }

  sddata->SetAttribute(attr);
  sddata->SetValue(value);
  mol.SetData(sddata);
}

//output flexible residue conformers
void result_info::writeFlex(std::ostream& out, std::string& ext, int modelnum) {
  using namespace OpenBabel;
  OBMol mol;
  OBConversion outconv;
  OBFormat *format = outconv.FormatFromExt(ext);
  if(!format) {
    throw usage_error("Invalid formt "+ext);
  }
  //residue are pdbqt only currently
  outconv.SetInFormat("PDBQT");
  outconv.SetOutFormat(format);
  outconv.ReadString(&mol, flexstr); //otherwise keep orig mol
  mol.SetTitle(name); //same name as cmpd
  outconv.SetLast(false);
  outconv.SetOutputIndex(modelnum + 1); //for pdb multi model output, workaround OB bug with ignoring model 1

  outconv.Write(&mol, &out);
}

//output molecular data
//ideally, we will deal natively in sdf and only use openbabel to convert for alternative formats
void result_info::write(std::ostream& out, std::string& ext,
    bool include_atom_terms, const weighted_terms *wt, int modelnum) {
  using namespace OpenBabel;
  OBMol mol;
  OBConversion outconv;

  OBFormat *format = outconv.FormatFromExt(ext);

  if(!format) {
    throw usage_error("Invalid format: "+ext);
  }
  if (sdfvalid && strcmp(format->GetID(), "sdf") == 0) { //use native sdf
    out << molstr;
    //now sd data
    out << "> <minimizedAffinity>\n";
    out << std::fixed << std::setprecision(5) << energy << "\n\n";

    if (rmsd >= 0) {
      out << "> <minimizedRMSD>\n";
      out << std::fixed << std::setprecision(5) << rmsd << "\n\n";
    }

    if (cnnscore >= 0) {
      out << "> <CNNscore>\n";
      out << std::fixed << std::setprecision(10) << cnnscore << "\n\n";
    }

    if (cnnaffinity != 0) {
      out << "> <CNNaffinity>\n";
      out << std::fixed << std::setprecision(10) << cnnaffinity << "\n\n";
    }

    if (include_atom_terms) {
      std::stringstream astr;
      writeAtomValues(astr, wt);
      out << "> <atomic_interaction_terms>\n";
      out << astr.str() << "\n\n";
    }

    out << "$$$$\n";
  } else
    if (!sdfvalid && ext == ".pdbqt") {
      out << "MODEL " << boost::lexical_cast<std::string>(modelnum) << "\n";
      out << "REMARK minimizedAffinity "
          << boost::lexical_cast<std::string>((float) energy);
      if (rmsd >= 0)
        out << "REMARK minimizedRMSD "
            << boost::lexical_cast<std::string>((float) rmsd);
      if (cnnscore >= 0)
        out << "REMARK CNNscore "
            << boost::lexical_cast<std::string>((float) cnnscore);
      if (cnnaffinity != 0)
        out << "REMARK CNNaffinity "
            << boost::lexical_cast<std::string>((float) cnnaffinity);
      out << molstr;
      out << "ENDMDL\n";
    } else //convert with openbabel
    {
      if (sdfvalid)
        outconv.SetInFormat("SDF");
      else
        outconv.SetInFormat("PDBQT");
      outconv.SetOutFormat(format);

      outconv.ReadString(&mol, molstr); //otherwise keep orig mol
      mol.DeleteData(OBGenericDataType::PairData); //remove remarks

      setMolData(format, mol, "minimizedAffinity",
          boost::lexical_cast<std::string>((float) energy));

      if (rmsd >= 0) {
        setMolData(format, mol, "minimizedRMSD",
            boost::lexical_cast<std::string>((float) rmsd));
      }
      if (cnnscore >= 0) {
        setMolData(format, mol, "CNNscore",
            boost::lexical_cast<std::string>((float) cnnscore));
      }
      if (cnnaffinity != 0) {
        setMolData(format, mol, "CNNaffinity",
            boost::lexical_cast<std::string>((float) cnnaffinity));
      }

      if (include_atom_terms) {
        std::stringstream astr;
        writeAtomValues(astr, wt);
        setMolData(format, mol, "atomic_interaction_terms", astr.str());
      }
      mol.SetTitle(name); //otherwise lose space separated names
      outconv.SetLast(false);
      outconv.SetOutputIndex(modelnum + 1); //for pdb multi model output, workaround OB bug with ignoring model 1
      outconv.SetOutStream(&out, false);
      format->WriteMolecule(&mol, &outconv); //conv ignores isLast
    }
}
