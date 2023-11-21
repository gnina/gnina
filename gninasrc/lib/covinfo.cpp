#include "covinfo.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include <openbabel/elements.h>
#include <openbabel/obiter.h>
#include <openbabel/residue.h>

using namespace boost;
using namespace std;
using namespace OpenBabel;

CovInfo::CovInfo(const CovOptions& copt, tee &l) : 
    log(&l), bond_order(copt.bond_order), optlevel(copt.covalent_optimize_lig), fix_latom_pos(copt.covalent_fix_lig_atom_position) {

  if (copt.covalent_rec_atom.size() == 0)
    return; // not set

  // rec atom is either chain:resid[icode]:atomname, chain:resid:resname:atomname or
  // x,y,z
  if (!getXYZ(copt.covalent_rec_atom, ratom_xyz)) {
    vector<string> tokens;

    boost::regex expr{"([^:]+):(-?\\d+)(\\w?):([^:]+)(?::([^:]+))?"};
    boost::smatch what;

    if(boost::regex_search(copt.covalent_rec_atom, what, expr)) {
      if(what.str(1).length() > 0) ratom_chain = what.str(1)[0];
      if(what.str(1).length() > 1) {
        throw usage_error("Chain identifiers with more than one character are not supported in covalent_rec_atom: "+copt.covalent_rec_atom);
      }
      ratom_num = lexical_cast<int>(what.str(2));
      ratom_icode = what.str(3);
      if(what[5].length() > 0) {
        ratom_res = what.str(4); 
        ratom_name = what.str(5); 
      } else {
        ratom_name = what.str(4);
      }
    } else {
        throw usage_error("Could not parse covalent_rec_atom: " + copt.covalent_rec_atom);
    }
    
  }

  // lig atom pattern is smarts
  if (!latom_pat.Init(copt.covalent_lig_atom_pattern)) {
    throw usage_error("Could not parse covalent_lig_atom_pattern: " + copt.covalent_lig_atom_pattern);
  }

  // lig atom position is optional, and is x,y,z
  latom_pos_set = getXYZ(copt.covalent_lig_atom_position, latom_pos);

  if( fix_latom_pos && !latom_pos_set) {
    *log << "WARNING: covalent_fix_lig_atom_position set without specifying covalent_lig_position. Ignoring\n";
  }

  initialized = true;
}

// return true if a is the cov_rec atom
bool CovInfo::is_rec_atom(OpenBabel::OBAtom *a) const {
  if (!initialized)
    return false;
  if (ratom_name.size() > 0) {
    OBResidue *r = a->GetResidue();
    if (r->GetNum() == ratom_num && r->GetChain() == ratom_chain) {
      if (ratom_res.size() == 0 || r->GetName() == ratom_res) {
        if (trim_copy(r->GetAtomID(a)) == ratom_name) {
          if (ratom_icode.length() == 0 || r->GetInsertionCode() == ratom_icode[0]) {
            return true;
          }
        }
      }
    }
  } else { // use cartesian coordinates
    vec c(a->GetX(), a->GetY(), a->GetZ());
    c -= ratom_xyz;
    float distsq = c.norm_sqr();
    if (distsq < 0.05) {
      return true;
    }
  }
  return false;
}

OpenBabel::OBAtom *CovInfo::find_rec_atom(OpenBabel::OBMol &mol) const {
  OBAtom *covrec = nullptr;
  FOR_ATOMS_OF_MOL(a, mol) {
    if (is_rec_atom(&*a)) {
      covrec = &*a;
      break; //first in file
    }
  }
  return covrec;
}

// return string representing receptor atom
std::string CovInfo::rec_atom_string() const {
  if (!initialized)
    return "";
  stringstream ss;
  if (ratom_name.size()) {
    ss << ratom_chain << ":" << ratom_num << ":";
    if (ratom_res.size())
      ss << ratom_res << ":";
    ss << ratom_name;
  } else {
    ratom_xyz.print(ss);
  }
  return ss.str();
}

// return true if string can be parsed as Cartesian coords
bool CovInfo::getXYZ(const string &str, vec &xyz) const {
  vector<string> tokens;
  boost::split(tokens, str, boost::is_any_of(","));
  if (tokens.size() == 3) {
    try {
      for (unsigned i = 0; i < 3; i++) {
        xyz[i] = lexical_cast<float>(tokens[i]);
      }
      return true;
    } catch (const bad_lexical_cast &e) {
      return false;
    }
  }
  return false;
}

//get covalent radius adjusted for hybridization
static double get_cov_rad(OBAtom *a) {
  if(a->GetAtomicNum() == 6) {
    //https://en.wikipedia.org/wiki/Covalent_radius
    int hyb = a->GetHyb();
    if(hyb == 1) return 0.69;
    if(hyb == 2) return 0.73;
  }
  return OBElements::GetCovalentRad(a->GetAtomicNum());
}

// return position for ligand atom la covalently bonding to receptor atom ra
vec CovInfo::lig_atom_pos(OBAtom *ra, OBAtom *la) const {
  vec recpos(ra->GetX(), ra->GetY(), ra->GetZ());
  float cdist = get_cov_rad(ra) + get_cov_rad(la);

  if (latom_pos_set) {
    // override from user
    float dist = (latom_pos - recpos).norm();
    if (dist > 1.5 * cdist) {
      *log << "WARNING: Large covalent bond distance using specified "
              "covalent_lig_atom_position: "
           << dist << "\n";
    }
    return latom_pos;

  } else {
    vector3 ret;
    ra->GetNewBondVector(ret, cdist);
    return vec(ret.x(), ret.y(), ret.z());
  }
}

// return matchs for smarts on ligand
std::vector<std::vector<int>> CovInfo::get_matches(OpenBabel::OBMol &mol) {
  vector<std::vector<int>> ret;

  if (!latom_pat.Match(mol)) {
    return ret;
  }

  return latom_pat.GetUMapList();
}