#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

#include <tuple>

#include "box.h"
#include "obmolopener.h"
#include "tee.h"
#include "user_opts.h"

#include <unordered_map>

#ifndef GNINA_COVINFO_H
#define GNINA_COVINFO_H

struct CovOptions {
  std::string covalent_rec_atom;
  std::string covalent_lig_atom_pattern;
  covalent_optimization covalent_optimize_lig = None;
  int bond_order = 1;
  std::string covalent_lig_atom_position;
  bool covalent_fix_lig_atom_position = false;
};

/* Parse covalent docking options and store relevant info.
 */
class CovInfo {

  tee *log = nullptr;

  char ratom_chain;
  int ratom_num = 0;
  std::string ratom_res; // optional
  std::string ratom_name;

  vec ratom_xyz = vec(0, 0, 0);

  OpenBabel::OBSmartsPattern latom_pat;

  int bond_order = 1;
  covalent_optimization optlevel = None;
  vec latom_pos = vec(0, 0, 0); // where to position latom, optional
  bool latom_pos_set = false;
  bool fix_latom_pos = false;
  

  bool initialized = false;

public:
  CovInfo() {}

  CovInfo(const CovOptions& copts, tee &l);

  bool has_content() const { return initialized; }

  bool is_rec_atom(OpenBabel::OBAtom *a) const;
  OpenBabel::OBAtom *find_rec_atom(OpenBabel::OBMol &mol) const;

  std::string rec_atom_string() const;

  bool has_user_lig_atom_pos() const { return latom_pos_set; }
  vec lig_atom_pos(OpenBabel::OBAtom *ra, OpenBabel::OBAtom *la) const;

  std::vector<std::vector<int> > get_matches(OpenBabel::OBMol &mol);

  int get_bond_order() const {
    return bond_order;
  }

  covalent_optimization get_opt_level() const {
    return optlevel;
  }
private:
  bool getXYZ(const std::string &str, vec &xyz) const;
};

#endif /* GNINA_COVINFO_H */
