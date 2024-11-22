#include "flexinfo.h"
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/generic.h>

#include <exception>
#include <limits>

using namespace std;
using namespace OpenBabel;

std::map<std::string, int> FlexInfo::num_heavy_atoms_per_residue = {
    {"ARG", 12}, {"HIS", 11}, {"LYS", 10}, {"ASP", 9}, {"GLU", 10}, {"SER", 7},
    {"THR", 8},  {"ASN", 9},  {"GLN", 9},  {"CYS", 7}, {"SEC", 7},  {"GLY", 5},
    {"PRO", 8},  {"ALA", 6},  {"VAL", 8},  {"ILE", 9}, {"LEU", 9},  {"MET", 9},
    {"PHE", 12}, {"TYR", 13}, {"TRP", 15},
};

FlexInfo::FlexInfo(const std::string &flexres, double flexdist,
                   const std::string &ligand, int nflex_,
                   bool nflex_hard_limit_, bool full_flex_output_, tee &l)
    : flex_dist(flexdist), nflex(nflex_), nflex_hard_limit(nflex_hard_limit_),
      full_flex_output(full_flex_output_), log(l) {
  using namespace OpenBabel;
  using namespace std;
  // first extract comma separated list
  if (flexres.size() > 0) {
    vector<string> tokens;
    boost::split(tokens, flexres, boost::is_any_of(","));

    vector<string> chres;
    for (unsigned i = 0, n = tokens.size(); i < n; i++) {
      // each token may be either chain:resid or just resid
      string tok = tokens[i];
      boost::split(chres, tok, boost::is_any_of(":"));
      char chain = 0;
      int resid = 0;
      char icode = 0;
      if (chres.size() >= 2) {
        if (chres[0].size() != 1)
          log << "WARNING: chain specification not single character "
              << chres[0] << "\n";
        chain = chres[0][0]; // if empty will be null which is what is desired
        resid = boost::lexical_cast<int>(chres[1]);
        if (chres.size() == 3) {      // Insertion code is present
          if (chres[2].size() == 1) { // Check that icode is single char
            icode = chres[2][0];
          } else { // Invalid icode
            log << "WARNING: ignoring invalid chain:resid:icode specifier "
                << tok << "\n";
            continue;
          }
        }
      } else if (chres.size() == 1) {
        resid = boost::lexical_cast<int>(chres[0]);
      } else {
        log << "WARNING: ignoring invalid chain:resid specifier " << tok
            << "\n";
        continue;
      }

      residues.insert(tuple<char, int, char>(chain, resid, icode));
    }
  }

  if (ligand.size() > 0 && flex_dist > 0) {
    // next read ligand for distance checking
    obmol_opener opener;
    OBConversion conv;
    opener.openForInput(conv, ligand);
    conv.Read(&distligand); // first ligand only
  }
}

// return true if residue isn't flexible
static bool isInflexible(const string &resname) {
  return resname == "ALA" || resname == "GLY" || resname == "PRO";
}

// remove inflexible residues from residues set
// the receptor is needed because we don't store the residue name
void FlexInfo::sanitizeResidues(OpenBabel::OBMol &receptor) {
  using namespace OpenBabel;
  if (!has_content())
    return;

  OBResidueIterator ritr, rend;
  OBResidue *firstres = receptor.BeginResidue(ritr);
  char defaultch = ' ';
  if (firstres)
    defaultch = firstres->GetChain();

  // Iterate over all receptor residues
  for (ritr = receptor.BeginResidues(), rend = receptor.EndResidues();
       ritr != rend; ++ritr) {
    OBResidue *r = *ritr;

    char ch = r->GetChain();
    int resid = r->GetNum();
    char icode = r->GetInsertionCode();
    std::string resname = r->GetName();
    if (ch == 0)
      ch = defaultch; // substitute default chain for unspecified chain
    tuple<char, int, char> res(ch, resid, icode);
    if (residues.count(res) >
        0) { // Residue in user-specified flexible residues
      if (isInflexible(resname)) { // Residue can't be flexible
        residues.erase(res); // Remove residue from list of flexible residues
        log << "WARNING: Removing non-flexible residue " << resname;
        log << " " << ch << ":" << resid << ":" << icode << "\n";
      }
    }
  }
}

void FlexInfo::checkResidue(OpenBabel::OBResidue *residue) {
  using namespace OpenBabel;

  if (residue == nullptr) {
    return;
  }

  double n_atoms_with_buffer = 0.0;

  std::string resname = residue->GetName();

  if (num_heavy_atoms_per_residue.count(resname)) { // Check key exists
    n_atoms_with_buffer = 1.5 * num_heavy_atoms_per_residue[resname];
  } else { // Non standard residue name
    n_atoms_with_buffer = std::numeric_limits<double>::max();
  }

  // TODO: Update with OBResidue::GetNumHvyAtoms() [openbabel#2299]
  unsigned int num_hvy_atoms = 0;
  for (auto atom = residue->BeginAtoms(); atom != residue->EndAtoms(); atom++) {
    if ((*atom)->GetAtomicNum() != OBElements::Hydrogen) {
      num_hvy_atoms++;
    }
  }

  if (num_hvy_atoms > n_atoms_with_buffer) {
    char ic = residue->GetInsertionCode();

    log << "WARNING: Residue " << residue->GetChain() << ":"
        << residue->GetNum();
    if (ic != '\0') {
      log << ":" << ic;
    }
    log << " appears to have too many atoms.\n";
  }
}

void FlexInfo::keepNearestResidues(const OpenBabel::OBMol &receptor,
                                   std::vector<std::size_t> &residues_idxs) {

  using namespace OpenBabel;

  // Loop over residue list and compute minimum distances
  std::vector<std::pair<std::size_t, double>> min_distance(
      residues_idxs.size());
  for (std::size_t i{0}; i < residues_idxs.size(); i++) {

    min_distance[i].first = residues_idxs[i];
    min_distance[i].second = std::numeric_limits<double>::max();

    std::size_t residx = min_distance[i].first;
    OBResidue *res = receptor.GetResidue(residx);

    // Minimum distance between ligand and current residue atoms
    double d2;

    // Loop over ligand atoms
    FOR_ATOMS_OF_MOL(alig, distligand) {
      vector3 al = alig->GetVector();

      // Loop over residue atoms
      for (const auto ares : res->GetAtoms()) {

        if (!isSideChain(res->GetAtomID(ares)))
          continue; // skip rigid backbone

        vector3 ar = ares->GetVector();

        d2 = al.distSq(ar);

        if (d2 < min_distance[i].second) {
          min_distance[i].second = d2;
        }
      }
    }
  }

  // Sort by distance
  std::sort(min_distance.begin(), min_distance.end(),
            [](const std::pair<std::size_t, double> &left,
               std::pair<std::size_t, double> &right) {
              return left.second < right.second;
            });

  residues.clear();

  // Kep only nearest nflex residues
  std::size_t res_added = 0;
  for (const auto &resdist : min_distance) {
    if (res_added == nflex) {
      break;
    }

    OBResidue *residue = receptor.GetResidue(resdist.first);

    char ch = residue->GetChain();
    int resid = residue->GetNum();
    char icode = residue->GetInsertionCode();

    residues.insert(std::tuple<char, int, char>(ch, resid, icode));

    res_added++;
  }
}

bool FlexInfo::isSideChain(std::string aid) {
  boost::trim(aid);

  return aid != "CA" && aid != "N" && aid != "C" && aid != "O" && aid != "H" &&
         aid != "HN";
}

// Remove the specified residue r from rigid and put it in flexres
void FlexInfo::extract_residue(OpenBabel::OBMol &rigid, OpenBabel::OBResidue *r,
                               OpenBabel::OBMol &flex, bool fullres) {

  std::vector<OBAtom *> flexatoms; // rigid atom ptrs that should be flexible
  boost::unordered_map<OBAtom *, int>
      flexmap; // map rigid atom ptrs to atom indices in flex

  // make absolutely sure that CA is the first atom
  // first bond is rigid, so take both CA and C
  bool CAseen = false;
  for (OBAtomIterator aitr = r->BeginAtoms(), aend = r->EndAtoms();
       aitr != aend; ++aitr) {
    OBAtom *a = *aitr;
    std::string aid = r->GetAtomID(a);
    boost::trim(aid);
    if (aid == "CA" || aid == "C") {
      flexatoms.push_back(a);
      flex.AddAtom(*a);
      flexmap[a] = flex.NumAtoms(); // after addatom since indexed by
    }
    if (aid == "CA") {
      if(CAseen) {
        stringstream msg;
        msg << "Multiple copies of residue " << r->GetChain() << ":" << r->GetNum() << ".  I can't handle this situation.";
        throw runtime_error(msg.str());
      }
      CAseen = true;
    }
  }

  for (OBAtomIterator aitr = r->BeginAtoms(), aend = r->EndAtoms();
       aitr != aend; ++aitr) {
    OBAtom *a = *aitr;
    std::string aid = r->GetAtomID(a);
    boost::trim(aid);
    bool isflex = true;
    if(aid == "CA" || aid == "C") {
      isflex = false;
    } 
    else if(aid == "N" || aid == "O" || aid == "HN" || a->IsNonPolarHydrogen()) {
      if(!fullres) {
        isflex = false;
      } else {
        if(aid != "N" || a->GetExplicitValence() == 1) { // Nterminal special case of flexible backbone
          OBPairData *dp = new OBPairData;
          dp->SetAttribute("Fixed");
          a->SetData(dp);
        }
      }
    }    

    if(aid == "H") {
      //added hydrogens will all be H, only want to ignore backbone as polar hydrogens are needd for typing
      isflex = true;
      FOR_NBORS_OF_ATOM(neigh, *a){
        std::string nid = r->GetAtomID(&*neigh);
        boost::trim(nid);
        if(nid == "N") {
          isflex = false;
          break;
        }
      }
    }

    if (isflex) {
      flexatoms.push_back(a);
      flex.AddAtom(*a);
      flexmap[a] = flex.NumAtoms(); // after addatom since indexed by
      if (boost::starts_with(aid, "S") &&
          boost::starts_with(r->GetName(), "CY") && !fullres) { //no warning for covalent
        if (a->GetHvyDegree() > 1) {
          log << "WARNING: Disulfide bonds are ignored in flexible cysteine "
                 "residues.\n";
        }
      }
    }
  }

  // now add bonds - at some point I should add functions to openbabel to do
  // this..
  for (unsigned i = 0, n = flexatoms.size(); i < n; i++) {
    OBAtom *a = flexatoms[i];
    FOR_BONDS_OF_ATOM(b, a) {
      OBBond &bond = *b;
      // just do one direction, if first atom is a
      // and second atom is a flexatom need to add bond
      if (a == bond.GetBeginAtom() && flexmap.count(bond.GetEndAtom())) {
        flex.AddBond(flexmap[a], flexmap[bond.GetEndAtom()],
                     bond.GetBondOrder(), bond.GetFlags());
      }
    }
  }

  flex.AddResidue(*r);
  OBResidue *newres = flex.GetResidue(0);
  if (newres) {
    // add all atoms with proper atom ids
    for (unsigned i = 0, n = flexatoms.size(); i < n; i++) {
      OBAtom *origa = flexatoms[i];
      OBAtom *newa = flex.GetAtom(flexmap[origa]);
      newres->RemoveAtom(origa);
      newres->AddAtom(newa);
      newa->SetResidue(newres);
      std::string aid = r->GetAtomID(origa);
      newres->SetAtomID(newa, aid);
    }
  }

  flex.SetChainsPerceived(true);
  rigid.SetChainsPerceived(true);

  // remove flexatoms from rigid
  for (unsigned i = 0, n = flexatoms.size(); i < n; i++) {
    OBAtom *a = flexatoms[i];
    rigid.DeleteAtom(a);
  }
  rigid.SetChainsPerceived(true);

}

void FlexInfo::extractFlex(OpenBabel::OBMol &receptor, OpenBabel::OBMol &rigid,
                           std::string &flexpdbqt) {
  rigid = receptor;
  rigid.SetChainsPerceived(true); // OB bug workaround
  flexpdbqt.clear();

  // identify residues close to distligand here
  Box b;
  b.add_ligand_box(distligand);
  b.expand(flex_dist);
  double flsq = flex_dist * flex_dist;

  std::vector<std::size_t> residues_idxs;
  FOR_ATOMS_OF_MOL(a, rigid) {
    if (a->GetAtomicNum() == 1)
      continue;                    // heavy atoms only
    if (a->GetResidue() != NULL) { // Check if residue exists
      if (!isSideChain(a->GetResidue()->GetAtomID(&(*a))))
        continue; // skip backbone
    } else {
      throw std::runtime_error(
          "Missing residue information needed by FlexInfo::extractFlex.");
    }

    vector3 v = a->GetVector();
    if (b.ptIn(v.x(), v.y(), v.z())) {
      // in box, see if any atoms are close enough
      FOR_ATOMS_OF_MOL(b, distligand) {
        vector3 bv = b->GetVector();
        if (v.distSq(bv) < flsq) {
          // process residue
          OBResidue *residue = a->GetResidue();
          if (residue) {
            checkResidue(residue);

            char ch = residue->GetChain();
            int resid = residue->GetNum();
            char icode = residue->GetInsertionCode();
            if (!isInflexible(residue->GetName())) {
              // Store index of flexible residue for retrival
              std::tuple<char, int, char> key(ch, resid, icode);
              if( residues.count(key) == 0) {
                residues_idxs.push_back(residue->GetIdx());
                residues.insert(std::tuple<char, int, char>(ch, resid, icode));
              }
            }
          }
          break;
        }
      }
    }
  }

  // replace any empty chains with first chain in mol
  OBResidueIterator ritr;
  OBResidue *firstres = rigid.BeginResidue(ritr);
  if (firstres)
    defaultch = firstres->GetChain();

  sanitizeResidues(receptor);

  if (nflex > -1 && residues.size() > nflex && nflex_hard_limit) {
    throw std::runtime_error(
        "Number of flexible residues found is higher than --flex_limit.");
  } else if (nflex > -1 && residues.size() > nflex) {
    log << "WARNING: Only the flex_max residues closer to the ligand are "
           "considered as flexible.\n";

    keepNearestResidues(receptor, residues_idxs);
  }

  std::vector<std::tuple<char, int, char>> sortedres(residues.begin(),
                                                     residues.end());
  for (unsigned i = 0, n = sortedres.size(); i < n; i++) {
    if (get<0>(sortedres[i]) == 0)
      get<0>(sortedres[i]) = defaultch;
  }

  sort(sortedres.begin(), sortedres.end());

  // reinsert residues now with default chain
  residues.clear();
  residues.insert(sortedres.begin(), sortedres.end());

  rigid.BeginModify();
  int flexcnt = 0;
  boost::unordered_set<std::tuple<char, int, char> > seen;

  // identify atoms that have to be in flexible component
  // this is the side chain and CA, but _not_ the C and N
  for (OBResidueIterator ritr = rigid.BeginResidues(),
                         rend = rigid.EndResidues();
       ritr != rend; ++ritr) {
    OBResidue *r = *ritr;
    char ch = r->GetChain();
    int resid = r->GetNum();
    char icode = r->GetInsertionCode();
    std::tuple<char, int, char> chres(ch, resid, icode);

    if (residues.count(chres)) {
      if(seen.count(chres)) {
        stringstream msg;
        msg << "Multiple copies of residue " << r->GetChain() << ":" << r->GetNum() << ".  I can't handle this situation.";
        throw runtime_error(msg.str());
      }
      flexcnt++;
      seen.insert(chres);
      // create a separate molecule for each flexible residue
      OBMol flex;
      extract_residue(rigid, r, flex);
      flexpdbqt += residue_to_pdbqt(flex);
    } // end if residue
  }

  if (flexcnt != residues.size()) {
    log << "WARNING: Only identified " << flexcnt << " of " << residues.size()
        << " requested flexible residues.\n";
  }
  rigid.EndModify(false); //don't nuke chain perception

}

void FlexInfo::print_flex() const {

  // Residues are stored as unordered_set
  // Sort before printing
  std::vector<std::tuple<char, int, char>> sortedres(residues.begin(),
                                                     residues.end());
  sort(sortedres.begin(), sortedres.end());

  if (sortedres.size() > 0) {
    log << "Flexible residues:";
    for (unsigned i = 0, n = sortedres.size(); i < n; i++) {
      log << " " << get<0>(sortedres[i]) << ":" << get<1>(sortedres[i]);
      if (get<2>(sortedres[i]) > 0)
        log << ":" << get<2>(sortedres[i]);
    }
    log << "\n";
  }
}

std::string FlexInfo::residue_to_pdbqt(OpenBabel::OBMol &flex) {
  OBConversion conv;
  conv.SetOutFormat("PDBQT");
  conv.AddOption("s", OBConversion::OUTOPTIONS); // flexible residue

  OBResidue *newres = flex.GetResidue(0);
  std::string flexres = conv.WriteString(&flex);
  std::string flexpdbqt;

  if (newres) {
    // the pdbqt writing code breaks flex into fragments, in the process it
    // loses all residue information so we rewrite the strings..
    std::vector<std::string> tokens;
    std::string resn = newres->GetName();
    while (resn.size() < 3)
      resn += " ";

    std::string resnum = boost::lexical_cast<std::string>(newres->GetNum());
    while (resnum.size() < 4)
      resnum = " " + resnum;

    char ch = newres->GetChain();
    char icode = newres->GetInsertionCode();
    boost::split(tokens, flexres, boost::is_any_of("\n"));
    for (unsigned i = 0, n = tokens.size(); i < n; i++) {
      std::string line = tokens[i];
      if (line.size() > 25) {
        // replace UNL with resn
        for (unsigned p = 0; p < 3; p++) {
          line[17 + p] = resn[p];
        }
        // resid
        for (unsigned p = 0; p < 4; p++) {
          line[22 + p] = resnum[p];
        }
        if (icode > 0)
          line[26] = icode;
        line[21] = ch;
      }

      if (line.size() > 0) {
        flexpdbqt += line;
        flexpdbqt += "\n";
      }
    }
  } else {
    flexpdbqt += flexres;
  }
  return flexpdbqt;
}

// remove specified residue (presumably covalent residue) from requested residues
// return true if removed
bool FlexInfo::omit_residue(std::tuple<char, int, char> r) {
  if(residues.count(r)) {
    residues.erase(r);
    return true;
  }
  return false;
}
