#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include <tuple>

#include "obmolopener.h"
#include "box.h"
#include "tee.h"

#include <unordered_map>

#ifndef SMINA_FLEXINFO_H
#define SMINA_FLEXINFO_H

/* Store information for identifying flexible residues in receptor and
 * provide routines for extracting these residues as needed.
 */
class FlexInfo {
    double flex_dist;
    boost::unordered_set<std::tuple<char, int, char> > residues; //chain, resid, insertion code
    OpenBabel::OBMol distligand;
    tee& log;

    int nflex = 0;;
    bool nflex_hard_limit = false;
    bool full_flex_output = false;

    char defaultch = ' '; // Default chain

    static std::map<std::string, int> num_heavy_atoms_per_residue;

  public:
    FlexInfo(tee& l)
        : flex_dist(0), log(l) {
    }
    FlexInfo(const std::string& flexres, double flexdist,
        const std::string& ligand, int nflex, bool nflex_hard_limit, 
        bool full_flex_output, tee& l);
    bool has_content() const {
      return residues.size() > 0 || flex_dist > 0;
    }

    void extract_residue(OpenBabel::OBMol& rigid,  OpenBabel::OBResidue *r, OpenBabel::OBMol& flexres, bool fullres=false);
    void extractFlex(OpenBabel::OBMol& receptor, OpenBabel::OBMol& rigid,
        std::string& flexpdbqt);

    bool omit_residue(std::tuple<char, int, char> r);
    std::string residue_to_pdbqt(OpenBabel::OBMol& flex);

    void print_flex() const;

    bool full_output() const { return full_flex_output;}

  private:
    void sanitizeResidues(OpenBabel::OBMol& receptor); //remove inflexible residues from residues set
    void keepNearestResidues(
      const OpenBabel::OBMol& 
      rigid, std::vector<std::size_t>& residues_idxs);
    bool isSideChain(std::string aid);

    void checkResidue(OpenBabel::OBResidue* residue); // Check number of atoms per residue is reasonable
};

#endif /* SMINA_FLEXINFO_H */
