#include <openbabel/mol.h>
#include <set>
#include "cnn_scorer.h"
#include "molgetter.h"

struct vis_options
{
  std::string ligand_name;
  std::string receptor_name;
  std::string receptor_output;
  std::string ligand_output;

  bool frags_only;
  bool atoms_only;
  bool verbose;

  float box_size;

  //vis_options(): box_size(23.5), frags_only(false), atoms_only(false), verbose(false) {}
};

class cnn_visualization
{
    public:
    cnn_visualization(const vis_options &visopts, const cnn_options &cnnopts, FlexInfo &finfo, tee &log, const vec &center);
    void color();
    void print();

    private:
    std::string rec_string, lig_string, rec_pdb_string, lig_pdb_string; 
    OpenBabel::OBMol lig_mol, rec_mol;
    float size, original_score;
    float cenCoords [3];
    vis_options visopts;
    cnn_options cnnopts;
    FlexInfo* finfo;
    tee* log;
    const vec* center;
    bool frags_only, atoms_only,  verbose;

    void process_molecules();
    float remove_and_score(std::vector<bool> removeList, bool isRec);
    void ligCenter();
    float score(const std::string &molString, bool isRec);
    void write_scores(std::vector<float> scoreList, bool isRec);
    bool check_in_range(std::set<int> atomList);
    std::vector<float> transform(std::vector<float> inList);
    void remove_residues();
    void remove_each_atom();
};
