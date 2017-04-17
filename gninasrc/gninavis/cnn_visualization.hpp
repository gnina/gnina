#include <openbabel/mol.h>
#include <set>
#include <unordered_set>
#include "cnn_scorer.h"
#include "molgetter.h"
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

struct vis_options
{
  std::string ligand_name;
  std::string receptor_name;
  std::string receptor_output;
  std::string ligand_output;
  std::string additivity;

  bool frags_only;
  bool atoms_only;
  bool verbose;
  bool output_files;
  bool skip_bound_check;
  int gpu;

  bool outputdx;
  float eps;
  float box_size;

  vis_options(): frags_only(false), atoms_only(false), verbose(false),
		  output_files(false), skip_bound_check(false), outputdx(false), gpu(0), box_size(23.5) {}
};

class cnn_visualization
{
    public:
    cnn_visualization(const vis_options &visopts, const cnn_options &cnnopts, const vec &center);
    void lrp();
    void removal();
    void print();

    private:
    std::string rec_string, lig_string, rec_pdb_string, lig_pdb_string; 
    OpenBabel::OBMol lig_mol, rec_mol;
    float size, original_score;
    float cenCoords [3];
    vis_options visopts;
    cnn_options cnnopts;
    const vec* center;
    model unmodified_receptor;
    model unmodified_ligand;
    bool frags_only, atoms_only,  verbose;

    void process_molecules();
    std::string modify_pdbqt(std::vector<int> atoms_to_remove, bool isRec);
    float score_modified_receptor(const std::string &modified_rec_string);
    float score_modified_ligand(const std::string &modified_lig_string);

    std::vector<std::string> rec_map;
    std::vector<std::string> lig_map;
    float score(const std::string &molString, bool isRec);
    void write_scores(std::vector<float> scoreList, bool isRec, bool removal);
    bool check_in_range(std::unordered_set<int> atomList);
    float transform_score_diff(float diff_val);
    void remove_residues();
    std::vector<float> remove_each_atom();
    void output_modified_string(const std::string &modified_string, const std::vector<int> &atoms_removed,
                                  bool receptor);
    void write_additivity(std::vector<float> single_score_diffs, std::vector<float> frag_score_diffs);
    std::vector<float> remove_fragments(int size);
    void remove_ligand_atoms();
    void add_adjacent_hydrogens(std::vector<int> &atoms_to_remove, bool isRec);
    int get_openbabel_index(double x_coordinate);
    void print_vector(const std::vector<int> &atoms_to_remove);

};
