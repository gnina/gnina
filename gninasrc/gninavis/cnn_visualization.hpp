#include <openbabel/mol.h>
#include <set>
#include <unordered_set>
#include "cnn_scorer.h"
#include "molgetter.h"
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>

struct vis_options {
    std::string ligand_name;
    std::string receptor_name;
    bool skip_receptor_output;
    bool skip_ligand_output;
    std::string additivity;
    std::string layer_to_ignore;
    std::string target;

    bool frags_only;
    bool atoms_only;
    bool verbose;
    bool output_files;
    bool skip_bound_check;
    bool zero_values;
    int gpu;

    bool outputdx;
    float box_size;
    double score_scale;

    vis_options()
        : frags_only(false), atoms_only(false), verbose(false),
            output_files(false), skip_bound_check(false), outputdx(false),
            gpu(0), box_size(23.5), score_scale(10) {
    }
};

class cnn_visualization {
  public:
    cnn_visualization(const vis_options &visopts, const cnn_options &cnnopts,
        const vec &center);
    void setup();
    void lrp();
    void gradient_vis();
    void masking();
    void print();

  private:
    std::string rec_string, lig_string, rec_pdb_string, lig_pdb_string;
    OpenBabel::OBMol lig_mol, rec_mol;
    float size, original_score;
    float cenCoords[3];
    vis_options visopts;
    cnn_options cnnopts;
    vec center;
    model unmodified_receptor;
    model unmodified_ligand;
    bool frags_only, atoms_only, verbose;
    double score_scale;

    std::unordered_map<std::string, int> rec_map;
    std::unordered_map<std::string, int> lig_map;

    std::string original_rec_string;
    std::string original_lig_string;

    void populate_coordinate_map(const std::string& molstring,
        std::unordered_map<std::string, int>& map);
    void process_molecules();
    std::string modify_pdbqt(
        const std::unordered_set<std::string> &atoms_to_remove, bool isRec);
    float score_modified_receptor(const std::string &modified_rec_string);
    float score_modified_ligand(const std::string &modified_lig_string);

    float score(const std::string &molString, bool isRec);
    void write_scores(const std::unordered_map<std::string, float> scores,
        bool isRec, std::string method);
    bool check_in_range(const std::unordered_set<std::string> &atom_xyzs);
    float transform_score_diff(float diff_val);
    void remove_residues();
    std::unordered_map<std::string, float> remove_each_atom();
    void write_additivity(
        const std::unordered_map<std::string, float> &single_score_diffs,
        const std::unordered_map<std::string, float> &frag_score_diffs);
    std::unordered_map<std::string, float> remove_fragments(int size);
    void remove_ligand_atoms();
    void add_adjacent_hydrogens(
        std::unordered_set<std::string> &atoms_to_remove, bool isRec);
    int get_openbabel_index(const std::string &xyz, bool rec);
    std::string get_xyz_from_index(int index, bool rec);
    std::string get_xyz(const std::string &line);
    void print_vector(const std::vector<int> &atoms_to_remove);

};
