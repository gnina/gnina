#pragma once
#include "common.h"
#include "random.h"
#include <string>

//enum options and their parsers
enum ApproxType
{
  LinearApprox, SplineApprox, Exact, GPU
};

std::istream& operator>>(std::istream& in, ApproxType& type);

enum pose_sort_order {
  CNNscore,
  CNNaffinity,
  Energy
};

//for reading in as a commandline option
std::istream& operator>>(std::istream &in, pose_sort_order &sort_order);

enum cnn_scoring_level {
  CNNnone, //don't use CNN
  CNNrescore, // use CNN only for final scoring and ranking
  CNNrefinement, // use CNN only for minimization
  CNNall // use CNN everywhere
};

std::istream& operator>>(std::istream &in, cnn_scoring_level &cnn_level);


struct cnn_options {
    //stores options associated with cnn scoring
    std::vector<std::string> cnn_models; //path(s) to model file
    std::vector<std::string> cnn_weights; //weights for model
    std::string cnn_recmap; //optional file specifying receptor atom typing to channel map
    std::string cnn_ligmap; //optional file specifying ligand atom typing to channel map
    std::vector<std::string> cnn_model_names; // name of builtin model
    vec cnn_center;
    fl resolution; //this isn't specified in model file, so be careful about straying from default
    unsigned cnn_rotations; //do we want to score multiple orientations?
    cnn_scoring_level cnn_scoring;
    double subgrid_dim;
    bool outputdx;
    bool outputxyz;
    bool gradient_check;
    bool move_minimize_frame;  //recenter with every scoring evaluation
    bool fix_receptor;
    bool verbose;
    std::string xyzprefix;
    unsigned seed; //random seed

    cnn_options()
        : cnn_center(NAN, NAN, NAN),
           resolution(0.5), cnn_rotations(0), cnn_scoring(CNNrescore),
            subgrid_dim(0.0), outputdx(false),
            outputxyz(false), gradient_check(false), move_minimize_frame(false),
            fix_receptor(false), verbose(false), seed(0) {
    }

    bool moving_receptor() const {
      //doesn't make sense to accumulate transformation gradient with moving center
      if(move_minimize_frame) return false;
      if(fix_receptor) return false; //just aren't doing it
      return true;
    }
};

//just a collection of user-specified configurations
struct user_settings {
    sz num_modes;
    fl out_min_rmsd;
    fl forcecap;
    int seed;
    int verbosity;
    int cpu;
    int device; //gpu number

    int exhaustiveness;
    int num_mc_steps;
    pose_sort_order sort_order;

    bool score_only;
    bool randomize_only;
    bool local_only;
    bool dominimize;
    bool include_atom_info;
    bool gpu_on;

    cnn_options cnnopts;

    //reasonable defaults
    user_settings()
        :  num_modes(9), out_min_rmsd(1), forcecap(1000),
            seed(auto_seed()), verbosity(1), cpu(1), device(0),
            exhaustiveness(10), num_mc_steps(0), sort_order(CNNscore), score_only(false),
            randomize_only(false), local_only(false), dominimize(false),
            include_atom_info(false), gpu_on(false) {

    }
};

