#pragma once
#include "common.h"
#include <string>

struct cnn_options {
    //stores options associated with cnn scoring
    std::string cnn_model; //path to model file
    std::string cnn_weights; //weights for model
    std::string cnn_recmap; //optional file specifying receptor atom typing to channel map
    std::string cnn_ligmap; //optional file specifying ligand atom typing to channel map
    std::string cnn_model_name; // name of builtin model
    vec cnn_center;
    fl resolution; //this isn't specified in model file, so be careful about straying from default
    unsigned cnn_rotations; //do we want to score multiple orientations?
    double subgrid_dim;
    bool cnn_scoring; //if true, do cnn_scoring of final pose
    bool cnn_refinement;
    bool outputdx;
    bool outputxyz;
    bool gradient_check;
    bool move_minimize_frame;  //recenter with every scoring evaluation
    bool fix_receptor;
    bool verbose;
    std::string xyzprefix;
    unsigned seed; //random seed

    cnn_options()
        : cnn_model_name("default2017"), cnn_center(NAN, NAN, NAN), resolution(0.5), cnn_rotations(0),
            subgrid_dim(0.0), cnn_scoring(false), cnn_refinement(false), outputdx(false),
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
    fl energy_range;
    sz num_modes;
    fl out_min_rmsd;
    fl forcecap;
    int seed;
    int verbosity;
    int cpu;
    int device; //gpu number

    int exhaustiveness;
    int num_mc_steps;
    bool score_only;
    bool randomize_only;
    bool local_only;
    bool dominimize;
    bool include_atom_info;
    bool gpu_on;

    cnn_options cnnopts;

    //reasonable defaults
    user_settings()
        : energy_range(2.0), num_modes(9), out_min_rmsd(1), forcecap(1000),
            seed(auto_seed()), verbosity(1), cpu(1), device(0),
            exhaustiveness(10), num_mc_steps(0), score_only(false),
            randomize_only(false), local_only(false), dominimize(false),
            include_atom_info(false), gpu_on(false) {

    }
};

