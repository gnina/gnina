/*
 * cnn_scorer.h
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#ifndef SRC_LIB_CNN_SCORER_H_
#define SRC_LIB_CNN_SCORER_H_

#include "caffe/caffe.hpp"
#include "caffe/net.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/molgrid_data_layer.hpp"
#include "boost/thread/mutex.hpp"

#include "nngridder.h"
#include "model.h"

struct cnn_options {
    //stores options associated with cnn scoring
    std::string cnn_model; //path to model file
    std::string cnn_weights; //weights for model
    std::string cnn_recmap; //optional file specifying receptor atom typing to channel map
    std::string cnn_ligmap; //optional file specifying ligand atom typing to channel map
    vec cnn_center;
    fl resolution; //this isn't specified in model file, so be careful about straying from default
    unsigned cnn_rotations; //do we want to score multiple orientations?
    bool cnn_scoring; //if true, do cnn_scoring of final pose
    bool outputdx;
    bool outputxyz;
    bool gradient_check;
    bool move_minimize_frame;  //recenter with every scoring evaluation
    bool fix_receptor;
    bool verbose;
    std::string xyzprefix;
    unsigned seed; //random seed

    cnn_options(): cnn_center(NAN,NAN,NAN), resolution(0.5), cnn_rotations(0), cnn_scoring(false), outputdx(false), outputxyz(false), gradient_check(false), move_minimize_frame(false), fix_receptor(true), verbose(false), seed(0) {}
};

/* This class evaluates protein-ligand poses according to a provided
 * Caffe convolutional neural net (CNN) model.
 */
class CNNScorer {
    typedef float Dtype;
    caffe::shared_ptr<caffe::Net<Dtype> > net;
    caffe::MolGridDataLayer<Dtype> *mgrid;
    caffe::MolGridDataParameter *mgridparam;
    cnn_options cnnopts;

    caffe::shared_ptr<boost::mutex> mtx; //todo, enable parallel scoring

    //scratch vectors to avoid memory reallocation
    vector<float3> gradient;
    vector<float4> atoms;
    vector<short> channels;

public:
    CNNScorer(): mgrid(NULL), mtx(new boost::mutex) {}
    virtual ~CNNScorer() {}

    CNNScorer(const cnn_options& opts, const model& m);

    bool initialized() const { return net.get(); }

    bool has_affinity() const; //return true if can predict affinity

    float score(model& m); //score only - no gradient
    float score(model& m, bool compute_gradient, float& affinity, float& loss);

    void outputDX(const string& prefix, double scale = 1.0, bool relevance = false, string layer_to_ignore = "", bool zero_values = false);
    void outputXYZ(const string& base, const vector<float4>& atoms,
               const vector<short>& whichGrid, const vector<float3>& gradient);
    std::unordered_map<string, float> get_scores_per_atom(bool receptor, bool relevance = false);

    void lrp(const model& m, const string& layer_to_ignore = "", bool zero_values = false);
    void gradient_setup(const model& m, const string& recname, const string& ligname, const string& layer_to_ignore = "");

    //readjust center
    bool set_center_from_model(model &m);

    const cnn_options& options() const { return cnnopts; }

    vec get_center() const { return mgrid->getCenter(); }
    fl get_grid_dim() const { return mgrid->getDimension(); }
    fl get_grid_res() const { return mgrid->getResolution(); }

protected:
  void get_net_output(Dtype& score, Dtype& aff, Dtype& loss);
  void check_gradient();
  friend void test_set_atom_gradients();
};

#endif /* SRC_LIB_CNN_SCORER_H_ */
