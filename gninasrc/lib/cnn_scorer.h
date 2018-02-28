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
    fl resolution; //this isn't specified in model file, so be careful about straying from default
    unsigned cnn_rotations; //do we want to score multiple orientations?
    bool cnn_scoring; //if true, do cnn_scoring of final pose
    bool outputdx;
    bool outputxyz;
    bool gradient_check;
    std::string xyzprefix;
    unsigned seed; //random seed

    cnn_options(): resolution(0.5), cnn_rotations(0), cnn_scoring(false), outputdx(false), outputxyz(false), gradient_check(false), seed(0) {}
};

/* This class evaluates protein-ligand poses according to a provided
 * Caffe convolutional neural net (CNN) model.
 */
class CNNScorer {
    typedef float Dtype;
    caffe::shared_ptr<caffe::Net<Dtype> > net;
    caffe::MolGridDataLayer<Dtype> *mgrid;
    caffe::MolGridDataParameter *mgridparam;
    unsigned rotations;
    unsigned seed;
    bool outputdx;
    bool outputxyz;
    bool gradient_check;
    mutable bool reset_center; //potential hack for debugging gradients
    std::string xyzprefix;

    caffe::shared_ptr<boost::mutex> mtx; //todo, enable parallel scoring

    //scratch vectors to avoid memory reallocation
    vector<float3> gradient;
    vector<float4> atoms;
    vector<short> channels;

public:
    CNNScorer(): mgrid(NULL), rotations(0), outputdx(false), outputxyz(false), gradient_check(false), reset_center(true), mtx(new boost::mutex) {}
    virtual ~CNNScorer() {}

    CNNScorer(const cnn_options& cnnopts, const vec& center, const model& m);

    bool initialized() const { return net.get(); }

    bool has_affinity() const; //return true if can predict affinity

    float score(model& m); //score only - no gradient
   // float score(model& m, float& affinity); //scores only - no gradient
    float score(model& m, bool compute_gradient, float& affinity, float& loss, bool silent = true);

    void outputDX(const string& prefix, double scale = 1.0, bool relevance = false, string layer_to_ignore = "", bool zero_values = false);
    void outputXYZ(const string& base, const vector<float4>& atoms,
               const vector<short>& whichGrid, const vector<float3>& gradient);
    std::unordered_map<string, float> get_scores_per_atom(bool receptor, bool relevance = false);

    void lrp(const model& m, const string& layer_to_ignore = "", bool zero_values = false);
    void gradient_setup(const model& m, const string& recname, const string& ligname, const string& layer_to_ignore = "");

    bool adjust_center() const;

protected:
  void get_net_output(Dtype& score, Dtype& aff, Dtype& loss);
  void check_gradient();

};

#endif /* SRC_LIB_CNN_SCORER_H_ */
