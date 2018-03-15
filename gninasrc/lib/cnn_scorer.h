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

    // caffe::shared_ptr<boost::mutex> mtx; //todo, enable parallel scoring

    //scratch vectors to avoid memory reallocation
    vector<float3> gradient;
    vector<float4> atoms;
    vector<short> channels;

public:
    CNNScorer(): mgrid(NULL), rotations(0), outputdx(false), outputxyz(false), gradient_check(false), reset_center(true)//, mtx(new boost::mutex) 
    {}
    virtual ~CNNScorer() {}

    CNNScorer(const cnn_options& cnnopts, const vec& center, const model& m);

    CNNScorer(const CNNScorer& src, caffe::MolGridDataLayer<Dtype>* src_grid) : 
        mgrid(src_grid), mgridparam(src.mgridparam), rotations(src.rotations), 
        seed(src.seed), outputdx(src.outputdx), outputxyz(src.outputxyz), 
        gradient_check(src.gradient_check), reset_center(src.reset_center), 
        xyzprefix(src.xyzprefix), gradient(src.gradient), atoms(src.atoms), 
        channels(src.channels) {
    }

    bool initialized() const { return net.get(); }

    bool has_affinity() const; //return true if can predict affinity

    float score(model& m, bool silent = true);
    float score(model& m, bool compute_gradient, float& affinity, bool silent = true);

    void outputDX(const string& prefix, double scale = 1.0, bool relevance = false, string layer_to_ignore = "", bool zero_values = false);
    void outputXYZ(const string& base, const vector<float4>& atoms,
               const vector<short>& whichGrid, const vector<float3>& gradient);
    std::unordered_map<string, float> get_scores_per_atom(bool receptor, bool relevance = false);

    void lrp(const model& m, const string& layer_to_ignore = "", bool zero_values = false);
    void gradient_setup(const model& m, const string& recname, const string& ligname, const string& layer_to_ignore = "");

    bool adjust_center() const;

protected:
  void get_net_output(Dtype& score, Dtype& aff);
  void check_gradient();

};

#endif /* SRC_LIB_CNN_SCORER_H_ */
