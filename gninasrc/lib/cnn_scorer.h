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
#include "cnn_data.h"

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
    CNNScorer()
        : mgrid(NULL), mtx(new boost::mutex) {
    }
    virtual ~CNNScorer() {
    }

    CNNScorer(const cnn_options& opts, const model& m);

    bool initialized() const {
      return net.get();
    }

    bool has_affinity() const; //return true if can predict affinity

    float score(model& m); //score only - no gradient
    float score(model& m, bool compute_gradient, float& affinity, float& loss);

    void outputDX(const string& prefix, double scale = 1.0, bool relevance =
        false, string layer_to_ignore = "", bool zero_values = false);
    void outputXYZ(const string& base, const vector<float4>& atoms,
        const vector<short>& whichGrid, const vector<float3>& gradient);
    std::unordered_map<string, float> get_scores_per_atom(bool receptor,
        bool relevance = false);

    void lrp(const model& m, const string& layer_to_ignore = "",
        bool zero_values = false);
    void gradient_setup(const model& m, const string& recname,
        const string& ligname, const string& layer_to_ignore = "");

    //readjust center
    bool set_center_from_model(model &m);

    const cnn_options& options() const {
      return cnnopts;
    }

    vec get_center() const {
      return mgrid->getCenter();
    }
    fl get_grid_dim() const {
      return mgrid->getDimension();
    }
    fl get_grid_res() const {
      return mgrid->getResolution();
    }

  protected:
    void get_net_output(Dtype& score, Dtype& aff, Dtype& loss);
    void check_gradient();
    friend void test_set_atom_gradients();
};

#endif /* SRC_LIB_CNN_SCORER_H_ */
