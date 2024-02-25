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
#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <vector>

#include "dl_scorer.h"
#include "model.h"
#include "cnn_data.h"


/* This class evaluates protein-ligand poses according to a provided
 * Caffe convolutional neural net (CNN) model.
 */

class CNNScorer : public DLScorer {
  public:
    typedef float Dtype;
  private:
    std::vector<caffe::shared_ptr<caffe::Net<Dtype> > > nets;
    std::vector<caffe::MolGridDataLayer<Dtype> *> mgrids;
    std::vector<caffe::MolGridDataParameter *> mgridparams;

    //scratch vectors to avoid memory reallocation
    std::vector<gfloat3> gradient, gradient_rec, gradient_lig;
    std::vector<gfloat3> atoms;
    std::vector<short> channels;

    void getGradient(caffe::MolGridDataLayer<Dtype> *mgrid);

  public:
    CNNScorer() {
    }
    virtual ~CNNScorer() {
    }

    CNNScorer(const cnn_options& opts);

    bool initialized() const {
      return nets.size() > 0;
    }

    bool has_affinity() const; //return true if can predict affinity

    float score(model& m,float& variance); //score only - no gradient
    float score(model& m, bool compute_gradient, float& affinity, float& loss, float& variance);

    void outputDX(const std::string& prefix, double scale = 1.0, bool relevance =
        false, std::string layer_to_ignore = "", bool zero_values = false);
    void outputXYZ(const std::string& base, const std::vector<gfloat3>& atoms,
        const std::vector<short>& whichGrid, const std::vector<gfloat3>& gradient);
    std::unordered_map<std::string, float> get_gradient_norm_per_atom(bool receptor);
    std::unordered_map<std::string, float> get_relevance_per_atom(bool receptor);

    void lrp(const model& m, const std::string& layer_to_ignore = "",
        bool zero_values = false);
    void gradient_setup(const model& m, const std::string& recname,
        const std::string& ligname, const std::string& layer_to_ignore = "");

    const cnn_options& options() const {
      return cnnopts;
    }



    fl get_grid_dim() const {
      if(mgrids.size() == 0) throw usage_error("CNN network not initialized in get_grid_dim");
      return mgrids[0]->getDimension();
    }
    fl get_grid_res() const {
      if(mgrids.size() == 0) throw usage_error("CNN network not initialized in get_grid_res");
      return mgrids[0]->getResolution();
    }

    void set_bounding_box(grid_dims& box) const;

    caffe::MolGridDataLayer<Dtype> * get_mgrid(unsigned which=0) {
      if(which >= mgrids.size()) throw usage_error("CNN network doesn't exist");
      return mgrids[which];
    }

     std::shared_ptr<DLScorer> fresh_copy() const {
      return std::make_shared<CNNScorer>(cnnopts);
     }

  protected:
    void get_net_output(caffe::shared_ptr<caffe::Net<Dtype> >& net, Dtype& score, Dtype& aff, Dtype& loss);
    void check_gradient(caffe::shared_ptr<caffe::Net<Dtype> >& net);

};

#endif /* SRC_LIB_CNN_SCORER_H_ */
