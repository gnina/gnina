/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"

#include "caffe/common.hpp"
#include "caffe/layer.hpp"
#include "caffe/net.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/layers/ndim_data_layer.hpp"

using namespace caffe;
using namespace std;

//initialize from commandline options
//throw error if missing required info
CNNScorer::CNNScorer(const cnn_options& cnnopts, const model& m) {
	if(cnnopts.cnn_scoring) {
		//load cnn model
	  NetParameter param;
	  ReadNetParamsFromTextFileOrDie(cnnopts.cnn_model, &param);
	  param.mutable_state()->set_phase(TEST);
	  net.reset(new Net<float>(param));

		//load weights
	  net->CopyTrainedLayersFrom(cnnopts.cnn_weights);

	  //check that network matches our expectations

	  //the first layer must be NDimData
	  const vector<shared_ptr<Layer<Dtype> > >& layers = net->layers();

	  //we also need an output layer
	  if(layers.size() < 1) {
	  	throw usage_error("No layers in model!");
	  }

	  NDimDataLayer<Dtype> *data = dynamic_cast<NDimDataLayer<Dtype>*>(layers[0].get());
	  if(data == NULL) {
	  	throw usage_error("First layer of model must be NDimData.");
	  }

	  if(net->num_outputs() > 1) {
	  	throw usage_error("Model must produce single output layer.");
	  }

		//initialize receptor part of grid
	}

}
