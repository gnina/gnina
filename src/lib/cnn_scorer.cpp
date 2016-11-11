/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"
#include "gridoptions.h"
#include "nngridder.h"

#include "boost/multi_array.hpp"
#include "caffe/layer.hpp"
#include "caffe/net.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/layers/ndim_data_layer.hpp"
#include "caffe/layers/conv_layer.hpp"

using namespace caffe;
using namespace std;

//initialize from commandline options
//throw error if missing required info
CNNScorer::CNNScorer(const cnn_options& cnnopts, const vec& center,
		const model& m) :
		rotations(cnnopts.cnn_rotations),
		 mtx(new boost::mutex)
{

	if (cnnopts.cnn_scoring)
	{
		if (cnnopts.cnn_model.size() == 0)
		{
			throw usage_error("Missing model for cnn scoring.");
		}
		if (cnnopts.cnn_weights.size() == 0)
		{
			throw usage_error("Missing weights for cnn scoring.");
		}

		//load cnn model
		NetParameter param;
		ReadNetParamsFromTextFileOrDie(cnnopts.cnn_model, &param);
		param.mutable_state()->set_phase(TEST);
		LayerParameter *first = param.mutable_layer(0);
		//must be ndim
		MolGridDataParameter *mgridparam = first->mutable_molgrid_data_param();
		if (mgridparam == NULL)
		{
			throw usage_error("First layer of model must be MolGridData.");
		}
		mgridparam->set_inmemory(true);

		//set batch size to 1
		unsigned bsize = 1;
		//unless we have rotations, in which case do them all at once, which turns out isn't actually faster :-(
		if (cnnopts.cnn_rotations > 0)
		{
			//let user specify rotations
			unsigned nrot = min(24U, cnnopts.cnn_rotations);
			mgridparam->set_rotate(nrot);
			//BUT it turns out this isn't actually faster
			//bsize = nrot;
		}
		if (cnnopts.cnn_gradient)
		{
			param.set_force_backward(true);
		}

		net.reset(new Net<float>(param));

		//load weights
		net->CopyTrainedLayersFrom(cnnopts.cnn_weights);

		//check that network matches our expectations

		//the first layer must be MolGridLayer
		const vector<caffe::shared_ptr<Layer<Dtype> > >& layers = net->layers();

		//we also need an output layer
		if (layers.size() < 1)
		{
			throw usage_error("No layers in model!");
		}

		mgrid = dynamic_cast<MolGridDataLayer<Dtype>*>(layers[0].get());
		if (mgrid == NULL)
		{
			throw usage_error("First layer of model must be MolGridDataLayer.");
		}

		if (!net->has_blob("output"))
		{
			throw usage_error("Model must have output layer named \"output\".");
		}
		if (net->blob_by_name("output")->count() != 2 * bsize)
		{
			throw usage_error(
					"Model output layer does not have exactly two outputs.");
		}
	}

}

//return score of model, assumes receptor has not changed from initialization
float CNNScorer::score(const model& m)
{
	boost::lock_guard<boost::mutex> guard(*mtx);
	if (!initialized())
		return -1.0;

	mgrid->setReceptor<atom>(m.get_fixed_atoms());
	mgrid->setLigand<atom,vec>(m.get_movable_atoms(),m.coordinates());

	double score = 0.0;
	const caffe::shared_ptr<Blob<Dtype> > outblob = net->blob_by_name("output");
	const caffe::shared_ptr<Blob<Dtype> > inblob = net->blob_by_name("data");

	unsigned cnt = 0;
	for (unsigned r = 0, n = max(rotations, 1U); r < n; r++)
	{
		net->Forward(); //do all rotations at once if requested

		caffe::shared_ptr<Blob<Dtype> > pooled = net->blob_by_name("unit1_pool");
		if (true) //TODO get from cnn_opts
		{
			int c = 34, s = 24; //TODO get these from params
			net->Backward();
			boost::multi_array<float, 4> gradient(boost::extents[c][s][s][s]);
			memcpy(gradient.origin(), pooled->cpu_diff(), c*s*s*s*sizeof(float));
			NNGridder::outputDX("DEBUG", gradient, m);
		}

		const Dtype* out = outblob->cpu_data();
		score += out[1];
		cout << "#Rotate " << out[1] << "\n";
		cnt++;
	}

	return score / cnt;
}
