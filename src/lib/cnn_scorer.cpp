/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"
#include "gridoptions.h"

#include "caffe/layer.hpp"
#include "caffe/net.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/layers/ndim_data_layer.hpp"

using namespace caffe;
using namespace std;

//initialize from commandline options
//throw error if missing required info
CNNScorer::CNNScorer(const cnn_options& cnnopts, const vec& center,
		const model& m) :
		rotations(cnnopts.cnn_rotations)
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
		NDimDataParameter *ndimparam = first->mutable_ndim_data_param();
		if (ndimparam == NULL)
		{
			throw usage_error("First layer of model must be NDimData.");
		}
		ndimparam->set_inmemory(true);

		//set batch size to 1
		unsigned bsize = 1;
		//unless we have rotations, in which case do them all at once, which turns out isn't actually faster :-(
		if (cnnopts.cnn_rotations > 0)
		{
			//let user specify rotations
			unsigned nrot = min(24U, cnnopts.cnn_rotations);
			ndimparam->set_rotate(nrot);
			//BUT it turns out this isn't actually faster
			//bsize = nrot;
		}

		ndimparam->set_batch_size(bsize);

		net.reset(new Net<float>(param));

		//load weights
		net->CopyTrainedLayersFrom(cnnopts.cnn_weights);

		//check that network matches our expectations

		//the first layer must be NDimData
		const vector<shared_ptr<Layer<Dtype> > >& layers = net->layers();

		//we also need an output layer
		if (layers.size() < 1)
		{
			throw usage_error("No layers in model!");
		}

		ndim = dynamic_cast<NDimDataLayer<Dtype>*>(layers[0].get());
		if (ndim == NULL)
		{
			throw usage_error("First layer of model must be NDimData.");
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

		//initialize receptor part of grid
		BlobShape shape = ndim->layer_param().ndim_data_param().shape();
		if (shape.dim_size() != 4)
		{
			throw usage_error(
					"Input data layer does not have correct number of dimensions.");
		}
		unsigned nchannels = shape.dim(0);
		unsigned dim = shape.dim(1);
		if (dim != shape.dim(2) || dim != shape.dim(3))
		{
			throw usage_error(
					"Input data layer does not have cubic dimensions.");
		}

		gridoptions gopt;
		gopt.res = cnnopts.resolution;
		gopt.dim = round((dim - 1) * gopt.res);

		grid.initialize(gopt);
	}

}

//return score of model, assumes receptor has not changed from initialization
float CNNScorer::score(const model& m)
{
	if (!initialized())
		return -1.0;
	grid.setModel(m);

	grid.outputMem(ndim->getMemoryData());
	ndim->memoryIsSet();

	double score = 0.0;
	const shared_ptr<Blob<Dtype> > outblob = net->blob_by_name("output");

	unsigned num = 1;
	if (rotations > 0)
	{
		num = (rotations / (outblob->count() / 2));
	}

	unsigned cnt = 0;
	for (unsigned r = 0; r < num; r++)
	{
		net->Forward(); //do all rotations at once if requested

		//take average of all rotations if requested
		for (unsigned i = 1, n = outblob->count(); i < n; i += 2)
		{
			const Dtype* out = outblob->cpu_data();
			score += out[i];
			cout << "#Rotate " << out[i] << "\n";
			cnt++;
		}
	}

	return score / cnt;

}
