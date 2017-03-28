/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"
#include "gridoptions.h"

#include "caffe/proto/caffe.pb.h"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include "cnn_data.h"

using namespace caffe;
using namespace std;

//initialize from commandline options
//throw error if missing required info
CNNScorer::CNNScorer(const cnn_options& cnnopts, const vec& center,
		const model& m) :
		rotations(cnnopts.cnn_rotations), seed(cnnopts.seed),
		 outputdx(cnnopts.outputdx), mtx(new boost::mutex)
{

	if (cnnopts.cnn_scoring)
	{
		NetParameter param;

		//load cnn model
		if (cnnopts.cnn_model.size() == 0)
		{
			google::protobuf::io::ArrayInputStream  modeldata(cnn_default_model, strlen(cnn_default_model));
			bool success = google::protobuf::TextFormat::Parse(&modeldata, &param);
			if(!success)
				throw usage_error("Error with default cnn model.");
			UpgradeNetAsNeeded("default",&param);
		}
		else
		{
			ReadNetParamsFromTextFileOrDie(cnnopts.cnn_model, &param);
		}

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
			unsigned nrot = cnnopts.cnn_rotations;
			mgridparam->set_random_rotation(true);
			//I think there's a bug in the axial rotations - they aren't all distinct
			//BUT it turns out this isn't actually faster
			//bsize = nrot;
		}
		param.set_force_backward(true);

		net.reset(new Net<float>(param));

		//load weights
		if (cnnopts.cnn_weights.size() == 0)
		{
			NetParameter wparam;
			google::protobuf::io::ArrayInputStream  weightdata(cnn_default_weights, cnn_default_weights_len);
			google::protobuf::io::CodedInputStream strm(&weightdata);
			strm.SetTotalBytesLimit(INT_MAX, 536870912);
			bool success = wparam.ParseFromCodedStream(&strm);
			if(!success)
				throw usage_error("Error with default weights.");

			net->CopyTrainedLayersFrom(wparam);
		}
		else
		{
			net->CopyTrainedLayersFrom(cnnopts.cnn_weights);
		}

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

//has an affinity prediction layer
bool CNNScorer::has_affinity() const
{
	return (bool)net->blob_by_name("predaff");
}

//return score and affinity of model, assumes receptor has not changed from initialization
float CNNScorer::score(const model& m, float& aff)
{
	boost::lock_guard<boost::mutex> guard(*mtx);
	if (!initialized())
		return -1.0;
	
	caffe::Caffe::set_random_seed(seed); //same random rotations for each ligand..

	mgrid->setReceptor<atom>(m.get_fixed_atoms());
	mgrid->setLigand<atom,vec>(m.get_movable_atoms(),m.coordinates());

	double score = 0.0;
	double affinity = 0.0;
	const caffe::shared_ptr<Blob<Dtype> > outblob = net->blob_by_name("output");
	const caffe::shared_ptr<Blob<Dtype> > affblob = net->blob_by_name("predaff");
	unsigned cnt = 0;
	for (unsigned r = 0, n = max(rotations, 1U); r < n; r++)
	{
		net->Forward(); //do all rotations at once if requested

		const Dtype* out = outblob->cpu_data();
		score += out[1];

		if(affblob)
		{
			//has affinity prediction
			const Dtype* aff = affblob->cpu_data();
			affinity += aff[0];
			//TODO: use the log
			cout << "#Rotate " << out[1] << " " << aff[0] << "\n";
		}
		else
		{
			cout << "#Rotate " << out[1] << "\n";
		}
		cnt++;

	}

	if(outputdx) {
		const caffe::shared_ptr<Blob<Dtype> > datablob = net->blob_by_name("data");
		if(datablob) {
			net->Backward();
			mgrid->dumpDiffDX(m.get_name(), datablob.get());
		}
	}
	aff = affinity/cnt;
	return score / cnt;

}


//return only score
float CNNScorer::score(const model& m)
{
	float aff = 0;
	return score(m,aff);
}
