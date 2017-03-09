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
		 mtx(new boost::mutex)
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

		//backpropagate to atoms
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

//return score of model, assumes receptor has not changed from initialization
//if compute_gradient is set, also adds cnn atom gradient to m.minus_forces
float CNNScorer::score(model& m, bool compute_gradient)
{
	boost::lock_guard<boost::mutex> guard(*mtx);
	if (!initialized())
		return -1.0;

	caffe::Caffe::set_random_seed(seed); //same random rotations for each ligand..

	mgrid->setReceptor<atom>(m.get_fixed_atoms());
	mgrid->setLigand<atom,vec>(m.get_movable_atoms(), m.coordinates());

	double score = 0.0;
	const caffe::shared_ptr<Blob<Dtype> > outblob = net->blob_by_name("output");

	unsigned cnt = 0;
	for (unsigned r = 0, n = max(rotations, 1U); r < n; r++)
	{
		net->Forward();
		score += outblob->cpu_data()[1];
		if (compute_gradient)
		{
			net->Backward();
			m.add_minus_forces(mgrid->getLigandGradient(0)); //TODO divide by cnt
		}
		cnt++;
	}
	return score / cnt;
}

void CNNScorer::outputXYZ(const string& base, const vector<float4> atoms, const vector<short> whichGrid, const vector<float3> gradient)
{
	const char* sym[] = {"C", "C", "C", "C", "Ca", "Fe", "Mg", "N", "N", "N", "N", "O", "O", "P", "S", "Zn",
			     "C", "C", "C", "C", "Br", "Cl", "F",  "N", "N", "N", "N", "O", "O", "O", "P", "S", "S", "I"};

	const string& fname = base + ".xyz";
	ofstream out(fname.c_str());
	out.precision(5);

	out << atoms.size() << "\n\n";
	for (unsigned i = 0, n = atoms.size(); i < n; ++i)
	{
		out << sym[whichGrid[i]] << " ";
		out << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z << " ";
		out << gradient[i].x << " " << gradient[i].y << " " << gradient[i].z;
		if (i + 1 < n) out << "\n";
	}
}

