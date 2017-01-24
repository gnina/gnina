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
		if (cnnopts.cnn_gradient)
		{
			param.set_force_backward(true);
			compute_gradient = true;
		}

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
float CNNScorer::score(const model& m)
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
			outputXYZ("DEBUG", m, mgrid->getReceptorAtomGradient(0), mgrid->getLigandAtomGradient(0));
		}

		cnt++;
	}

	return score / cnt;
}

void CNNScorer::outputXYZ(const string& base, const model& m, const vector<float3> rec_grad, const vector<float3> lig_grad)
{
	const char* el_sym[] = {"H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Met"};

	const string& rec_fname = base + "_rec.xyz";
	ofstream rec_out(rec_fname.c_str());
	rec_out.precision(5);
	const atomv& rec_atoms = m.get_fixed_atoms();
	rec_out << rec_atoms.size() << "\n\n";
	for (unsigned i = 0, n = rec_atoms.size(); i < n; ++i)
	{
		rec_out << el_sym[smina_atom_type::default_data[rec_atoms[i].sm].el] << " ";
		rec_out << rec_atoms[i].coords[0] << " " << rec_atoms[i].coords[1] << " " << rec_atoms[i].coords[2] << " ";
		rec_out << rec_grad[i].x << " " << rec_grad[i].y << " " << rec_grad[i].z;
		if (i + 1 < n)
			rec_out << "\n";
	}

	const string& lig_fname = base + "_lig.xyz";
	ofstream lig_out(lig_fname.c_str());
	lig_out.precision(5);
	const atomv& lig_atoms = m.get_movable_atoms();
	const vector<vec> lig_coords = m.coordinates();
	lig_out << lig_atoms.size() << "\n\n";
	for (unsigned i = 0, n = lig_atoms.size(); i < n; ++i)
	{
		lig_out << el_sym[smina_atom_type::default_data[lig_atoms[i].sm].el] << " ";
		lig_out << lig_coords[i][0] << " " << lig_coords[i][1] << " " << lig_coords[i][2] << " ";
		lig_out << lig_grad[i].x << " " << lig_grad[i].y << " " << lig_grad[i].z;
		if (i + 1 < n)
			lig_out << "\n";
	}
}

