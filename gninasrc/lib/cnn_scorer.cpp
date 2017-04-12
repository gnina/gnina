/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"
#include "gridoptions.h"

#include "caffe/proto/caffe.pb.h"
#include "caffe/layers/pooling_layer.hpp"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include "nngridder.h"

#include "cnn_data.h"

using namespace caffe;
using namespace std;

//initialize from commandline options
//throw error if missing required info
CNNScorer::CNNScorer(const cnn_options& cnnopts, const vec& center,
        const model& m) :
        rotations(cnnopts.cnn_rotations), seed(cnnopts.seed),
        outputdx(cnnopts.outputdx), outputxyz(cnnopts.outputxyz),
        mtx(new boost::mutex) {

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

        net.reset(new Net<Dtype>(param));

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


//returns relevance scores per atom
std::vector<float> CNNScorer::get_relevances(bool receptor)
{
    vector<float4> atoms;

    if (receptor)
    {
        mgrid->getReceptorAtoms(0, atoms);
        mgrid->getReceptorGradient(0,gradient);
    }
    else
    {
        mgrid->getLigandAtoms(0, atoms);
        mgrid->getLigandGradient(0, gradient);
    }

    std::vector<float> relevances(atoms.size());

    for (unsigned i = 0, n = gradient.size(); i < n; ++i)
    {
        //relevances summed to x field of gradient, will be fixed
        relevances[i] = gradient[i].x;
        std::cout << "OUTPUT GRADIENT: " << gradient[i].x << '\n';
    }

    return relevances;
}


void CNNScorer::lrp(const model& m, const string& recname, const string& ligname)
{
    boost::lock_guard<boost::mutex> guard(*mtx);
    
    caffe::Caffe::set_random_seed(seed); //same random rotations for each ligand..

    mgrid->setReceptor<atom>(m.get_fixed_atoms());
    mgrid->setLigand<atom,vec>(m.get_movable_atoms(),m.coordinates());

    net->Forward();
    net->Backward_relevance();
    
    //outputXYZ("LRP_rec", mgrid->getReceptorAtoms(0), mgrid->getReceptorChannels(0), mgrid->getReceptorGradient(0));
	mgrid->getLigandGradient(0, gradient);
	mgrid->getLigandAtoms(0, atoms);
	mgrid->getLigandChannels(0, channels);
    outputXYZ(ligname, atoms, channels, gradient);

	mgrid->getReceptorGradient(0, gradient);
	mgrid->getReceptorAtoms(0, atoms);
	mgrid->getReceptorChannels(0, channels);
    outputXYZ(recname, atoms, channels, gradient);
}

//has an affinity prediction layer
bool CNNScorer::has_affinity() const
{
	return (bool)net->blob_by_name("predaff");
}

//return score of model, assumes receptor has not changed from initialization
//if compute_gradient is set, also adds cnn atom gradient to m.minus_forces
float CNNScorer::score(model& m, bool compute_gradient, float& aff)
{
	boost::lock_guard<boost::mutex> guard(*mtx);
	if (!initialized())
		return -1.0;

	caffe::Caffe::set_random_seed(seed); //same random rotations for each ligand..

	mgrid->setReceptor<atom>(m.get_fixed_atoms());
	mgrid->setLigand<atom,vec>(m.get_movable_atoms(), m.coordinates());

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
		if (affblob)
		{
			//has affinity prediction
			const Dtype* aff = affblob->cpu_data();
			affinity += aff[0];
			cout << "#Rotate " << out[1] << " " << aff[0] << "\n";
		}
		else
		{
			cout << "#Rotate " << out[1] << "\n";
		}
		if (compute_gradient)
		{
			net->Backward();
			mgrid->getLigandGradient(0, gradient);
			m.add_minus_forces(gradient); //TODO divide by cnt?
		}
		cnt++;
	}

	if (outputdx) {
		outputDX(m.get_name());
	}
	if (outputxyz) {
		if (!compute_gradient)
			net->Backward();
		const string& ligname = m.get_name() + "_lig";
		const string& recname = m.get_name() + "_rec";
		mgrid->getLigandGradient(0, gradient);
		mgrid->getLigandAtoms(0, atoms);
		mgrid->getLigandChannels(0, channels);
		outputXYZ(ligname, atoms, channels, gradient);

		mgrid->getReceptorGradient(0, gradient);
		mgrid->getReceptorAtoms(0, atoms);
		mgrid->getReceptorChannels(0, channels);
		outputXYZ(recname, atoms, channels, gradient);
	}

	//TODO m.scale_minus_forces(1 / cnt);
	aff = affinity / cnt;
	return score / cnt;
}


//return only score
float CNNScorer::score(model& m)
{
	float aff = 0;
	return score(m, false, aff);
}


//dump dx files of the diff
void CNNScorer::outputDX(const string& prefix, double scale, bool relevance)
{
	const caffe::shared_ptr<Blob<Dtype> > datablob = net->blob_by_name("data");
	const vector<caffe::shared_ptr<Layer<Dtype> > >& layers = net->layers();
	if(datablob) {
		//this is a big more fragile than I would like.. if there is a pooling layer before
		//the first convoluational of fully connected layer and it is a max pooling layer,
		//change it to average before the backward to avoid a discontinuous map
		PoolingLayer<Dtype> *pool = NULL;
		for(unsigned i = 1, nl = layers.size(); i < nl; i++)
		{
			pool = dynamic_cast<PoolingLayer<Dtype>*>(layers[i].get());
			if(pool)
				break; //found it
			else if(layers[i]->type() == string("Convolution"))
				break; //give up
			else if(layers[i]->type() == string("InnerProduct"))
				break;
		}
		if(pool) {
			if(pool->pool() == PoolingParameter_PoolMethod_MAX) {
				pool->set_pool(PoolingParameter_PoolMethod_AVE);
			} else {
				pool = NULL; //no need to reset to max
			}
		}

		//must redo backwards with average pooling
		if(relevance)
			net->Backward_relevance();
		else
			net->Backward();

		string p = prefix;
		if(p.length() == 0) p = "dx";
		mgrid->dumpDiffDX(p, datablob.get(), scale);

		if(pool) {
			pool->set_pool(PoolingParameter_PoolMethod_MAX);
		}

	}
}


void CNNScorer::outputXYZ(const string& base, const vector<float4>& atoms, const vector<short>& whichGrid, const vector<float3>& gradient)
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

