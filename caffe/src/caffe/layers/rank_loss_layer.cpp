#include <algorithm>
#include <cmath>
#include <vector>

#include "caffe/layers/rank_loss_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template<typename Dtype>
void RankLossLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
		const vector<Blob<Dtype>*>& top) {
	LossLayer<Dtype>::Reshape(bottom, top);

	//check single labels
	const vector<int>& labels = bottom[1]->shape();
	for (unsigned i = 1, n = labels.size(); i < n; i++) {
		//skip batch size (i = 0)
		CHECK_EQ(labels[i], 1);
	}

	//must have only one or two classes
	const vector<int>& preds = bottom[0]->shape();
	int numclasses = 1;
	for (unsigned i = 1, n = preds.size(); i < n; i++) {
		//skip batch size (i = 0)
		numclasses *= preds[i];
	}
	CHECK_LE(numclasses,2)<< "RankLossLayer requires one or two classes";
}

//compute pairwise rank loss
template<typename Dtype>
void RankLossLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
		const vector<Blob<Dtype>*>& top) {

	int num = bottom[0]->num();
	Dtype loss = 0;

	if (allpairs) { // all n^2/2 pairs
		unsigned npairs = 0;
		for (unsigned i = 0; i < num; i++) {
			for (unsigned j = i + 1; j < num; j++) {
				loss += compute_pair_loss(bottom, i, j);
				npairs += 1;
			}
		}
		loss /= npairs;
	} else {
		//adjacent pairs only
		for (int i = 0; i < num; i += 2) {
			loss += compute_pair_loss(bottom, i, i + 1);
		}
		loss /= (num / 2);
	}

	top[0]->mutable_cpu_data()[0] = loss;
}

template<typename Dtype>
void RankLossLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
		const vector<bool>& propagate_down,
		const vector<Blob<Dtype>*>& bottom) {
	if (propagate_down[1]) {
		LOG(FATAL)<< this->type()
		<< " Layer cannot backpropagate to label inputs.";
	}
	if (propagate_down[0]) {
	    Dtype* bottom_diff = bottom[0]->mutable_cpu_diff();
	    caffe_set(bottom[0]->count(), Dtype(0), bottom_diff);

		int num = bottom[0]->num();

		if(allpairs) { // all n^2/2 pairs
			for(unsigned i = 0; i < num; i++) {
				for(unsigned j = i+1; j < num; j++) {
					compute_pair_gradient(top, bottom, i, j);
				}
			}
			/*LOG(INFO) <<"RANKGRADIENTS";
			for(unsigned i = 0; i < num; i++) {
				 LOG(INFO) << i << ": " << bottom_diff[i] << "\n";
			}*/
		} else {
			//adjacent pairs only
			for (int i = 0; i < num; i+=2) {
				compute_pair_gradient(top,bottom, i, i+1);
			}
		}
	}
}

template<typename Dtype>
Dtype RankLossLayer<Dtype>::compute_pair_loss(
		const vector<Blob<Dtype>*>& bottom, unsigned i, unsigned j) {
	//bottom[0] should be scores
	//bottom[1] should be labels
	//scores can have one or two value (in which case binary classification is assumed)
	const Dtype* bottom_data = bottom[0]->cpu_data();
	const Dtype* bottom_label = bottom[1]->cpu_data();
	int num = bottom[0]->num();
	int dim = bottom[0]->count() / bottom[0]->num();
	Dtype loss = 0;
	Dtype si = 0, sj = 0;

	CHECK_LT(i, num)<< "RankLoss: invalid index i";
	CHECK_LT(j, num)<< "RankLoss: invalid index j";

	if (dim == 1) {
		//single scores
		si = bottom_data[i];
		sj = bottom_data[j];
	} else if (dim == 2) {
		//binary classification, take score of second class
		si = bottom_data[i*dim+1];
		sj = bottom_data[j*dim+1];
	} else {
		CHECK(false) << "RankLoss requires single score per example or binary class probabilities.";
	}

	Dtype Li = bottom_label[i];
	Dtype Lj = bottom_label[j];
	Dtype diff = si - sj;
	Dtype ediff = exp(diff);
	Dtype Pij =  1.0;
	if(std::isfinite(ediff)) 
		Pij = ediff / (1 + ediff);

	if (Li > Lj) {
		Pij = std::max(Pij, Dtype(kLOG_THRESHOLD));
		loss = -log(Pij);
	} else {
		Dtype val = 1-Pij;
		val = std::max(val, Dtype(kLOG_THRESHOLD));
		loss = -log(val);
	}
//LOG(INFO) << "RANKPAIRLOSS " << i <<","<<j<<" "<<loss<< " " << Pij << " " << ediff << " s: " << si << "," << sj << " L: " << Li << "," << Lj;
	return loss;
}

template<typename Dtype>
void RankLossLayer<Dtype>::compute_pair_gradient(
		const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom,
		unsigned i, unsigned j) {
	const Dtype* bottom_data = bottom[0]->cpu_data();
	const Dtype* bottom_label = bottom[1]->cpu_data();
    Dtype* bottom_diff = bottom[0]->mutable_cpu_diff();

	int num = bottom[0]->num();
	int dim = bottom[0]->count() / bottom[0]->num();
	Dtype si = 0, sj = 0;

	CHECK_LT(i, num)<< "RankLoss: invalid index i";
	CHECK_LT(j, num)<< "RankLoss: invalid index j";

	if (dim == 1) {
		//single scores
		si = bottom_data[i];
		sj = bottom_data[j];
	} else if (dim == 2) {
		//binary classification, take score of second class
		si = bottom_data[i*dim+1];
		sj = bottom_data[j*dim+1];
	} else {
		CHECK(false) << "RankLoss requires single score per example or binary class probabilities.";
	}

	Dtype Li = bottom_label[i];
	Dtype Lj = bottom_label[j];
	Dtype ediff = exp(si-sj);
	Dtype d = 0;
	if(std::isfinite(ediff))
		d = - 1.0 / (1.0+exp(si-sj));

	if(Li < Lj) {
		//reverse sign
		d = -d;
	}

	//scale by num
	d /= num;
	//also by total loss
	Dtype scale = top[0]->cpu_diff()[0];
	d *= scale; //not

	if(dim == 1) {
		bottom_diff[i] += d;
		bottom_diff[j] += -d;
		//LOG(INFO) << "RANKLOSS1 " << i << " L:" << Li << " s:" << si << " "<< d << "\t" << j << " L:" << Lj << " s:" << sj << " " << -d << "\n";

	} else if(dim == 2) {
		//compensate for applying gradient to both classes
		d /= 2.0;
		bottom_diff[i*2] += -d;
		bottom_diff[i*2+1] += d;
		bottom_diff[j*2] += d;
		bottom_diff[j*2+1] += -d;

		//LOG(INFO) << "RANKLOSS " << i << " L:" << Li << " s:" << si << " "<< d << "\t" << j << " L:" << Lj << " s:" << sj << " " << -d << "\n";
		//LOG(INFO) << "RANKINPUTS " << i << "(" << Li <<"): " << bottom_data[i*dim] << "," << bottom_data[i*dim+1] << "   " << j << "(" << Lj <<"): " << bottom_data[j*dim] << "," << bottom_data[j*dim+1] << "\n";
	}
}

INSTANTIATE_CLASS(RankLossLayer);
REGISTER_LAYER_CLASS(RankLoss);

}  // namespace caffe
