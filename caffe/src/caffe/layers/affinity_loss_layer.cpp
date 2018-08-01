#include <vector>

#include "caffe/layers/affinity_loss_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template<typename Dtype>
void AffinityLossLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  LossLayer<Dtype>::Reshape(bottom, top);
  CHECK_EQ(bottom[0]->count(1), bottom[1]->count(1))<< "Inputs must have the same dimension.";
  diff_.ReshapeLike(*bottom[0]);
}

template<typename Dtype>
void AffinityLossLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  int count = bottom[0]->count();
  Dtype sum = 0.0;
  Dtype gap = this->layer_param_.affinity_loss_param().gap() / 2.0;
  Dtype delta = this->layer_param_.affinity_loss_param().delta();
  Dtype penalty = this->layer_param_.affinity_loss_param().penalty();
  bool huber = this->layer_param_.affinity_loss_param().pseudohuber();
  Dtype ranklossm = this->layer_param_.affinity_loss_param().ranklossmult();
  bool fractional_gap = this->layer_param_.affinity_loss_param().fractional_gap();
  Dtype defaultzero = this->layer_param_.affinity_loss_param().diff_for_zero();
  bool hasweights = this->layer_param_.affinity_loss_param().weighted();

  Dtype delta2 = delta*delta;
  const Dtype *labels = bottom[1]->cpu_data();
  const Dtype *preds = bottom[0]->cpu_data();
  const Dtype *weights = NULL;
  if(hasweights) weights = bottom[2]->cpu_data();

  Dtype *d = diff_.mutable_cpu_data();

  for (unsigned i = 0; i < count; i++) {
    Dtype label = labels[i];
    Dtype pred = preds[i];
    Dtype diff = 0.0;
    Dtype weight = hasweights ? weights[i] : 1.0;
    if (label > 0) { //normal euclidean
      diff = pred - label;
      if(fractional_gap && gap > 0) {
        gap = fabs(diff)/gap;
      }
      if (diff < 0) {
        diff = std::min(diff + gap, Dtype(0));
      } else {
        diff = std::max(diff - gap, Dtype(0));
      }
    } else if (label < 0 && pred > -label) { //hinge like
      diff = pred + label + penalty;
      if(fractional_gap && gap > 0) {
        gap = fabs(diff)/gap;
      }
      diff = std::max(diff - gap, Dtype(0));
      
    } else { //ignore
      diff = defaultzero;
    }

    d[i] = weight*diff;

    if(huber) {
      //https://en.wikipedia.org/wiki/Huber_loss
      Dtype hval = diff/delta;
      sum += weight*delta2*(sqrt(1+hval*hval) - 1.0);
    }
    else {
      sum += weight*diff * diff;
    }

  }

  Dtype rankloss = 0;
  if(ranklossm > 0) {
    bool includeneg = this->layer_param_.affinity_loss_param().ranklossneg();
    //compute a rank loss
    nranklosspairs = 0; //save for backwards
    int num = bottom[0]->num();
    for (unsigned i = 0; i < num; i++) {
      Dtype labeli = labels[i];
      if(labeli > 0) {
        for (unsigned j = i + 1; j < num; j++) {
          Dtype labelj = labels[j];
          //correctly rank good poses
          if(labelj > 0) {
            rankloss += compute_pair_loss(bottom, i, j, ranklossm);
            nranklosspairs++;
          }
        }
      } else if(includeneg) { //negative, positive should be before
        if(i > 0 && labels[i-1] > 0 && labels[i-1] == -labels[i]) {
          rankloss += compute_pair_loss(bottom, i, i-1, ranklossm);
          nranklosspairs++;
        }
      }
    }
    rankloss /= nranklosspairs;
  }

  Dtype loss = sum / bottom[0]->num() / Dtype(2) + rankloss;
  top[0]->mutable_cpu_data()[0] = loss;
}

template<typename Dtype>
void AffinityLossLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {

  bool huber = this->layer_param_.affinity_loss_param().pseudohuber();
  Dtype scale = this->layer_param_.affinity_loss_param().scale();
  Dtype delta = this->layer_param_.affinity_loss_param().delta();
  Dtype maxgrad = this->layer_param_.affinity_loss_param().maxgrad();
  Dtype ranklossm = this->layer_param_.affinity_loss_param().ranklossmult();
  const Dtype *labels = bottom[1]->cpu_data();

  for (int i = 0; i < 2; ++i) {
    if (propagate_down[i]) {
      const Dtype sign = (i == 0) ? 1 : -1;
      if(huber) {
        //x/(1+(x/delta)^2)
        const Dtype *diff = diff_.cpu_data();
        Dtype *out = bottom[i]->mutable_cpu_diff();
        Dtype mult = sign * scale * top[0]->cpu_diff()[0] / bottom[i]->num();
        for(unsigned j = 0, n = bottom[i]->count(); j < n; j++) {
          Dtype x = diff[j];
          Dtype val = x/delta;
          out[j] = mult * x/(1.0 + val*val);
        }
      } else {
        const Dtype sign = (i == 0) ? 1 : -1;
        const Dtype alpha = sign * scale * top[0]->cpu_diff()[0] / bottom[i]->num();
        caffe_cpu_axpby(bottom[i]->count(),              // count
            alpha,                              // alpha
            diff_.cpu_data(),                   // a
            Dtype(0),                           // beta
            bottom[i]->mutable_cpu_diff());  // b
      }
    }
  }


  if(ranklossm > 0) {
    bool includeneg = this->layer_param_.affinity_loss_param().ranklossneg();
    //compute a rank loss
    int num = bottom[0]->num();
    for (unsigned i = 0; i < num; i++) {
      Dtype labeli = labels[i];
      if(labeli > 0)
      {
        for (unsigned j = i + 1; j < num; j++) {
            Dtype labelj = labels[j];
            //correctly rank good poses
            if(labelj > 0)
            {
                compute_pair_gradient(top, bottom, i, j, ranklossm);
            }
        }
      } else if(includeneg) { //negative, positive should be before
        if(i > 0 && labels[i-1] > 0 && labels[i-1] == -labels[i]) {
          compute_pair_gradient(top, bottom, i, i-1, ranklossm);
        }
      }
    }
  }


  if(maxgrad > 0) {
      //cap the maximum value of the gradient
      for(unsigned i = 0, n = bottom[0]->num(); i < n; i++) {
          Dtype x = bottom[0]->cpu_diff()[i];
          Dtype sign = x < 0 ? -1 : 1;
          x = fabs(x);
          bottom[0]->mutable_cpu_diff()[i] = sign*maxgrad*x/(x+maxgrad);
      }
  }
 /* 
   for(unsigned i = 0, n = bottom[0]->num(); i < n; i++) {
   LOG(INFO) << "AFFGRAD " << i << " " << bottom[0]->cpu_diff()[i];
   }
   */
   
}

template<typename Dtype>
Dtype AffinityLossLayer<Dtype>::compute_pair_loss(
        const vector<Blob<Dtype>*>& bottom, unsigned i, unsigned j, Dtype mult) {
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

    if (Li > Lj) { //if one is negative, then should be less regardless
        Pij = std::max(Pij, Dtype(kLOG_THRESHOLD));
        loss = -log(Pij);
    } else {
        Dtype val = 1-Pij;
        val = std::max(val, Dtype(kLOG_THRESHOLD));
        loss = -log(val);
    }
//LOG(INFO) << "RANKPAIRLOSS " << i <<","<<j<<" "<<loss<< " " << Pij << " " << ediff << " s: " << si << "," << sj << " L: " << Li << "," << Lj;
    return mult*loss;
}
    
template <typename Dtype>
void AffinityLossLayer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom)
{
    //Backward_cpu(top, propagate_down, bottom);
    bottom[0]->mutable_cpu_diff()[0] = bottom[0]->cpu_data()[0];
}

template<typename Dtype>
void AffinityLossLayer<Dtype>::compute_pair_gradient(
        const vector<Blob<Dtype>*>& top, const vector<Blob<Dtype>*>& bottom,
        unsigned i, unsigned j, Dtype mult) {
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
    } else {
        CHECK(false) << "RankLoss requires single score per example or binary class probabilities.";
    }

    Dtype Li = bottom_label[i];
    Dtype Lj = bottom_label[j];
    Dtype diff = si-sj;
    Dtype ediff = exp(diff);
    Dtype d = 0;
    if(std::isfinite(ediff))
        d = - 1.0 / (1.0+exp(diff));

    if(Li < Lj) { //we only care about order, so okay if one is negative
        //reverse sign
        d = -d;
    }

    //scale by number of pairs, which is computed in forward
    CHECK_GT(nranklosspairs, 0) << "Invalid nranklosspairs";
    d /= nranklosspairs; //and also by batch size
    //also by total loss
    Dtype scale = top[0]->cpu_diff()[0];
    d *= scale*mult; //and mult

    if(dim == 1) {
        bottom_diff[i] += d;
        bottom_diff[j] += -d;
        //LOG(INFO) << "RANKLOSS1 " << i << " L:" << Li << " s:" << si << " "<< d << "\t" << j << " L:" << Lj << " s:" << sj << " " << -d << "\n";

    }
}
INSTANTIATE_CLASS(AffinityLossLayer);
REGISTER_LAYER_CLASS(AffinityLoss);

}  // namespace caffe
