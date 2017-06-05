#include <vector>

#include "caffe/layers/conv_layer.hpp"

namespace caffe {

template <typename Dtype>
void ConvolutionLayer<Dtype>::compute_output_shape() {
  const int* kernel_shape_data = this->kernel_shape_.cpu_data();
  const int* stride_data = this->stride_.cpu_data();
  const int* pad_data = this->pad_.cpu_data();
  const int* dilation_data = this->dilation_.cpu_data();
  this->output_shape_.clear();
  for (int i = 0; i < this->num_spatial_axes_; ++i) {
    // i + 1 to skip channel axis
    const int input_dim = this->input_shape(i + 1);
    const int kernel_extent = dilation_data[i] * (kernel_shape_data[i] - 1) + 1;
    const int output_dim = (input_dim + 2 * pad_data[i] - kernel_extent)
        / stride_data[i] + 1;
    this->output_shape_.push_back(output_dim);
  }
}

template <typename Dtype>
void ConvolutionLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  const Dtype* weight = this->blobs_[0]->cpu_data();
  for (int i = 0; i < bottom.size(); ++i) {
    const Dtype* bottom_data = bottom[i]->cpu_data();
    Dtype* top_data = top[i]->mutable_cpu_data();
    for (int n = 0; n < this->num_; ++n) {
      this->forward_cpu_gemm(bottom_data + n * this->bottom_dim_, weight,
          top_data + n * this->top_dim_);
      if (this->bias_term_) {
        const Dtype* bias = this->blobs_[1]->cpu_data();
        this->forward_cpu_bias(top_data + n * this->top_dim_, bias);
      }
    }
  }
}

template <typename Dtype>
void ConvolutionLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
      const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
  const Dtype* weight = this->blobs_[0]->cpu_data();
  Dtype* weight_diff = this->blobs_[0]->mutable_cpu_diff();
  for (int i = 0; i < top.size(); ++i) {
    const Dtype* top_diff = top[i]->cpu_diff();
    const Dtype* bottom_data = bottom[i]->cpu_data();
    Dtype* bottom_diff = bottom[i]->mutable_cpu_diff();
    // Bias gradient, if necessary.
    if (this->bias_term_ && this->param_propagate_down_[1]) {
      Dtype* bias_diff = this->blobs_[1]->mutable_cpu_diff();
      for (int n = 0; n < this->num_; ++n) {
        this->backward_cpu_bias(bias_diff, top_diff + n * this->top_dim_);
      }
    }
    if (this->param_propagate_down_[0] || propagate_down[i]) {
      for (int n = 0; n < this->num_; ++n) {
        // gradient w.r.t. weight. Note that we will accumulate diffs.
        if (this->param_propagate_down_[0]) {
          this->weight_cpu_gemm(bottom_data + n * this->bottom_dim_,
              top_diff + n * this->top_dim_, weight_diff);
        }
        // gradient w.r.t. bottom data, if necessary.
        if (propagate_down[i]) {
          this->backward_cpu_gemm(top_diff + n * this->top_dim_, weight,
              bottom_diff + n * this->bottom_dim_);
        }
      }
    }
  }
}

template <typename Dtype>
void ConvolutionLayer<Dtype>::Backward_relevance(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom, const float eps) {

    //recalculate z_ij, otherwise relu is applied
    Forward_cpu(bottom, top);

    float top_sum = 0;

    for(int i = 0; i < top[0]->count(); ++i)
    {
        top_sum += top[0]->cpu_diff()[i];
    }
    std::cout << "CONV TOP SUM: " << top_sum << '\n';

    const Dtype* weight = this->blobs_[0]->cpu_data();

    for (int i = 0; i < top.size(); ++i)
    {
        const Dtype* top_diff = top[i]->cpu_diff();
        const Dtype* bottom_data = bottom[i]->cpu_data();
        Dtype* bottom_diff = bottom[i]->mutable_cpu_diff();
        caffe_set(bottom[i]->count(), Dtype(0.0), bottom_diff);

        const Dtype* top_data = top[i]->cpu_data();

        //int Mfull = this->num_output_;

        const int first_spatial_axis = this->channel_axis_ + 1;
        int N = bottom[i]->count(first_spatial_axis);
        //int K = this->blobs_[0]->count(1);

        Blob<Dtype> top_data_with_eps((top[i])->shape());

        int outcount = top_data_with_eps.count();

        Dtype* relevance = top_data_with_eps.mutable_cpu_data();
        caffe_copy<Dtype>(outcount, top_diff, relevance);

        Blob<Dtype> bias_removed (top[i]->shape());
        Dtype* bias_removed_data = bias_removed.mutable_cpu_data();

        //stores data
        Blob<Dtype> x_ij (bottom[i]->shape());
        Dtype* x_ij_data = x_ij.mutable_cpu_data();

        //will store data * weights
        Blob<Dtype> z_ij (top[i]->shape());
        Dtype* z_ij_data = z_ij.mutable_cpu_data();

        //copy bottom_data into x_ij_data
        caffe_copy<Dtype>(x_ij.count(), bottom_data, x_ij_data);

        //std::cout << "n: " << this->num_ << '\n';
        for (int n = 0; n < this->num_; ++n) 
        {
            this->forward_cpu_gemm(x_ij_data + n * this->bottom_dim_, weight,
                z_ij_data + n * this->top_dim_);
        }

        for(int c = 0; c < outcount; ++c)
        {
            Dtype bias = this->blobs_[1]->cpu_data()[c/N];
            Dtype val = top_data[c] - bias;
            if(val > 0)
            {
              relevance[c] /= val + eps;
            }
            else if(val < 0)
            {
              relevance[c] /= val - eps;
            }

            //calculate w_ij from top for sanity check
            bias_removed_data[c] = val;
        }

        long double z_ij_total = 0;
        long double bias_removed_total = 0;
        long double top_data_total = 0;

        //std::cout << "shapes: " << top[i]->count() << '|' << z_ij.count() << '\n';

        for(int c = 0; c < outcount; c++)
        {
            //std::cout << z_ij_data[c] << "|" << bias_removed_data[c] << '|' << relevance[c] << '\n';
            //std::cout << bias_removed_data[c] - z_ij_data[c] << '\n';
            z_ij_total += z_ij_data[c];
            bias_removed_total = bias_removed_data[c];
            top_data_total += top_data[c];
        }
        std::cout << z_ij_total << "|" << bias_removed_total <<  '|' << top_data_total << '\n';

        for (int n = 0; n < this->num_; ++n)
        {
            this->backward_cpu_gemm(relevance + n * this->top_dim_, weight, bottom_diff + n * this->bottom_dim_);

            for(int d = 0; d < this->bottom_dim_; ++d)
            {
                bottom_diff[d + n * this->bottom_dim_] 
                        *= bottom_data[d + n * this->bottom_dim_];
            }
        }
    }

    float bottom_sum = 0;
    for (int i = 0; i < bottom[0]->count(); ++i)
    {
        bottom_sum += bottom[0]->cpu_diff()[i];
    }
    std::cout << "CONV BOTTOM SUM: " << bottom_sum << '\n';

}

#ifdef CPU_ONLY
STUB_GPU(ConvolutionLayer);
#endif

INSTANTIATE_CLASS(ConvolutionLayer);

}  // namespace caffe
