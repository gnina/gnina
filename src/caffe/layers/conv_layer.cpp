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
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {

        const float eps = .00001;

        /*
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

            int Mfull = this->num_output_;

            const int first_spatial_axis = this->channel_axis_ + 1;
            int N = bottom[i]->count(first_spatial_axis);
            int K = this->blobs_[0]->count(1);

            Blob<Dtype> top_data_with_eps((top[i])->shape());

            int outcount = top_data_with_eps.count();

            Dtype* top_data_with_eps_data = top_data_with_eps.mutable_cpu_data();
            caffe_copy<Dtype>(outcount, top_diff, top_data_with_eps_data);

            for(int c = 0; c < outcount; ++c)
            {
                if(top_data[c] > 0)
                {
                    top_data_with_eps_data[c] /= top_data[c] + eps;
                }
                else if(top_data[c] < 0)
                {
                    top_data_with_eps_data[c] /= top_data[c] - eps;
                }
            }

            for (int n = 0; n < this->num_; ++n)
            {
                this->backward_cpu_gemm(top_data_with_eps_data + n * this->top_dim_, weight, bottom_diff + n * this->bottom_dim_);

                for(int d = 0; d < this->bottom_dim_; ++d)
                {
                    bottom_diff[d + n * this->bottom_dim_] 
                            *= bottom_data[d + n * this->bottom_dim_];
                }
            }
            }
        */
        float top_sum = 0;
        for (int i = 0; i < top[0]->count(); ++i)
        {
            top_sum += top[0]->cpu_diff()[i];
            //std::cout << top[0]->cpu_diff()[i] << "|";
            //std::cout << bottom[0]->cpu_diff()[i] << "|";
        }
        std::cout << "CONV TOP SUM: " << top_sum << '\n';


       // std::cout <<  top[0]->count() << "\n";
        //std::cout <<  bottom[0]->count() << "\n";

        int i = 0; //assume only using top[0]
        const Dtype* weight = this->blobs_[i]->cpu_data();
        const Dtype* top_diff = top[i]->cpu_diff();
        const Dtype* bottom_data = bottom[i]->cpu_data();
        Dtype* bottom_diff = bottom[i]->mutable_cpu_diff();
        caffe_set(bottom[i]->count(), Dtype(0.), bottom_diff);
        for (int n = 0; n < this->num_; ++n)
        {
            Dtype* alphabetas = this->alphabeta(top_diff+n *this->top_dim_,
                            weight, bottom_data + n * this->bottom_dim_,
                            bottom_diff + n * this->bottom_dim_);

            //std::cout << "AB CONV DATA: " << '\n';
            //for (int i = 0; i < 1000; ++i)
            //{
            //    std::cout << alphabetas[i] << "|";
            //}
            caffe_copy(bottom[i]->count(), alphabetas, bottom_diff);
        }



        float bottom_sum = 0;
        for (int i = 0; i < bottom[0]->count(); ++i)
        {
            bottom_sum += bottom[0]->cpu_diff()[i];
            //std::cout << bottom[0]->cpu_diff()[i] << "|";
        }
        std::cout << "CONV BOTTOM SUM: " << bottom_sum << '\n';

}

#ifdef CPU_ONLY
STUB_GPU(ConvolutionLayer);
#endif

INSTANTIATE_CLASS(ConvolutionLayer);

}  // namespace caffe
