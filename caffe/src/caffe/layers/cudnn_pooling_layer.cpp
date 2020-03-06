#ifdef USE_CUDNN
#include <vector>

#include "caffe/layers/cudnn_pooling_layer.hpp"

namespace caffe {


template<typename Dtype>
void CuDNNPoolingLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  PoolingLayer<Dtype>::LayerSetUp(bottom, top);
  CUDNN_CHECK(cudnnCreate(&handle_));

  //check pooling type
  switch (this->layer_param_.pooling_param().pool()) {
  case PoolingParameter_PoolMethod_MAX:
    mode_ = CUDNN_POOLING_MAX_DETERMINISTIC;
    break;
  case PoolingParameter_PoolMethod_AVE:
    mode_ = CUDNN_POOLING_AVERAGE_COUNT_EXCLUDE_PADDING;
    break;
  default:
    LOG(FATAL)<< "Pooling method not supported by CUDNN.";
  }

  //allocate memory
  cudnn::createTensorDesc<Dtype>(&bottom_desc_);
  cudnn::createTensorDesc<Dtype>(&top_desc_);
  CUDNN_CHECK(cudnnCreatePoolingDescriptor(&pooling_desc_));

  cudnnSetPoolingNdDescriptor(pooling_desc_, mode_, CUDNN_PROPAGATE_NAN,
      this->num_spatial_axes_, this->kernel_shape_.cpu_data(),
      this->pad_.cpu_data(), this->stride_.cpu_data());

  handles_setup_ = true;
}

template <typename Dtype>
void CuDNNPoolingLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  PoolingLayer<Dtype>::Reshape(bottom, top);

  if(this->global_pooling_) {
    //kernel is size of input blob
    int* kernel_shape_data = this->kernel_shape_.mutable_cpu_data();
    for (int i = 0; i < this->num_spatial_axes_; ++i) {
      kernel_shape_data[i] = bottom[0]->shape(this->channel_axis_ + i+1);
    }
    cudnnSetPoolingNdDescriptor(pooling_desc_, mode_, CUDNN_PROPAGATE_NAN,
        this->num_spatial_axes_, this->kernel_shape_.cpu_data(),
        this->pad_.cpu_data(), this->stride_.cpu_data());
  }
  cudnn::setTensorNdDesc<Dtype>(&bottom_desc_, bottom[0]->shape());
  cudnn::setTensorNdDesc<Dtype>(&top_desc_,top[0]->shape());
}

template <typename Dtype>
CuDNNPoolingLayer<Dtype>::~CuDNNPoolingLayer() {
  // Check that handles have been setup before destroying.
  if (!handles_setup_) { return; }

  cudnnDestroyTensorDescriptor(bottom_desc_);
  cudnnDestroyTensorDescriptor(top_desc_);
  cudnnDestroyPoolingDescriptor(pooling_desc_);
  cudnnDestroy(handle_);
}

INSTANTIATE_CLASS(CuDNNPoolingLayer);

}   // namespace caffe
#endif
