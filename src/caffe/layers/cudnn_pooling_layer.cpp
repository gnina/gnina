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
    mode_ = CUDNN_POOLING_MAX;
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

  cudnnSetTensorNdDescriptorEx(bottom_desc_, CUDNN_TENSOR_NCHW, cudnn::dataType<Dtype>::type, this->input_shape_.count(), this->input_shape_.cpu_data());
  cudnnSetTensorNdDescriptorEx(top_desc_, CUDNN_TENSOR_NCHW, cudnn::dataType<Dtype>::type, this->output_shape_.count(), this->output_shape_.cpu_data());
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
