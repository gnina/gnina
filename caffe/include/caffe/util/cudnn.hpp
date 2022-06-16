#ifndef CAFFE_UTIL_CUDNN_H_
#define CAFFE_UTIL_CUDNN_H_
#ifdef USE_CUDNN

#include <cudnn.h>
#include <vector>
#include <gflags/gflags.h>

#include "caffe/common.hpp"
#include "caffe/proto/caffe.pb.h"

DECLARE_bool(use_tensor_core);

#define CUDNN_VERSION_MIN(major, minor, patch) \
    (CUDNN_VERSION >= (major * 1000 + minor * 100 + patch))

#define CUDNN_CHECK(condition) \
do { \
  cudnnStatus_t status = condition; \
  if(status != CUDNN_STATUS_SUCCESS) \
    throw std::runtime_error(std::string("CUDNN Error (out of memory?): ")+cudnnGetErrorString(status)); \
} while (0)


namespace caffe {

namespace cudnn {

template <typename Dtype> class dataType;
template<> class dataType<float>  {
 public:
  static const cudnnDataType_t type = CUDNN_DATA_FLOAT;
  static float oneval, zeroval;
  static const void *one, *zero;
};
template<> class dataType<double> {
 public:
  static const cudnnDataType_t type = CUDNN_DATA_DOUBLE;
  static double oneval, zeroval;
  static const void *one, *zero;
};

template <typename Dtype>
inline void createTensorDesc(cudnnTensorDescriptor_t* desc) {
  CUDNN_CHECK(cudnnCreateTensorDescriptor(desc));
}

template <typename Dtype>
inline void setTensor4dDesc(cudnnTensorDescriptor_t* desc,
    int n, int c, int h, int w,
    int stride_n, int stride_c, int stride_h, int stride_w) {
  CUDNN_CHECK(cudnnSetTensor4dDescriptorEx(*desc, dataType<Dtype>::type,
        n, c, h, w, stride_n, stride_c, stride_h, stride_w));
}

template <typename Dtype>
inline void setTensorNdDesc(cudnnTensorDescriptor_t* desc, const std::vector<int>& shape_) {
  std::vector<int> shape(shape_);
  while(shape.size() < 4) {
    //apparently tensors can't be too small (not documented)
    shape.push_back(1);
  }
  CUDNN_CHECK(cudnnSetTensorNdDescriptorEx(*desc, CUDNN_TENSOR_NCHW, dataType<Dtype>::type, shape.size(), &shape[0]));
}

template <typename Dtype>
inline void setTensor4dDesc(cudnnTensorDescriptor_t* desc,
    int n, int c, int h, int w) {
  const int stride_w = 1;
  const int stride_h = w * stride_w;
  const int stride_c = h * stride_h;
  const int stride_n = c * stride_c;
  setTensor4dDesc<Dtype>(desc, n, c, h, w,
                         stride_n, stride_c, stride_h, stride_w);
}

template <typename Dtype>
inline void createFilterDesc(cudnnFilterDescriptor_t* desc, const std::vector<int>& shape) {
  CUDNN_CHECK(cudnnCreateFilterDescriptor(desc));
  CUDNN_CHECK(cudnnSetFilterNdDescriptor(*desc, dataType<Dtype>::type,
      CUDNN_TENSOR_NCHW, shape.size(), &shape[0]));
}

template <typename Dtype>
inline void createConvolutionDesc(cudnnConvolutionDescriptor_t* conv) {
  CUDNN_CHECK(cudnnCreateConvolutionDescriptor(conv));
  if(FLAGS_use_tensor_core) {
    CUDNN_CHECK(cudnnSetConvolutionMathType(*conv, CUDNN_TENSOR_OP_MATH));
  }
}

template<typename Dtype>
inline void setConvolutionDesc(cudnnConvolutionDescriptor_t* conv,
    cudnnTensorDescriptor_t bottom, cudnnFilterDescriptor_t filter,
    const std::vector<int>& pad, const std::vector<int>& stride) {
  CHECK_EQ(pad.size(),stride.size())<< "Mismatched pad and stride in setConvolutionDesc";
  unsigned n = pad.size();
  std::vector<int> dilation(n, 1);

  CUDNN_CHECK(cudnnSetConvolutionNdDescriptor(*conv, n, &pad[0], &stride[0], &dilation[0],
          CUDNN_CROSS_CORRELATION,  dataType<Dtype>::type));
}


template <typename Dtype>
inline void createActivationDescriptor(cudnnActivationDescriptor_t* activ_desc,
    cudnnActivationMode_t mode) {
  CUDNN_CHECK(cudnnCreateActivationDescriptor(activ_desc));
  CUDNN_CHECK(cudnnSetActivationDescriptor(*activ_desc, mode,
                                           CUDNN_PROPAGATE_NAN, Dtype(0)));
}

}  // namespace cudnn

}  // namespace caffe

#endif  // USE_CUDNN
#endif  // CAFFE_UTIL_CUDNN_H_
