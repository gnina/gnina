#include <string>
#include <vector>

#include "caffe/blob.hpp"
#include "caffe/common.hpp"
#include "caffe/filler.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template<typename Dtype>
template<AccessPattern apat>
void Flex_LSTMLayer::Forward_cpu {
  // unclear if there is a payoff to not just doing data duplication here
}

template<typename Dtype>
template<AccessPattern apat>
void Flex_LSTMLayer::Forward_gpu {
  // rather than explicitly unrolling the net and doing separate kernel
  // launches for each dummy layer, we'll handle everything in a single kernel
  // launch and share data via warp shuffles and shared memory whenever
  // possible
}
}  // namespace caffe
