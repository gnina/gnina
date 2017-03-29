#include <stdint.h>
#include <vector>


#include "caffe/layers/molgrid_data_layer.hpp"


//gridding is implemented in gridmaker
#include "gninasrc/lib/gridmaker.cu"

namespace caffe {

template <typename Dtype>
void MolGridDataLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
		const vector<Blob<Dtype>*>& top)
{
	forward(bottom, top, true);
}

template <typename Dtype>
void MolGridDataLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& bottom,
		const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& top)
{
	backward(bottom, top, true);
}

INSTANTIATE_LAYER_GPU_FORWARD(MolGridDataLayer);
INSTANTIATE_LAYER_GPU_BACKWARD(MolGridDataLayer);

}  // namespace caffe
