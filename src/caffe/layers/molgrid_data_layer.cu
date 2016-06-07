#include <stdint.h>
#include <vector>


#include "caffe/layers/molgrid_data_layer.hpp"


//gridding is implemented in gridmaker
#include "gnina/src/lib/gridmaker.cu"

namespace caffe {

template <typename Dtype>
void MolGridDataLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
		const vector<Blob<Dtype>*>& top)
{
	forward(bottom, top, true);
}

INSTANTIATE_LAYER_GPU_FORWARD(MolGridDataLayer);

}  // namespace caffe
