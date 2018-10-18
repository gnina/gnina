#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"

namespace caffe {

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  const FlexLSTMParameter& param = this->layer_param_.flex_lstm_param();
  const MolGridDataParameter& mgrid_param = this->layer_param_.molgrid_data_param();
  cube_stride = param.stride();
  if (cube_stride) {
    pattern = strided_cube;
  }
  num_timesteps = bottom[0]->shape(0);
  batch_size = bottom[0]->shape(1);
  ntypes = bottom[0]->shape(2);
  dim = bottom[0]->shape(3);
  unsigned resolution = mgrid_param.resolution();
  unsigned subgrid_dim_in_angstroms = mgrid_param.subgrid_dim();
  subgrid_dim = ::round(subgrid_dim_in_angstroms / resolution) + 1;
  example_size = ntypes * dim * dim * dim;
  current_timestep = 0;
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  //bottom is data, current_x; top is current_x
  const int num_instances = bottom[0]->shape(1);
  const int num_channels = bottom[0]->shape(2);
  //current_x 1xBxCxSdimxSdimxSdim
  CHECK_EQ(6, bottom[1]->num_axes());
  CHECK_EQ(1, bottom[1]->shape(0));
  CHECK_EQ(num_instances, bottom[1]->shape(1));
  CHECK_EQ(num_channels, bottom[1]->shape(2));
  CHECK_EQ(subgrid_dim, bottom[1]->shape(3));
  CHECK_EQ(subgrid_dim, bottom[1]->shape(4));
  CHECK_EQ(subgrid_dim, bottom[1]->shape(5));
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  //set up current_x = GetData<apat>
  if (pattern == AccessPattern::strided_cube) {
    strided_cube_data_handler handler;
    handler.GetData(bottom[0]->cpu_data(), top[0]->mutable_cpu_data(), batch_size, ntypes, 
        subgrid_dim, dim, current_timestep, cube_stride, example_size);
  }

  if (current_timestep != num_timesteps - 1)
    ++current_timestep;
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
  //update the data blob contents to be correct for the previous timestep
  if (pattern == AccessPattern::strided_cube) {
    strided_cube_data_handler handler;
    handler.GetData(bottom[0]->cpu_data(), top[0]->mutable_cpu_data(), batch_size, ntypes, 
        subgrid_dim, dim, current_timestep, cube_stride, example_size);
    //also accumulate gradients for the current timestep in the right location
    handler.AccumulateDiff(top[0]->cpu_diff(), bottom[0]->mutable_cpu_diff(), batch_size, 
        ntypes, subgrid_dim, dim, current_timestep, cube_stride, example_size);
  }
  if (current_timestep != 0)
    --current_timestep;
}

#ifdef CPU_ONLY
// STUB_GPU(LSTMDataGetterLayer);
#endif

INSTANTIATE_CLASS(LSTMDataGetterLayer);
REGISTER_LAYER_CLASS(LSTMDataGetter);
} // namespace caffe
