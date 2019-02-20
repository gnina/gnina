#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"

namespace caffe {

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  //Infer the data access pattern from the parameters and update relevant
  //data members. 
  const LSTMDataGetterParameter& param = this->layer_param_.lstm_datagetter_param();
  cube_stride = param.stride();
  if (cube_stride) {
    pattern = AccessPatterns::strided_cube;
    Dtype subgrid_dim_in_angstroms = param.subgrid_dim();
    Dtype resolution = param.resolution();
    subgrid_dim = ::round(subgrid_dim_in_angstroms / resolution) + 1;
    dim = bottom[0]->shape(2);
    unsigned slices_per_dim = ((dim - subgrid_dim) / cube_stride) + 1;
    num_timesteps = slices_per_dim * slices_per_dim * slices_per_dim;
    example_size = ntypes * dim * dim * dim;
  }
  else {
    std::cerr << "Flex LSTM layer currently only supports a strided cube access pattern";
    exit(-1);
  }
  batch_size = bottom[0]->shape(0);
  ntypes = bottom[0]->shape(1);
  current_timestep = this->layer_param_.lstm_datagetter_param().timestep();
  Reshape(bottom, top);
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  //bottom is data, top is current_x
  vector<int> current_x_shape;
  //if access_pattern == strided_cube, current_x 1xBxCxSdimxSdimxSdim
  switch(pattern) {
    case AccessPatterns::strided_cube:
      {
        //for the strided cube pattern, the top and bottom dims are the same
        //for each timestep; the top is a shared buffer whose contents change
        //but the shape stays the same. so the first datagetter sets the shape
        //and the subsequent layers just check that the shape doesn't deviate
        //from what they expect
        if (current_timestep == 0) {
          current_x_shape.push_back(1);
          current_x_shape.push_back(batch_size);
          current_x_shape.push_back(ntypes);
          current_x_shape.push_back(subgrid_dim);
          current_x_shape.push_back(subgrid_dim);
          current_x_shape.push_back(subgrid_dim);
          top[0]->Reshape(current_x_shape);
        }
        else {
          CHECK_EQ(6, top[0]->num_axes());
          CHECK_EQ(1, top[0]->shape(0));
          CHECK_EQ(batch_size, top[0]->shape(1));
          CHECK_EQ(ntypes, top[0]->shape(2));
          CHECK_EQ(subgrid_dim, top[0]->shape(3));
          CHECK_EQ(subgrid_dim, top[0]->shape(4));
          CHECK_EQ(subgrid_dim, top[0]->shape(5));
        }
        break;
      }
    default:
      {
        CHECK_LT(pattern, AccessPatterns::num_patterns) << "Invalid access pattern " << pattern;
      }
  }
}

template <typename Dtype>
void GetData(const Dtype* src, Dtype* dest, unsigned dim, unsigned subgrid_dim, 
    unsigned x_offset, unsigned y_offset, unsigned z_offset, unsigned batch_size, 
    unsigned ntypes) {
  //strided cube version:
  //extract a single subcube corresponding to the correct stride,
  //starting at our properly offset (x,y,z) indices
  unsigned overall_size = dim * dim * dim;
  unsigned example_size = ntypes * dim * dim * dim;
  for (unsigned batch_idx=0; batch_idx < batch_size; ++batch_idx) {
    for (unsigned grid=0; grid < ntypes; ++grid) {
      for (unsigned i=0; i<subgrid_dim; ++i) {
        for (unsigned j=0; j<subgrid_dim; ++j) {
          for (unsigned k=0; k<subgrid_dim; ++k) {
            dest[(((batch_idx * ntypes + grid) * subgrid_dim + i) * subgrid_dim + j) * 
              subgrid_dim + k] = src[batch_idx * example_size + grid * overall_size + 
              x_offset * dim * dim + y_offset * dim + z_offset + 
              ((i * dim) + j) * dim + k];
          }
        }
      }
    }
  }
}

template <typename Dtype>
void AccumulateDiff(const Dtype* src, Dtype* dest, unsigned dim, unsigned subgrid_dim, 
    unsigned x_offset, unsigned y_offset, unsigned z_offset, unsigned batch_size, 
    unsigned ntypes) {
  //strided cube version:
  //accumulate diff for subcube into the correct location in the full grid
  unsigned overall_size = dim * dim * dim;
  unsigned example_size = ntypes * dim * dim * dim;
  for (unsigned batch_idx=0; batch_idx < batch_size; ++batch_idx) {
    for (unsigned grid=0; grid < ntypes; ++grid) {
      for (unsigned i=0; i<subgrid_dim; ++i) {
        for (unsigned j=0; j<subgrid_dim; ++j) {
          for (unsigned k=0; k<subgrid_dim; ++k) {
            dest[batch_idx * example_size + grid * overall_size + 
               x_offset * dim * dim + y_offset * dim + z_offset + 
               ((i * dim) + j) * dim + k] += 
            src[(((batch_idx * ntypes + grid) * subgrid_dim + i) * subgrid_dim + j) * 
              subgrid_dim + k];
          }
        }
      }
    }
  }
}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
  //Update blob that will be used as input to the RNN at this timestep
  //according to the chosen access pattern
  const Dtype* src = bottom[0]->cpu_data();
  Dtype* dest = top[0]->mutable_cpu_data();
  switch(pattern) {
    case AccessPatterns::strided_cube:
      {
        //use the current_timestep to find the location of the first value in
        //the subcube we're going to use at this timestep; this is our starting
        //offset
        unsigned factor = (((dim - subgrid_dim) / cube_stride) + 1);
        unsigned x_offset = ((current_timestep / (factor * factor)) % factor) * cube_stride;
        unsigned y_offset = ((current_timestep / factor) % factor) * cube_stride;
        unsigned z_offset = (current_timestep % factor) * cube_stride;
        GetData(src, dest, dim, subgrid_dim, x_offset, y_offset, z_offset, batch_size, ntypes);
        break;
      }
    default:
      {
        CHECK_LT(pattern, AccessPatterns::num_patterns) << "Invalid access pattern " << pattern;
      }
  }

}

template <typename Dtype>
void LSTMDataGetterLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
  //- use diff computed for the per-timestep blob to compute the relevant part
  //of the diff we're building up for the full input
  //
  //- also update the data blob contents to be correct for the *previous* timestep - 
  //by the time the DataGetter layer is hit during backward, the blobs that need
  //current_x to be set to the correct contents for its timestep have already
  //computed their diffs with it, so now we set up the contents to work for the
  //layers before it
  Dtype* total_diff = bottom[0]->mutable_cpu_diff();
  //if we're just starting backward, zero grid diff...is this necessary?
  if (current_timestep == num_timesteps-1) {
    memset(total_diff, 0, num_timesteps * batch_size * ntypes * dim * dim * dim);
  }
  switch(pattern) {
    case AccessPatterns::strided_cube:
      {
        unsigned factor = (((dim - subgrid_dim) / cube_stride) + 1);
        unsigned x_offset = ((current_timestep / (factor * factor)) % factor) * cube_stride;
        unsigned y_offset = ((current_timestep/ factor) % factor) * cube_stride;
        unsigned z_offset = (current_timestep % factor) * cube_stride;
        //accumulate gradients for the current timestep in the right location
        AccumulateDiff(top[0]->cpu_diff(), total_diff, dim, subgrid_dim, 
            x_offset, y_offset, z_offset, batch_size, ntypes);
        if (current_timestep > 0) {
          x_offset -= cube_stride;
          y_offset -= cube_stride;
          z_offset -= cube_stride;
          GetData(bottom[0]->cpu_data(), top[0]->mutable_cpu_data(), dim, subgrid_dim, 
              x_offset, y_offset, z_offset, batch_size, ntypes);
        }
        break;
      }
  default:
      {
        CHECK_LT(pattern, AccessPatterns::num_patterns) << "Invalid access pattern " << pattern;
      }
  }
}

#ifdef CPU_ONLY
// STUB_GPU(LSTMDataGetterLayer);
#endif

INSTANTIATE_CLASS(LSTMDataGetterLayer);
REGISTER_LAYER_CLASS(LSTMDataGetter);
} // namespace caffe
