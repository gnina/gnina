#include <string>
#include <vector>

#include "caffe/blob.hpp"
#include "caffe/common.hpp"
#include "caffe/filler.hpp"
#include "caffe/layer.hpp"
#include "caffe/layers/flex_lstm_layer.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void FlexLSTMLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  // this differs from the base class because the data blob doesn't need to be
  // TxBx..., since the data getter will actually set up the data used at each
  // timestep; we don't assume anything about it and instead use the shape of
  // the labels to determine T_ and N_
  CHECK_EQ(bottom[1]->num_axes(), 2)
      << "bottom[1] must have exactly 2 axes -- (#timesteps, #streams)";
  this->T_ = bottom[1]->shape(0);
  this->N_ = bottom[1]->shape(1);
  LOG(INFO) << "Initializing recurrent layer: assuming input batch contains "
            << this->T_ << " timesteps of " << this->N_ << " independent streams.";

  // If expose_hidden is set, we take as input and produce as output
  // the hidden state blobs at the first and last timesteps.
  this->expose_hidden_ = this->layer_param_.recurrent_param().expose_hidden();

  // Get (recurrent) input/output names.
  vector<string> output_names;
  this->OutputBlobNames(&output_names);
  vector<string> recur_input_names;
  this->RecurrentInputBlobNames(&recur_input_names);
  vector<string> recur_output_names;
  this->RecurrentOutputBlobNames(&recur_output_names);
  const int num_recur_blobs = recur_input_names.size();
  CHECK_EQ(num_recur_blobs, recur_output_names.size());

  // If provided, bottom[2] is a static input to the recurrent net.
  const int num_hidden_exposed = this->expose_hidden_ * num_recur_blobs;
  this->static_input_ = (bottom.size() > 2 + num_hidden_exposed);
  if (this->static_input_) {
    CHECK_GE(bottom[2]->num_axes(), 1);
    CHECK_EQ(this->N_, bottom[2]->shape(0));
  }

  // Create a NetParameter; setup the inputs that aren't unique to particular
  // recurrent architectures.
  NetParameter net_param;

  LayerParameter* input_layer_param = net_param.add_layer();
  input_layer_param->set_type("Input");
  InputParameter* input_param = input_layer_param->mutable_input_param();
  input_layer_param->add_top("x");
  BlobShape input_shape;
  for (int i = 0; i < bottom[0]->num_axes(); ++i) {
    input_shape.add_dim(bottom[0]->shape(i));
  }
  input_param->add_shape()->CopyFrom(input_shape);

  input_shape.Clear();
  for (int i = 0; i < bottom[1]->num_axes(); ++i) {
    input_shape.add_dim(bottom[1]->shape(i));
  }
  input_layer_param->add_top("cont");
  input_param->add_shape()->CopyFrom(input_shape);

  if (this->static_input_) {
    input_shape.Clear();
    for (int i = 0; i < bottom[2]->num_axes(); ++i) {
      input_shape.add_dim(bottom[2]->shape(i));
    }
    input_layer_param->add_top("x_static");
    input_param->add_shape()->CopyFrom(input_shape);
  }

  // Call the child's FillUnrolledNet implementation to specify the unrolled
  // recurrent architecture.
  this->FillUnrolledNet(&net_param);

  // Prepend this layer's name to the names of each layer in the unrolled net.
  const string& layer_name = this->layer_param_.name();
  if (layer_name.size()) {
    for (int i = 0; i < net_param.layer_size(); ++i) {
      LayerParameter* layer = net_param.mutable_layer(i);
      layer->set_name(layer_name + "_" + layer->name());
    }
  }

  // Add "pseudo-losses" to all outputs to force backpropagation.
  // (Setting force_backward is too aggressive as we may not need to backprop to
  // all inputs, e.g., the sequence continuation indicators.)
  vector<string> pseudo_losses(output_names.size());
  for (int i = 0; i < output_names.size(); ++i) {
    LayerParameter* layer = net_param.add_layer();
    pseudo_losses[i] = output_names[i] + "_pseudoloss";
    layer->set_name(pseudo_losses[i]);
    layer->set_type("Reduction");
    layer->add_bottom(output_names[i]);
    layer->add_top(pseudo_losses[i]);
    layer->add_loss_weight(1);
  }

  // Create the unrolled net.
  this->unrolled_net_.reset(new Net<Dtype>(net_param));
  this->unrolled_net_->set_debug_info(
      this->layer_param_.recurrent_param().debug_info());

  // Setup pointers to the inputs.
  this->x_input_blob_ = CHECK_NOTNULL(this->unrolled_net_->blob_by_name("x").get());
  this->cont_input_blob_ = CHECK_NOTNULL(this->unrolled_net_->blob_by_name("cont").get());
  if (this->static_input_) {
    this->x_static_input_blob_ =
        CHECK_NOTNULL(this->unrolled_net_->blob_by_name("x_static").get());
  }

  // Setup pointers to paired recurrent inputs/outputs.
  this->recur_input_blobs_.resize(num_recur_blobs);
  this->recur_output_blobs_.resize(num_recur_blobs);
  for (int i = 0; i < recur_input_names.size(); ++i) {
    this->recur_input_blobs_[i] =
        CHECK_NOTNULL(this->unrolled_net_->blob_by_name(recur_input_names[i]).get());
    this->recur_output_blobs_[i] =
        CHECK_NOTNULL(this->unrolled_net_->blob_by_name(recur_output_names[i]).get());
  }

  // Setup pointers to outputs.
  CHECK_EQ(top.size() - num_hidden_exposed, output_names.size())
      << "OutputBlobNames must provide an output blob name for each top.";
  this->output_blobs_.resize(output_names.size());
  for (int i = 0; i < output_names.size(); ++i) {
    this->output_blobs_[i] =
        CHECK_NOTNULL(this->unrolled_net_->blob_by_name(output_names[i]).get());
  }

  // We should have 2 inputs (x and cont), plus a number of recurrent inputs,
  // plus maybe a static input.
  CHECK_EQ(2 + num_recur_blobs + this->static_input_,
           this->unrolled_net_->input_blobs().size());

  // This layer's parameters are any parameters in the layers of the unrolled
  // net. We only want one copy of each parameter, so check that the parameter
  // is "owned" by the layer, rather than shared with another.
  this->blobs_.clear();
  for (int i = 0; i < this->unrolled_net_->params().size(); ++i) {
    if (this->unrolled_net_->param_owners()[i] == -1) {
      LOG(INFO) << "Adding parameter " << i << ": "
                << this->unrolled_net_->param_display_names()[i];
      this->blobs_.push_back(this->unrolled_net_->params()[i]);
    }
  }
  // Check that param_propagate_down is set for all of the parameters in the
  // unrolled net; set param_propagate_down to true in this layer.
  for (int i = 0; i < this->unrolled_net_->layers().size(); ++i) {
    for (int j = 0; j < this->unrolled_net_->layers()[i]->blobs().size(); ++j) {
      CHECK(this->unrolled_net_->layers()[i]->param_propagate_down(j))
          << "param_propagate_down not set for layer " << i << ", param " << j;
    }
  }
  this->param_propagate_down_.clear();
  this->param_propagate_down_.resize(this->blobs_.size(), true);

  // Set the diffs of recurrent outputs to 0 -- we can't backpropagate across
  // batches.
  for (int i = 0; i < this->recur_output_blobs_.size(); ++i) {
    caffe_set(this->recur_output_blobs_[i]->count(), Dtype(0),
              this->recur_output_blobs_[i]->mutable_cpu_diff());
  }

  // Check that the last output_names.size() layers are the pseudo-losses;
  // set last_layer_index so that we don't actually run these layers.
  const vector<string>& layer_names = this->unrolled_net_->layer_names();
  this->last_layer_index_ = layer_names.size() - 1 - pseudo_losses.size();
  for (int i = this->last_layer_index_ + 1, j = 0; i < layer_names.size(); ++i, ++j) {
    CHECK_EQ(layer_names[i], pseudo_losses[j]);
  }
}

template <typename Dtype>
void FlexLSTMLayer<Dtype>::RecurrentInputShapes(vector<BlobShape>* shapes) const {
  const int num_output = this->layer_param_.recurrent_param().num_output();
  const int num_blobs = 2;
  shapes->resize(num_blobs);
  //h and c
  for (int i = 0; i < num_blobs; ++i) {
    (*shapes)[i].Clear();
    (*shapes)[i].add_dim(1);  // a single timestep
    (*shapes)[i].add_dim(this->N_);
    (*shapes)[i].add_dim(num_output);
  }
}

template <typename Dtype>
void FlexLSTMLayer<Dtype>::FillUnrolledNet(NetParameter* net_param) const {
  const int num_output = this->layer_param_.recurrent_param().num_output();
  const int stride = this->layer_param_.flex_lstm_param().stride();
  const float subgrid_dim = this->layer_param_.flex_lstm_param().subgrid_dim();
  const float resolution = this->layer_param_.flex_lstm_param().resolution();
  CHECK_GT(num_output, 0) << "num_output must be positive";
  const FillerParameter& weight_filler =
      this->layer_param_.recurrent_param().weight_filler();
  const FillerParameter& bias_filler =
      this->layer_param_.recurrent_param().bias_filler();

  // Add generic LayerParameters (without bottoms/tops) of layer types we'll
  // use to save redundant code.
  LayerParameter hidden_param;
  hidden_param.set_type("InnerProduct");
  hidden_param.mutable_inner_product_param()->set_num_output(num_output * 4);
  hidden_param.mutable_inner_product_param()->set_bias_term(false);
  hidden_param.mutable_inner_product_param()->set_axis(2);
  hidden_param.mutable_inner_product_param()->
      mutable_weight_filler()->CopyFrom(weight_filler);

  LayerParameter biased_hidden_param(hidden_param);
  biased_hidden_param.mutable_inner_product_param()->set_bias_term(true);
  biased_hidden_param.mutable_inner_product_param()->
      mutable_bias_filler()->CopyFrom(bias_filler);

  LayerParameter sum_param;
  sum_param.set_type("Eltwise");
  sum_param.mutable_eltwise_param()->set_operation(
      EltwiseParameter_EltwiseOp_SUM);

  LayerParameter scale_param;
  scale_param.set_type("Scale");
  scale_param.mutable_scale_param()->set_axis(0);

  LayerParameter slice_param;
  slice_param.set_type("Slice");
  slice_param.mutable_slice_param()->set_axis(0);

  LayerParameter split_param;
  split_param.set_type("Split");

  vector<BlobShape> input_shapes;
  LSTMLayer<Dtype>::RecurrentInputShapes(&input_shapes);
  CHECK_EQ(2, input_shapes.size());

  LayerParameter* input_layer_param = net_param->add_layer();
  input_layer_param->set_type("Input");
  InputParameter* input_param = input_layer_param->mutable_input_param();

  input_layer_param->add_top("c_0");
  input_param->add_shape()->CopyFrom(input_shapes[0]);

  input_layer_param->add_top("h_0");
  input_param->add_shape()->CopyFrom(input_shapes[1]);

  LayerParameter* cont_slice_param = net_param->add_layer();
  cont_slice_param->CopyFrom(slice_param);
  cont_slice_param->set_name("cont_slice");
  cont_slice_param->add_bottom("cont");
  cont_slice_param->mutable_slice_param()->set_axis(0);

  if (this->static_input_) {
    // Add layer to transform x_static to the gate dimension.
    //     W_xc_x_static = W_xc_static * x_static
    LayerParameter* x_static_transform_param = net_param->add_layer();
    x_static_transform_param->CopyFrom(hidden_param);
    x_static_transform_param->mutable_inner_product_param()->set_axis(1);
    x_static_transform_param->set_name("W_xc_x_static");
    x_static_transform_param->add_param()->set_name("W_xc_static");
    x_static_transform_param->add_bottom("x_static");
    x_static_transform_param->add_top("W_xc_x_static_preshape");
    x_static_transform_param->add_propagate_down(true);

    LayerParameter* reshape_param = net_param->add_layer();
    reshape_param->set_type("Reshape");
    BlobShape* new_shape =
         reshape_param->mutable_reshape_param()->mutable_shape();
    new_shape->add_dim(1);  // One timestep.
    // Should infer this->N as the dimension so we can reshape on batch size.
    new_shape->add_dim(-1);
    new_shape->add_dim(
        x_static_transform_param->inner_product_param().num_output());
    reshape_param->set_name("W_xc_x_static_reshape");
    reshape_param->add_bottom("W_xc_x_static_preshape");
    reshape_param->add_top("W_xc_x_static");
  }

  LayerParameter output_concat_layer;
  output_concat_layer.set_name("h_concat");
  output_concat_layer.set_type("Concat");
  output_concat_layer.add_top("h");
  output_concat_layer.mutable_concat_param()->set_axis(0);

  for (int t = 1; t <= this->T_; ++t) {
    string tm1s = format_int(t - 1);
    string ts = format_int(t);

    cont_slice_param->add_top("cont_" + ts);

    //Add layer to choose the data for this timestep
    LayerParameter* datagetter_param = net_param->add_layer();
    datagetter_param->set_type("LSTMDataGetter");
    datagetter_param->set_name("datagetter_" + tm1s);
    datagetter_param->add_bottom("x");
    datagetter_param->add_top("current_x");
    datagetter_param->mutable_lstm_datagetter_param()->set_timestep(t-1);
    datagetter_param->mutable_lstm_datagetter_param()->set_stride(stride);
    datagetter_param->mutable_lstm_datagetter_param()->set_subgrid_dim(subgrid_dim);
    datagetter_param->mutable_lstm_datagetter_param()->set_resolution(resolution);

    // Add layers to flush the hidden state when beginning a new
    // sequence, as indicated by cont_t.
    //     h_conted_{t-1} := cont_t * h_{t-1}
    //
    // Normally, cont_t is binary (i.e., 0 or 1), so:
    //     h_conted_{t-1} := h_{t-1} if cont_t == 1
    //                       0   otherwise
    {
      LayerParameter* cont_h_param = net_param->add_layer();
      cont_h_param->CopyFrom(scale_param);
      cont_h_param->set_name("h_conted_" + tm1s);
      cont_h_param->add_bottom("h_" + tm1s);
      cont_h_param->add_bottom("cont_" + ts);
      cont_h_param->add_top("h_conted_" + tm1s);
    }

    // Add layer to transform current piece of x to the hidden state dimension.
    //     W_xc_x = W_xc * current_x + b_c
    LayerParameter* current_x_transform_param = net_param->add_layer();
    current_x_transform_param->CopyFrom(biased_hidden_param);
    current_x_transform_param->set_name("x_transform_" + ts);
    current_x_transform_param->add_bottom("current_x");
    current_x_transform_param->add_top("W_xc_x");
    current_x_transform_param->add_propagate_down(true);
  
    // Add layer to compute
    //     W_hc_h_{t-1} := W_hc * h_conted_{t-1}
    {
      LayerParameter* w_param = net_param->add_layer();
      w_param->CopyFrom(hidden_param);
      w_param->set_name("transform_" + ts);
      w_param->add_param()->set_name("W_hc");
      w_param->add_bottom("h_conted_" + tm1s);
      w_param->add_top("W_hc_h_" + tm1s);
      w_param->mutable_inner_product_param()->set_axis(2);
    }

    // Add the outputs of the linear transformations to compute the gate input.
    //     gate_input_t := W_hc * h_conted_{t-1} + W_xc * x_t + b_c
    //                   = W_hc_h_{t-1} + W_xc_x_t + b_c
    {
      LayerParameter* input_sum_layer = net_param->add_layer();
      input_sum_layer->CopyFrom(sum_param);
      input_sum_layer->set_name("gate_input_" + ts);
      input_sum_layer->add_bottom("W_hc_h_" + tm1s);
      input_sum_layer->add_bottom("W_xc_x");
      if (this->static_input_) {
        input_sum_layer->add_bottom("W_xc_x_static");
      }
      input_sum_layer->add_top("gate_input_" + ts);
    }

    // Add LSTMUnit layer to compute the cell & hidden vectors c_t and h_t.
    // Inputs: c_{t-1}, gate_input_t = (i_t, f_t, o_t, g_t), cont_t
    // Outputs: c_t, h_t
    //     [ i_t' ]
    //     [ f_t' ] := gate_input_t
    //     [ o_t' ]
    //     [ g_t' ]
    //         i_t := \sigmoid[i_t']
    //         f_t := \sigmoid[f_t']
    //         o_t := \sigmoid[o_t']
    //         g_t := \tanh[g_t']
    //         c_t := cont_t * (f_t .* c_{t-1}) + (i_t .* g_t)
    //         h_t := o_t .* \tanh[c_t]
    {
      LayerParameter* lstm_unit_param = net_param->add_layer();
      lstm_unit_param->set_type("LSTMUnit");
      lstm_unit_param->add_bottom("c_" + tm1s);
      lstm_unit_param->add_bottom("gate_input_" + ts);
      lstm_unit_param->add_bottom("cont_" + ts);
      lstm_unit_param->add_top("c_" + ts);
      lstm_unit_param->add_top("h_" + ts);
      lstm_unit_param->set_name("unit_" + ts);
    }
    output_concat_layer.add_bottom("h_" + ts);
  }  // for (int t = 1; t <= this->T_; ++t)

  {
    LayerParameter* c_T_copy_param = net_param->add_layer();
    c_T_copy_param->CopyFrom(split_param);
    c_T_copy_param->add_bottom("c_" + format_int(this->T_));
    c_T_copy_param->add_top("c_T");
  }
  net_param->add_layer()->CopyFrom(output_concat_layer);
}

template <typename Dtype>
void FlexLSTMLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  CHECK_EQ(bottom[1]->num_axes(), 2)
      << "bottom[1] must have exactly 2 axes -- (#timesteps, #streams)";
  CHECK_EQ(this->T_, bottom[1]->shape(0)) << "input number of timesteps changed";
  this->N_ = bottom[1]->shape(1);
  this->x_input_blob_->ReshapeLike(*bottom[0]);
  vector<int> cont_shape = bottom[1]->shape();
  this->cont_input_blob_->Reshape(cont_shape);
  if (this->static_input_) {
    this->x_static_input_blob_->ReshapeLike(*bottom[2]);
  }
  vector<BlobShape> recur_input_shapes;
  RecurrentInputShapes(&recur_input_shapes);
  CHECK_EQ(recur_input_shapes.size(), this->recur_input_blobs_.size());
  for (int i = 0; i < recur_input_shapes.size(); ++i) {
    this->recur_input_blobs_[i]->Reshape(recur_input_shapes[i]);
  }
  this->unrolled_net_->Reshape();
  this->x_input_blob_->ShareData(*bottom[0]);
  this->x_input_blob_->ShareDiff(*bottom[0]);
  this->cont_input_blob_->ShareData(*bottom[1]);
  if (this->static_input_) {
    this->x_static_input_blob_->ShareData(*bottom[2]);
    this->x_static_input_blob_->ShareDiff(*bottom[2]);
  }
  if (this->expose_hidden_) {
    const int bottom_offset = 2 + this->static_input_;
    for (int i = bottom_offset, j = 0; i < bottom.size(); ++i, ++j) {
      CHECK(this->recur_input_blobs_[j]->shape() == bottom[i]->shape())
          << "shape mismatch - recur_input_blobs_[" << j << "]: "
          << this->recur_input_blobs_[j]->shape_string()
          << " vs. bottom[" << i << "]: " << bottom[i]->shape_string();
      this->recur_input_blobs_[j]->ShareData(*bottom[i]);
    }
  }
  for (int i = 0; i < this->output_blobs_.size(); ++i) {
    top[i]->ReshapeLike(*(this->output_blobs_[i]));
    top[i]->ShareData(*(this->output_blobs_[i]));
    top[i]->ShareDiff(*(this->output_blobs_[i]));
  }
  if (this->expose_hidden_) {
    const int top_offset = this->output_blobs_.size();
    for (int i = top_offset, j = 0; i < top.size(); ++i, ++j) {
      top[i]->ReshapeLike(*(this->recur_output_blobs_[j]));
    }
  }
}


INSTANTIATE_CLASS(FlexLSTMLayer);
REGISTER_LAYER_CLASS(FlexLSTM);

}  // namespace caffe
