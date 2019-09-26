#pragma once
#include "loss.h"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include "cnn_scorer.h"
#include "tee.h"
#include "molgetter.h"
#include <boost/algorithm/string.hpp>

using namespace caffe;
typedef MolGridDataLayer<float> mgridT;

void setLigand(model& m, std::vector<float3>& coords, std::vector<smt>& smtypes);

void setReceptor(model& m, std::vector<float3>& coords, std::vector<smt>& smtypes);

void do_exact_vs(LayerParameter param, caffe::Net<float>& net, 
    std::string vsfile, std::vector<std::string>& ref_ligs,
    std::vector<caffe::shared_ptr<std::ostream> >& out, 
    bool gpu, std::string dist_method, float positive_threshold, float negative_threshold, 
    bool compute_cost=true);

void do_constant_fill(float* fillgrid, size_t gsize, float fillval);
