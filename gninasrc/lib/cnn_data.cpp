/*
 * cnn_data.cpp
 *
 *  Created on: Aug 13, 2018
 *      Author: dkoes
 *
 *  Initializes default models.  These are specified in def files in the models folder.
 *  Model weights are converted to C code with xxd -i
 */

#include "cnn_data.h"
#include <boost/assign.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <vector>

std::string default_model_name = "crossdock_default2018";

boost::unordered_map<std::string, cnn_model_def> cnn_models = boost::assign::map_list_of(
    "default2017", cnn_model_def(
    #include "models/default2017.def"
    ))(
    "crossdock_default2018", cnn_model_def(
    #include "models/default2018.def"
#include "models/weights/crossdock_2018_0.def"
    ))(
    "crossdock_default2018_1", cnn_model_def(
    #include "models/default2018.def"
#include "models/weights/crossdock_2018_1.def"
    ))(
    "crossdock_default2018_2", cnn_model_def(
    #include "models/default2018.def"
#include "models/weights/crossdock_2018_2.def"
    ))(
    "crossdock_default2018_3", cnn_model_def(
    #include "models/default2018.def"
#include "models/weights/crossdock_2018_3.def"
    ))(
    "crossdock_default2018_4", cnn_model_def(
    #include "models/default2018.def"
#include "models/weights/crossdock_2018_4.def"
    ))(
        "general_default2018", cnn_model_def(
        #include "models/default2018.def"
    #include "models/weights/general_2018_0.def"
        ))(
        "general_default2018_1", cnn_model_def(
        #include "models/default2018.def"
    #include "models/weights/general_2018_1.def"
        ))(
        "general_default2018_2", cnn_model_def(
        #include "models/default2018.def"
    #include "models/weights/general_2018_2.def"
        ))(
        "general_default2018_3", cnn_model_def(
        #include "models/default2018.def"
    #include "models/weights/general_2018_3.def"
        ))(
        "general_default2018_4", cnn_model_def(
        #include "models/default2018.def"
    #include "models/weights/general_2018_4.def"
        ))(
            "dense", cnn_model_def(
            #include "models/dense.def"
        #include "models/weights/crossdock_dense_0.def"
            ))(
            "dense_1", cnn_model_def(
            #include "models/dense.def"
        #include "models/weights/crossdock_dense_1.def"
            ))(
            "dense_2", cnn_model_def(
            #include "models/dense.def"
        #include "models/weights/crossdock_dense_2.def"
            ))(
            "dense_3", cnn_model_def(
            #include "models/dense.def"
        #include "models/weights/crossdock_dense_3.def"
            ))(
            "dense_4", cnn_model_def(
            #include "models/dense.def"
        #include "models/weights/crossdock_dense_4.def"
            ));

std::string builtin_cnn_models()
{
  std::vector<std::string> names;
  names.reserve(cnn_models.size());
  for (auto kv : cnn_models) {
    names.push_back(kv.first);
  }
  sort(names.begin(), names.end());
  return boost::algorithm::join(names, " ");
}


