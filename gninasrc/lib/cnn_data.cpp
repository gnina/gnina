/*
 * cnn_data.cpp
 *
 *  Created on: Aug 13, 2018
 *      Author: dkoes
 *
 *  Initializes default models.  These are specified in def files in the models folder.
 */

#include "cnn_data.h"
#include <boost/assign.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <vector>


boost::unordered_map<std::string, cnn_model_def> cnn_models = boost::assign::map_list_of(
    "default2017", cnn_model_def(
        #include "models/default2017.def"
        )
        );



std::string builtin_cnn_models()
{
  std::vector<std::string> names; names.reserve(cnn_models.size());
  for(auto kv: cnn_models) {
    names.push_back(kv.first);
  }
  sort(names.begin(), names.end());
  return boost::algorithm::join(names," ");
}
