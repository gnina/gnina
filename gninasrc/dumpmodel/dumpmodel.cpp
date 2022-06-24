/*
 * dumpmodel.cpp
 *
 * Output a built-in caffemodel.
 *
 *  Created on: Jun 24, 2022
 *      Author: dkoes
 */

#include <iostream>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>
#include "cnn_data.h"


int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cerr << "Need to specify default model name: " + builtin_cnn_models() << "\n";
    return -1;
  }
  std::string name = argv[1];
  std::string outfile = name+".caffemodel";
  if(boost::filesystem::exists(outfile)) {
    std::cerr << outfile << " already exists. Giving up.\n";
    return -1;
  }
  if(cnn_models.count(name) == 0) {
    std::cerr << name << " is not a valid model name:\n"+builtin_cnn_models() << "\n";
    return -1;
  }
  const cnn_model_def& model = cnn_models[name];

  std::ofstream fout;
  fout.open(outfile.c_str(), std::ios::binary | std::ios::out);
  fout.write((const char*)model.weights, model.num_bytes);
  fout.close();
}



