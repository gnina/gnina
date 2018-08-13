/*
 * cnn_data.h
 *
 *  Created on: Dec 21, 2016
 *      Author: dkoes
 *
 *  Model(s) and weight(s) for CNNs compiled into the code (used xxd -i )
 *  TODO: support multiple models generically indexed by name
 *  TODO: generate this file from the raw data automatically.
 *
 *
 *  Create c str from binary:
 *  xxd -i train_affinity_g0_p1_rec1_astrat0_b1_h4_iter_150000.caffemodel > x
 *
 */

#ifndef SRC_LIB_CNN_DATA_H_
#define SRC_LIB_CNN_DATA_H_

#include <boost/unordered_map.hpp>
#include <string>

struct cnn_model_def {
    const char *model;
    const char **recmap; //these are null terminated arrays of strings
    const char **ligmap;
    const unsigned char *weights;
    const int num_weights;

    cnn_model_def(const char *m, const char **r, const char **l, const unsigned char *w, int n):
      model(m), recmap(r), ligmap(l), weights(w), num_weights(n) {}
    cnn_model_def(): model(NULL), recmap(NULL), ligmap(NULL), weights(NULL), num_weights(0) {}
};


extern boost::unordered_map<std::string, cnn_model_def> cnn_models;
extern std::string builtin_cnn_models(); //return available names

#endif
