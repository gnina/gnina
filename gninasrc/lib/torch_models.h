/*
 * torch_models.h
 *
 *  Created on: Feb 26, 2024
 *      Author: dkoes
 *
 *  Model(s) and weight(s) for torch script models. 
 *  Binary model data is embedded using the linker:
 *  ld -r -b binary -o binary.o foo.bar  # then link in binary.o
 *  (see https://stackoverflow.com/questions/4158900/embedding-resources-in-executable-using-gcc)
 *  torch_models.cpp will be dynamically generated using make_model_cpp.py
 */

#ifndef SRC_LIB_TORCH_MODELS_H_
#define SRC_LIB_TORCH_MODELS_H_

#include <boost/unordered_map.hpp>
#include <string>

extern boost::unordered_map<std::string, std::pair<char*, char*> > torch_models;
extern std::string builtin_torch_models(); //return available names

#endif
