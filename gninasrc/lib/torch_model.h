/*
 * torch_model.h
 *
 *  Created on: Feb 14, 2024
 *      Author: dkoes
 */

#ifndef SRC_LIB_TORCH_MODEL_H_
#define SRC_LIB_TORCH_MODEL_H_

#include <iostream>
#include <jsoncpp/json/json.h>
#include <libmolgrid/atom_typer.h>
#include <libmolgrid/grid_maker.h>
#include <torch/script.h>

#include "tee.h"
#include "atom_type.h"
#include "gpu_math.h"

/** Wrapper class for reading in and executing a TorchScript Model */
template <bool isCUDA> class TorchModel {
  libmolgrid::GridMaker gmaker;
  torch::jit::script::Module module;
  std::shared_ptr<libmolgrid::AtomTyper> lig_typer;
  std::shared_ptr<libmolgrid::AtomTyper> rec_typer;
  std::vector<gfloat3> gradient_rec, gradient_lig;
  bool apply_logistic_loss = false; //more for debugging 
  bool skip_softmax = false;

public:
  TorchModel(std::istream &in, const std::string &name, tee *log);

  std::vector<float> forward(const std::vector<float3> &rec_coords, const std::vector<smt> &rec_types,
               const std::vector<float3> &lig_coords, const std::vector<smt> &lig_types, const vec& center,
               bool rotate, bool compute_gradient);

  //assumes forward was called with compute_gradient
  void getLigandGradient(std::vector<gfloat3>& grad);
  void getReceptorGradient(std::vector<gfloat3>& grad);

  float get_grid_dim() const { return gmaker.get_dimension(); }
  float get_grid_res() const { return gmaker.get_resolution(); }

  std::shared_ptr<libmolgrid::AtomTyper> get_lig_typer() const { return lig_typer; }
  std::shared_ptr<libmolgrid::AtomTyper> get_rec_typer() const { return rec_typer; }
};

#endif /* SRC_LIB_TORCH_MODEL_H_ */