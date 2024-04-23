/*
 * cnn_torch_scorer.cpp
 *
 *  Created on: Dec 21, 2023
 *      Author: dkoes
 */

#include "cnn_torch_scorer.h"
#include "common.h"
#include "gridoptions.h"
#include "torch_models.h"

#include "torch_model.h"
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <torch/script.h>

using namespace std;
using namespace boost::algorithm;

// initialize from commandline options
// throw error if missing required info
template <bool isCUDA> CNNTorchScorer<isCUDA>::CNNTorchScorer(const cnn_options &opts, tee *log) : DLScorer(opts) {
  if (cnnopts.cnn_scoring == CNNnone)
    return; // no cnn

  if(cnnopts.cnn_models.size() == 0) {
    if (cnnopts.cnn_model_names.size() == 0) {
      // not specified, use a default ensemble
      // this has been selected to provide the best docking performance
      // with a minimal amount of models
      cnnopts.cnn_model_names.push_back("dense_1_3");
      cnnopts.cnn_model_names.push_back("dense_1_3_PT_KD_3");
      cnnopts.cnn_model_names.push_back("crossdock_default2018_KD_4");
    } else if (cnnopts.cnn_model_names.size() == 1 && cnnopts.cnn_model_names[0] == "fast") {
      cnnopts.cnn_model_names[0] = "all_default_to_default_1_3_1";
    }
  }

  if (cnnopts.cnn_weights.size() > 0) {
    throw usage_error("Torch models do not have separate weight files.");
  }

  vector<string> model_names;
  // expand ensembles
  for (const auto &name : cnnopts.cnn_model_names) {
    if (ends_with(name, "_ensemble")) {
      // get everything that starts with the prefix before _ensemble
      string prefix = replace_last_copy(name, "_ensemble", "");
      for (const auto &item : torch_models) {
        if (starts_with(item.first, prefix)) {
          model_names.push_back(item.first);
        }
      }
    } else {
      model_names.push_back(name);
    }
  }

  // load built-in models
  for (const auto &name : model_names) {
    using namespace boost::iostreams;

    if (torch_models.count(name) == 0) {
      throw usage_error("Invalid model name: " + name);
    }

    const char *modelstart = torch_models[name].first;
    const char *modelend = torch_models[name].second;

    basic_array_source<char> model_source(modelstart, modelend - modelstart);
    stream<basic_array_source<char>> model_stream(model_source);
    models.push_back(std::make_shared<TorchModel<isCUDA>>(model_stream, name, log));
  }

  // load external models
  for (unsigned i = 0, n = cnnopts.cnn_models.size(); i < n; i++) {
    string fname(cnnopts.cnn_models[i]);
    ifstream in(fname.c_str());
    if (!in)
      throw usage_error("Could not open file " + fname);
    models.push_back(std::make_shared<TorchModel<isCUDA>>(in, fname, log));
  }

  // check that networks matches our expectations - how?
}

// has an affinity prediction layer
template <bool isCUDA> bool CNNTorchScorer<isCUDA>::has_affinity() const {
  return true; //??
}

// return score of model, assumes receptor has not changed from initialization
// also sets affinity (if available) and loss (for use with minimization)
// if compute_gradient is set, also adds cnn atom gradient to m.minus_forces
// if maintain center, it will not reposition the molecule
// ALERT: clears minus forces
template <bool isCUDA>
float CNNTorchScorer<isCUDA>::score(model &m, bool compute_gradient, float &affinity, float &loss, float &variance) {
  boost::lock_guard<boost::recursive_mutex> guard(*mtx);
  if (!initialized())
    return -1.0;

  // Get ligand atoms and coords from movable atoms
  setLigand(m);
  // Get receptor atoms and flex/inflex coordinats from movable atoms
  setReceptor(m);

  m.clear_minus_forces();
  // these variables will accumulate across models/rotations
  double score = 0.0;
  affinity = 0.0;
  loss = 0.0;
  unsigned cnt = 0;

  // TODO: need number of models
  unsigned nscores = models.size() * max(cnnopts.cnn_rotations, 1U);
  vector<float> affinities;
  if (nscores > 1)
    affinities.reserve(nscores);

  gradient.resize(receptor_coords.size() + ligand_coords.size());
  // loop over models
  for (auto &model : models) {
    torch::manual_seed(cnnopts.seed); // same random rotations for each ligand..
    libmolgrid::random_engine.seed(cnnopts.seed);

    if (!isnan(cnnopts.cnn_center[0])) {
      current_center = cnnopts.cnn_center;
    }

    for (unsigned r = 0, n = max(cnnopts.cnn_rotations, 1U); r < n; r++) {
      Dtype s = 0, a = 0, l = 0;

      vec grid_center = vec(NAN, NAN, NAN); // recalculate from ligand      
      auto output = model->forward(receptor_coords, receptor_smtypes, ligand_coords, ligand_smtypes, grid_center, r > 0,
                                   compute_gradient);
      if (output.size() > 0)
        s = output[0];
      if (output.size() > 1)
        a = output[1];
      if (output.size() > 2)
        l = output[2];

      score += s;
      if (nscores > 1)
        affinities.push_back(a);
      affinity += a;
      loss += l;

      if (cnnopts.cnn_rotations > 1) {
        if (cnnopts.verbose) {
          std::cout << "RotateScore: " << s << "\n";
          if (a)
            std::cout << "RotateAff: " << a << "\n";
        }
      }

      if (compute_gradient || cnnopts.outputxyz) {

        // Get gradient from mgrid into CNNScorer::gradient
        getGradient(model, gradient);

        // Update ligand (and flexible residues) gradient
        m.add_minus_forces(gradient);

      }
      cnt++;
    } // end rotations
  }

  // if there were multiple evaluations, scale appropriately
  if (cnt > 1) {
    m.scale_minus_forces(1.0 / cnt);
  }
  affinity /= cnt;
  loss /= cnt;
  score /= cnt; // mean
  variance = 0;
  if (affinities.size() > 1) {
    float sum = 0;
    for (float s : affinities) {
      float diff = affinity - s;
      diff *= diff;
      sum += diff;
    }
    variance = sum / affinities.size();
  }

  if (cnnopts.verbose)
    std::cout << std::fixed << std::setprecision(10) << "cnnscore " << score << "\n";

  return score;
}

// return only score
template <bool isCUDA> float CNNTorchScorer<isCUDA>::score(model &m, float &variance) {
  float aff = 0;
  float loss = 0;
  return score(m, false, aff, loss, variance);
}

template <bool isCUDA>
void CNNTorchScorer<isCUDA>::getGradient(std::shared_ptr<TorchModel<isCUDA> > model, std::vector<gfloat3> &gradient) {
  gradient.resize(receptor_map.size()+ligand_map.size());
  gradient_lig.reserve(ligand_map.size());

// Get ligand gradient
  model->getLigandGradient(gradient_lig);

// Get receptor gradient
  if (receptor_map.size() != 0)
  { // Optimization of flexible residues
    model->getReceptorGradient(gradient_rec);
  }

  for(sz i = 0, n = ligand_map.size(); i < n; i++) {
    gradient[ligand_map[i]] = gradient_lig[i];
  }
  for(sz i = 0, n = receptor_map.size(); i < n; i++) {
    gradient[receptor_map[i]] = gradient_rec[i];
  }

  VINA_CHECK(gradient.size() == ligand_map.size() + receptor_map.size());
}

template <bool isCUDA> void CNNTorchScorer<isCUDA>::set_bounding_box(grid_dims &box) const {
  vec center = get_center();
  fl dim = get_grid_dim();
  fl n = dim / get_grid_res();
  fl half = dim / 2.0;

  for (unsigned i = 0; i < 3; i++) {
    fl c = center[i];
    box[i].begin = c - half;
    box[i].end = c + half;
    box[i].n = n;
  }
}

template <bool isCUDA> fl CNNTorchScorer<isCUDA>::get_grid_dim() const {
  VINA_CHECK(models.size() > 0);
  return models[0]->get_grid_dim();
}

template <bool isCUDA> fl CNNTorchScorer<isCUDA>::get_grid_res() const {
  VINA_CHECK(models.size() > 0);
  return models[0]->get_grid_res();
}

// explicit instantiations
template class CNNTorchScorer<true>;
template class CNNTorchScorer<false>;