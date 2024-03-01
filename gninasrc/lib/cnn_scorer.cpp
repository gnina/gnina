/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"
#include "gridoptions.h"

#include "caffe/layers/pooling_layer.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/math_functions.hpp"
#include <boost/algorithm/string.hpp>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include "cnn_data.h"

using namespace caffe;
using namespace std;
using namespace boost::algorithm;

static void setup_mgridparm(caffe::MolGridDataParameter *mgridparam, const cnn_options &opts, const string &mname) {
  // initialize grid parameters
  if (mgridparam == NULL) {
    throw usage_error("First layer of model must be MolGridData.");
  }
  mgridparam->set_inmemory(true);
  mgridparam->set_subgrid_dim(opts.subgrid_dim);
  mgridparam->set_batch_size(1);

  if (mname.size() != 0) { // using built-in model, load reg maps
    const char *recmap = cnn_models[mname].recmap;
    const char *ligmap = cnn_models[mname].ligmap;
    mgridparam->set_mem_recmap(recmap);
    mgridparam->set_mem_ligmap(ligmap);
  }

  if (opts.cnn_rotations > 0) {
    mgridparam->set_random_rotation(true);
  } else {
    mgridparam->set_random_rotation(false);
  }
  mgridparam->set_random_translate(0);
}

// initialize from commandline options
// throw error if missing required info
CNNScorer::CNNScorer(const cnn_options &opts) : DLScorer(opts) {
  if (cnnopts.cnn_scoring == CNNnone)
    return; // no cnn

  if (cnnopts.cnn_model_names.size() == 0 && cnnopts.cnn_models.size() == 0) {
    // not specified, use a default ensemble
    // this has been selected to provide the best docking performance
    // with a minimal amount of models
    cnnopts.cnn_model_names.push_back("dense");
    cnnopts.cnn_model_names.push_back("general_default2018_3");
    cnnopts.cnn_model_names.push_back("dense_3");
    cnnopts.cnn_model_names.push_back("crossdock_default2018");
    cnnopts.cnn_model_names.push_back("redock_default2018_2");
  }

  if (cnnopts.cnn_models.size() != cnnopts.cnn_weights.size()) {
    throw usage_error("Different numbers of models and weights specified");
  }

  vector<string> model_names;
  // expand ensembles
  for (const auto &name : cnnopts.cnn_model_names) {
    if (ends_with(name, "_ensemble")) {
      // get everything that starts with the prefix before _ensemble
      string prefix = replace_last_copy(name, "_ensemble", "");
      for (const auto &item : cnn_models) {
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
    NetParameter param;
    if (cnn_models.count(name) == 0) {
      throw usage_error("Invalid model name: " + name);
    }

    const char *model = cnn_models[name].model;
    google::protobuf::io::ArrayInputStream modeldata(model, strlen(model));
    bool success = google::protobuf::TextFormat::Parse(&modeldata, &param);
    if (!success)
      throw usage_error("Error with built-in cnn model " + name);
    UpgradeNetAsNeeded("default", &param);

    param.mutable_state()->set_phase(TEST);

    LayerParameter *first = param.mutable_layer(0);
    mgridparams.push_back(first->mutable_molgrid_data_param());
    setup_mgridparm(mgridparams.back(), cnnopts, name);

    param.set_force_backward(true);
    auto net = caffe::shared_ptr<caffe::Net<Dtype>>(new Net<Dtype>(param));
    nets.push_back(net);

    // load weights
    NetParameter wparam;

    const unsigned char *weights = cnn_models[name].weights;
    unsigned int nbytes = cnn_models[name].num_bytes;

    google::protobuf::io::ArrayInputStream weightdata(weights, nbytes);
    google::protobuf::io::CodedInputStream strm(&weightdata);
#if GOOGLE_PROTOBUF_VERSION >= 3006000
    strm.SetTotalBytesLimit(INT_MAX);
#else
    strm.SetTotalBytesLimit(INT_MAX, 536870912);
#endif
    success = wparam.ParseFromCodedStream(&strm);
    if (!success)
      throw usage_error("Error with default weights.");

    net->CopyTrainedLayersFrom(wparam);
  }

  // load external models
  for (unsigned i = 0, n = cnnopts.cnn_models.size(); i < n; i++) {
    const string &mfile = cnnopts.cnn_models[i];
    const string &wfile = cnnopts.cnn_weights[i];
    NetParameter param;

    ReadNetParamsFromTextFileOrDie(mfile, &param);
    param.mutable_state()->set_phase(TEST);
    LayerParameter *first = param.mutable_layer(0);
    mgridparams.push_back(first->mutable_molgrid_data_param());
    setup_mgridparm(mgridparams.back(), cnnopts, "");

    param.set_force_backward(true);
    auto net = caffe::shared_ptr<caffe::Net<Dtype>>(new Net<Dtype>(param));
    nets.push_back(net);
    net->CopyTrainedLayersFrom(wfile);
  }

  // check that networks matches our expectations and set mgrids
  for (unsigned i = 0, n = nets.size(); i < n; i++) {
    auto net = nets[i];

    // the first layer must be MolGridLayer
    const vector<caffe::shared_ptr<Layer<Dtype>>> &layers = net->layers();
    auto mgrid = dynamic_cast<MolGridDataLayer<Dtype> *>(layers[0].get());
    if (mgrid == NULL) {
      throw usage_error("First layer of model must be MolGridDataLayer.");
    }
    mgrids.push_back(mgrid);

    // we also need an output layer
    if (layers.size() < 1) {
      throw usage_error("No layers in model!");
    }

    if (!net->has_blob("output")) {
      throw usage_error("Model must have output layer named \"output\".");
    }
    if (!net->has_blob("loss")) {
      throw usage_error(
          "Model must have loss calculation layer named \"loss\" (to compute gradient for optimization).");
    }
    if (net->blob_by_name("output")->count() != 2) {
      throw usage_error("Model output layer does not have exactly two outputs.");
    }
  }
}

// returns gradient scores per atom
// assumes necessary pass (backward or backward_relevance) has already been done
std::unordered_map<string, float> CNNScorer::get_gradient_norm_per_atom(bool receptor) {
  std::unordered_map<string, gfloat3> gradient;
  if (mgrids.size() != 1)
    throw usage_error("Gradient visualization does not support model ensembles yet");
  auto &mgrid = mgrids[0];
  if (receptor) {
    mgrid->getReceptorAtoms(0, atoms);
    mgrid->getMappedReceptorGradient(0, gradient);
  } else {
    mgrid->getLigandAtoms(0, atoms);
    mgrid->getMappedLigandGradient(0, gradient);
  }

  std::unordered_map<string, float> scores;

  for (std::pair<string, gfloat3> pair : gradient) {
    // sqrt(x^2 + y^2 + z^2)
    float x = pair.second.x;
    float y = pair.second.y;
    float z = pair.second.z;
    scores[pair.first] = sqrt(x * x + y * y + z * z);
  }

  return scores;
}

// assumes backwards_relevance has been called
std::unordered_map<string, float> CNNScorer::get_relevance_per_atom(bool receptor) {
  std::unordered_map<string, float> relevance;
  if (mgrids.size() != 1)
    throw usage_error("Relevance visualization does not support model ensembles yet");
  auto &mgrid = mgrids[0];
  if (receptor) {
    mgrid->getReceptorAtoms(0, atoms);
    mgrid->getMappedReceptorRelevance(0, relevance);
  } else {
    mgrid->getLigandAtoms(0, atoms);
    mgrid->getMappedLigandRelevance(0, relevance);
  }

  return relevance;
}

void CNNScorer::lrp(const model &m, const string &layer_to_ignore, bool zero_values) {
  boost::lock_guard<boost::recursive_mutex> guard(*mtx);
  if (mgrids.size() != 1)
    throw usage_error("Relevance visualization does not support model ensembles yet");
  auto &mgrid = mgrids[0];
  auto net = nets[0];
  caffe::Caffe::set_random_seed(cnnopts.seed); // same random rotations for each ligand..

  setLigand(m);
  setReceptor(m);

  mgrid->setReceptor(receptor_coords, receptor_smtypes);
  mgrid->setLigand(ligand_coords, ligand_smtypes);
  mgrid->setLabels(1); // for now pose optimization only
  mgrid->enableLigandGradients();
  mgrid->enableReceptorGradients();

  net->Forward();
  if (zero_values) {
    outputDX("zero_blob", 1.0, true, layer_to_ignore, zero_values);
  } else {
    if (layer_to_ignore == "") {
      net->Backward_relevance();
    } else {
      net->Backward_relevance(layer_to_ignore);
    }
  }
}

// do forward and backward pass for gradient visualization
void CNNScorer::gradient_setup(const model &m, const string &recname, const string &ligname,
                               const string &layer_to_ignore) {
  boost::lock_guard<boost::recursive_mutex> guard(*mtx);
  if (mgrids.size() != 1)
    throw usage_error("Gradient visualization does not support model ensembles yet");
  auto &mgrid = mgrids[0];
  auto net = nets[0];
  caffe::Caffe::set_random_seed(cnnopts.seed); // same random rotations for each ligand..

  setLigand(m);
  setReceptor(m);

  mgrid->setReceptor(receptor_coords, receptor_smtypes);
  mgrid->setLigand(ligand_coords, ligand_smtypes);
  mgrid->setLabels(1); // for now pose optimization only
  mgrid->enableLigandGradients();
  mgrid->enableReceptorGradients();

  net->Forward();

  if (layer_to_ignore.length() == 0) {
    net->Backward();
  } else // have to skip layer
  {
    net->Backward_ignore_layer(layer_to_ignore);
  }

  if (ligname.size() > 0) {
    mgrid->getLigandGradient(0, gradient);
    mgrid->getLigandAtoms(0, atoms);
    mgrid->getLigandChannels(0, channels);
    outputXYZ(ligname, atoms, channels, gradient);
  }

  if (recname.size() > 0) {
    mgrid->getReceptorGradient(0, gradient);
    mgrid->getReceptorAtoms(0, atoms);
    mgrid->getReceptorChannels(0, channels);
    outputXYZ(recname, atoms, channels, gradient);
  }
}
// has an affinity prediction layer
bool CNNScorer::has_affinity() const {
  if (nets.size() < 1)
    return false;
  return (bool)nets[0]->blob_by_name("predaff");
}

// populate score and aff with current network output
void CNNScorer::get_net_output(caffe::shared_ptr<caffe::Net<Dtype>> &net, Dtype &score, Dtype &aff, Dtype &loss) {
  const caffe::shared_ptr<Blob<Dtype>> outblob = net->blob_by_name("output");
  const caffe::shared_ptr<Blob<Dtype>> lossblob = net->blob_by_name("loss");
  const caffe::shared_ptr<Blob<Dtype>> affblob = net->blob_by_name("predaff");

  const Dtype *out = outblob->cpu_data();
  score = out[1];
  aff = 0.0;
  if (affblob) {
    aff = affblob->cpu_data()[0];
  }

  loss = lossblob->cpu_data()[0];
}

// Get ligand (and flexible receptor) gradient
void CNNScorer::getGradient(caffe::MolGridDataLayer<Dtype> *mgrid) {
  gradient.resize(num_atoms);
  gradient_lig.reserve(num_atoms);

  // Get ligand gradient
  mgrid->getLigandGradient(0, gradient_lig);

  // Get receptor gradient
  std::vector<gfloat3> gradient_rec;
  if (receptor_map.size() != 0) { // Optimization of flexible residues
    mgrid->getReceptorGradient(0, gradient_rec);
  }

  for (sz i = 0, n = ligand_map.size(); i < n; i++) {
    gradient[ligand_map[i]] = gradient_lig[i];
  }
  for (sz i = 0, n = receptor_map.size(); i < n; i++) {
    gradient[receptor_map[i]] = gradient_rec[i];
  }

  CHECK_EQ(gradient.size(), ligand_map.size() + receptor_map.size());
}

// return score of model, assumes receptor has not changed from initialization
// also sets affinity (if available) and loss (for use with minimization)
// if compute_gradient is set, also adds cnn atom gradient to m.minus_forces
// if maintain center, it will not reposition the molecule
// ALERT: clears minus forces
float CNNScorer::score(model &m, bool compute_gradient, float &affinity, float &loss, float &variance) {
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

  unsigned nscores = nets.size() * max(cnnopts.cnn_rotations, 1U);
  vector<float> affinities;
  if (nscores > 1)
    affinities.reserve(nscores);
  for (unsigned i = 0, n = nets.size(); i < n; i++) {
    caffe::Caffe::set_random_seed(cnnopts.seed); // same random rotations for each ligand..
    auto net = nets[i];
    auto mgrid = mgrids[i];

    if (!isnan(cnnopts.cnn_center[0])) {
      mgrid->setGridCenter(cnnopts.cnn_center);
      current_center = mgrid->getGridCenter();
    } else if (!isnan(current_center[0])) {
      mgrid->setGridCenter(current_center);
    }

    mgrid->setLigand(ligand_coords, ligand_smtypes, true);

    // don't move receptor
    mgrid->setReceptor(receptor_coords, receptor_smtypes);
    current_center = mgrid->getGridCenter(); // has been recalculated from ligand
    if (cnnopts.verbose) {
      std::cout << "current center: ";
      current_center.print(std::cout);
      std::cout << "\n";
    }

    if (compute_gradient || cnnopts.outputxyz) {
      mgrid->enableLigandGradients();
      if (cnnopts.outputxyz) {
        mgrid->enableReceptorGradients();
      } else if (receptor_map.size() != 0) {
        mgrid->enableReceptorGradients(); // rmeli: TODO flexres gradients only
      }
    }

    mgrid->setLabels(1); // for now pose optimization only
    for (unsigned r = 0, n = max(cnnopts.cnn_rotations, 1U); r < n; r++) {
      Dtype s = 0, a = 0, l = 0;
      net->Forward(); // do all rotations at once if requested

      get_net_output(net, s, a, l);
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
        net->Backward();
        // Get gradient from mgrid into CNNScorer::gradient
        getGradient(mgrid);

        // Update ligand (and flexible residues) gradient
        m.add_minus_forces(gradient);
      }
      cnt++;
    } // end rotations

    if (cnnopts.outputxyz) {
      const string &ligname = cnnopts.xyzprefix + "_lig.xyz";
      const string &recname = cnnopts.xyzprefix + "_rec.xyz";

      mgrid->getLigandGradient(0, gradient);
      mgrid->getLigandAtoms(0, atoms);
      mgrid->getLigandChannels(0, channels);
      outputXYZ(ligname, atoms, channels, gradient);

      mgrid->getReceptorGradient(0, gradient); // rmeli: TODO Full gradient or just flexres?
      mgrid->getReceptorAtoms(0, atoms);
      mgrid->getReceptorChannels(0, channels);
      outputXYZ(recname, atoms, channels, gradient);
    }

    if (cnnopts.gradient_check) {
      check_gradient(net);
    }

    if (cnnopts.outputdx && i == 0) {
      // DANGER! This modifies the values in the network
      outputDX(m.get_name());
    }
  } // end models loop

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
float CNNScorer::score(model &m, float &variance) {
  float aff = 0;
  float loss = 0;
  return score(m, false, aff, loss, variance);
}

// To aid in debugging, will compute the gradient at the
// grid level, apply it with different multiples, and evaluate
// the effect. Perhaps may evaluate atom gradients as well?
//
// IMPORTANT: assumes model is already setup
// Prints out the results
void CNNScorer::check_gradient(caffe::shared_ptr<caffe::Net<Dtype>> &net) {
  Dtype origscore = 0;
  Dtype origaff = 0;
  Dtype origloss = 0;
  Dtype newscore = 0.0;
  Dtype newaff = 0.0;
  Dtype newloss = 0.0;

  std::cout << std::scientific;
  Dtype lambda = 1.0;
  for (unsigned i = 0; i < 4; i++) {
    // score pose
    net->Forward();
    get_net_output(net, origscore, origaff, origloss);

    // backprop
    net->Backward();

    // get grid and diff blobs
    const caffe::shared_ptr<Blob<Dtype>> datablob = net->blob_by_name("data");
    Dtype *data = datablob->mutable_cpu_data();
    Dtype *diff = datablob->mutable_cpu_diff();

    // apply gradient
    caffe_cpu_axpby(datablob->count(), -lambda, diff, 1.0f, data); // sets data

    // propagate forward, starting _after_ molgrid
    net->ForwardFrom(1);

    // compare scores
    get_net_output(net, newscore, newaff, newloss);

    std::cout << "CHECK   " << origscore - newscore << "    LAMBDA: " << lambda << "   OLD: " << origscore << ","
              << origaff << "   NEW: " << newscore << "," << newaff << std::endl;

    // test a single channel
    unsigned channel = 16 + 15; // compile time constant !!!
    net->Forward();
    net->Backward();

    vector<int> inds;
    inds.push_back(0);
    inds.push_back(channel);
    unsigned off = datablob->offset(inds);
    data = datablob->mutable_cpu_data();
    diff = datablob->mutable_cpu_diff();

    unsigned n = datablob->count(2);
    for (unsigned i = 0; i < n; i++) {
      if (*(data + off + i) > 0) // only modify if density actually exists
        *(data + off + i) += -lambda * *(diff + off + i);
    }
    // caffe_cpu_axpby(datablob->count(2), -lambda, diff+off, 1.0f, data+off);
    net->ForwardFrom(1);
    get_net_output(net, newscore, newaff, newloss);

    std::cout << "CHECKch " << origscore - newscore << " " << channel << " LAMBDA: " << lambda
              << "   OLD: " << origscore << "," << origaff << "   NEW: " << newscore << "," << newaff << std::endl;

    // super expensive - evaluate every grid point for channel
    net->Forward();
    get_net_output(net, origscore, origaff, origloss);
    for (unsigned i = 0; i < n; i++) {
      data = datablob->mutable_cpu_data();
      diff = datablob->mutable_cpu_diff();
      float gval = *(data + off + i);
      float gdiff = *(diff + off + i);
      if (gval > 0) {
        *(data + off + i) += lambda;
        net->ForwardFrom(1);
        get_net_output(net, newscore, newaff, newloss);
        std::cout << "GRIDch " << i << " gval " << gval << " gdiff " << gdiff << " lambda " << lambda << " change "
                  << origscore - newscore << "\n";

        data = datablob->mutable_cpu_data();
        diff = datablob->mutable_cpu_diff();
        *(data + off + i) = max(gval - lambda, 0.0f);
        net->ForwardFrom(1);
        get_net_output(net, newscore, newaff, newloss);
        std::cout << "GRIDch " << i << " gval " << gval << " gdiff " << gdiff << " lambda " << -lambda << " change "
                  << origscore - newscore << "\n";
        *(data + off + i) = gval;
      }
    }

    lambda /= 10.0;
  }
}

// dump dx files of the diff, but only for first network
// zero_values: run backward relevance with only dead node values
void CNNScorer::outputDX(const string &prefix, double scale, bool lrp, string layer_to_ignore, bool zero_values)

{
  auto net = nets[0];
  auto mgrid = mgrids[0];

  const caffe::shared_ptr<Blob<Dtype>> datablob = net->blob_by_name("data");

  const vector<caffe::shared_ptr<Layer<Dtype>>> &layers = net->layers();
  //    if(datablob) {
  // this is a big more fragile than I would like.. if there is a pooling layer before
  // the first convoluational of fully connected layer and it is a max pooling layer,
  // change it to average before the backward to avoid a discontinuous map
  PoolingLayer<Dtype> *pool = NULL;
  for (unsigned i = 1, nl = layers.size(); i < nl; i++) {
    pool = dynamic_cast<PoolingLayer<Dtype> *>(layers[i].get());
    if (pool)
      break; // found it
    else if (layers[i]->type() == string("Convolution"))
      break; // give up
    else if (layers[i]->type() == string("InnerProduct"))
      break;
  }

  if (pool) {
    if (pool->pool() == PoolingParameter_PoolMethod_MAX) {
      pool->set_pool(PoolingParameter_PoolMethod_AVE);
    } else {
      pool = NULL; // no need to reset to max
    }
  }

  // must redo backwards with average pooling
  if (lrp) {
    net->Backward_relevance(layer_to_ignore, zero_values);
  } else
    net->Backward();

  string p = prefix;
  if (p.length() == 0)
    p = "dx";
  mgrid->dumpDiffDX(p, datablob.get(), scale);

  if (pool) {
    pool->set_pool(PoolingParameter_PoolMethod_MAX);
  }
}

void CNNScorer::outputXYZ(const string &base, const vector<gfloat3> &atoms, const vector<short> &whichGrid,
                          const vector<gfloat3> &gradient) {
  const char *sym[] = {"C", "C", "C", "C",  "Ca", "Fe", "Mg", "N", "N", "N", "N", "O", "O", "P", "S", "Zn", "C",
                       "C", "C", "C", "Br", "Cl", "F",  "N",  "N", "N", "N", "O", "O", "O", "P", "S", "S",  "I"};

  ofstream out(base.c_str());
  out.precision(5);

  // Count number of valid (>= 0) whichGrid entries
  size_t n_atoms{0};
  for (size_t i = 0, n = atoms.size(); i < n; ++i) {
    if (whichGrid[i] >= 0) {
      n_atoms++;
    }
  }
  out << n_atoms << "\n\n"; // xyz

  // Print coordinates and gradients
  for (unsigned i = 0, n = atoms.size(); i < n; ++i) {

    // Skip invalid channel
    if (whichGrid[i] < 0) {
      continue;
    }

    out << sym[whichGrid[i]] << " ";
    out << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z << " ";
    out << gradient[i].x << " " << gradient[i].y << " " << gradient[i].z;
    if (i + 1 < n)
      out << "\n";
  }
}

void CNNScorer::set_bounding_box(grid_dims &box) const {
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