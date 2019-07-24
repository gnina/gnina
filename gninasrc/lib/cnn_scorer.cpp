/*
 * cnn_scorer.cpp
 *
 *  Created on: May 2, 2016
 *      Author: dkoes
 */

#include "cnn_scorer.h"
#include "gridoptions.h"

#include "caffe/proto/caffe.pb.h"
#include "caffe/layers/pooling_layer.hpp"
#include "caffe/util/math_functions.hpp"
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

#include "cnn_data.h"

using namespace caffe;
using namespace std;

//initialize from commandline options
//throw error if missing required info
CNNScorer::CNNScorer(const cnn_options& opts)
    : mgrid(NULL), cnnopts(opts), mtx(new boost::recursive_mutex), current_center(NAN,NAN,NAN) {

  if (cnnopts.cnn_scoring || cnnopts.cnn_refinement) {
    NetParameter param;

    //load cnn model
    if (cnnopts.cnn_model.size() == 0) {
      if(cnn_models.count(cnnopts.cnn_model_name) == 0) {
        throw usage_error("Invalid model name: "+cnnopts.cnn_model_name);
      }

      const char *model = cnn_models[cnnopts.cnn_model_name].model;
      google::protobuf::io::ArrayInputStream modeldata(model, strlen(model));
      bool success = google::protobuf::TextFormat::Parse(&modeldata, &param);
      if (!success) throw usage_error("Error with built-in cnn model "+cnnopts.cnn_model_name);
      UpgradeNetAsNeeded("default", &param);
    } else {
      ReadNetParamsFromTextFileOrDie(cnnopts.cnn_model, &param);
    }

    param.mutable_state()->set_phase(TEST);

    LayerParameter *first = param.mutable_layer(0);
    mgridparam = first->mutable_molgrid_data_param();
    if (mgridparam == NULL) {
      throw usage_error("First layer of model must be MolGridData.");
    }
    mgridparam->set_inmemory(true);
    mgridparam->set_subgrid_dim(opts.subgrid_dim);

    if (cnnopts.cnn_model.size() == 0) {
      const char *recmap = cnn_models[cnnopts.cnn_model_name].recmap;
      const char *ligmap = cnn_models[cnnopts.cnn_model_name].ligmap;
      mgridparam->set_mem_recmap(recmap);
      mgridparam->set_mem_ligmap(ligmap);
    }

    //set batch size to 1
    unsigned bsize = 1;
    //unless we have rotations, in which case do them all at once, which turns out isn't actually faster :-(
    if (cnnopts.cnn_rotations > 0) {
      //let user specify rotations
      unsigned nrot = cnnopts.cnn_rotations;
      mgridparam->set_random_rotation(true);
      //I think there's a bug in the axial rotations - they aren't all distinct
      //BUT it turns out this isn't actually faster
      //bsize = nrot;
    } else {
      mgridparam->set_random_rotation(false);
      mgridparam->set_random_translate(0);
    }

    param.set_force_backward(true);

    net.reset(new Net<Dtype>(param));

    //load weights
    if (cnnopts.cnn_weights.size() == 0) {
      NetParameter wparam;

      const unsigned char *weights = cnn_models[cnnopts.cnn_model_name].weights;
      unsigned int nweights = cnn_models[cnnopts.cnn_model_name].num_weights;

      google::protobuf::io::ArrayInputStream weightdata(weights,nweights);
      google::protobuf::io::CodedInputStream strm(&weightdata);
      strm.SetTotalBytesLimit(INT_MAX, 536870912);
      bool success = wparam.ParseFromCodedStream(&strm);
      if (!success) throw usage_error("Error with default weights.");

      net->CopyTrainedLayersFrom(wparam);
    } else {
      net->CopyTrainedLayersFrom(cnnopts.cnn_weights);
    }

    //check that network matches our expectations
    //the first layer must be MolGridLayer
    const vector<caffe::shared_ptr<Layer<Dtype> > >& layers = net->layers();
    mgrid = dynamic_cast<MolGridDataLayer<Dtype>*>(layers[0].get());
    if (mgrid == NULL) {
      throw usage_error("First layer of model must be MolGridDataLayer.");
    }

    //we also need an output layer
    if (layers.size() < 1) {
      throw usage_error("No layers in model!");
    }

    if (!net->has_blob("output")) {
      throw usage_error("Model must have output layer named \"output\".");
    }
    if (net->blob_by_name("output")->count() != 2 * bsize) {
      throw usage_error(
          "Model output layer does not have exactly two outputs.");
    }
  }

}

//returns gradient scores per atom
//assumes necessary pass (backward or backward_relevance) has already been done
std::unordered_map<string, float> CNNScorer::get_gradient_norm_per_atom(bool receptor) {
  std::unordered_map<string, gfloat3> gradient;

  if (receptor) {
    mgrid->getReceptorAtoms(0, atoms);
    mgrid->getMappedReceptorGradient(0, gradient);
  } else {
    mgrid->getLigandAtoms(0, atoms);
    mgrid->getMappedLigandGradient(0, gradient);
  }

  std::unordered_map<string, float> scores;

  for (std::pair<string, gfloat3> pair : gradient) {
      //sqrt(x^2 + y^2 + z^2)
      float x = pair.second.x;
      float y = pair.second.y;
      float z = pair.second.z;
      scores[pair.first] = sqrt(x * x + y * y + z * z);
  }

  return scores;
}

//assumes backwards_relevance has been called
std::unordered_map<string, float> CNNScorer::get_relevance_per_atom(bool receptor) {
  std::unordered_map<string, float> relevance;

  if (receptor) {
    mgrid->getReceptorAtoms(0, atoms);
    mgrid->getMappedReceptorRelevance(0, relevance);
  } else {
    mgrid->getLigandAtoms(0, atoms);
    mgrid->getMappedLigandRelevance(0, relevance);
  }

  return relevance;
}

void CNNScorer::lrp(const model& m, const string& layer_to_ignore,
    bool zero_values) {
  boost::lock_guard<boost::recursive_mutex> guard(*mtx);

  caffe::Caffe::set_random_seed(cnnopts.seed); //same random rotations for each ligand..

  mgrid->setReceptor(m.get_fixed_atoms());
  mgrid->setLigand(m.get_movable_atoms(), m.coordinates()); // TODO: Movable atoms for receptor
  mgrid->setLabels(1); //for now pose optimization only
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

//do forward and backward pass for gradient visualization
void CNNScorer::gradient_setup(const model& m, const string& recname,
    const string& ligname, const string& layer_to_ignore) {
  boost::lock_guard<boost::recursive_mutex> guard(*mtx);

  caffe::Caffe::set_random_seed(cnnopts.seed); //same random rotations for each ligand..

  mgrid->setReceptor(m.get_fixed_atoms());
  mgrid->setLigand(m.get_movable_atoms(), m.coordinates());
  mgrid->setLabels(1); //for now pose optimization only
  mgrid->enableLigandGradients();
  mgrid->enableReceptorGradients();

  net->Forward();

  if (layer_to_ignore.length() == 0) {
    net->Backward();
  } else //have to skip layer
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
//has an affinity prediction layer
bool CNNScorer::has_affinity() const {
  return (bool) net->blob_by_name("predaff");
}

//reset center to be around ligand; reset receptor transformation
//call this before minimizing a ligand
void CNNScorer::set_center_from_model(model& m) {

  if (isfinite(current_center[0]) && cnnopts.moving_receptor()) {
    //when we recenter, we need to apply any receptor transformations to the ligand (inversed)
    //unless, of course, the center hasn't been set yet
    vec center = get_center();
    if (cnnopts.verbose) {
      std::cout << "CNN center " << center[0] << "," << center[1] << ","
          << center[2] << "\n";
      std::cout << "Rec transform ";
      m.rec_conf.print();
      std::cout << "\n";
    }
    vecv& coords = m.coordinates();
    gfloat3 c(center[0], center[1], center[2]);
    gfloat3 trans(-m.rec_conf.position[0], -m.rec_conf.position[1],
        -m.rec_conf.position[2]);
    qt rot = m.rec_conf.orientation.inverse();

    VINA_FOR_IN(i, coords) {
      gfloat3 coord = gfloat3(coords[i][0], coords[i][1], coords[i][2]) + trans-c; //inverse transformation
      gfloat3 pt = rot.rotate(coord.x,coord.y,coord.z);
      pt += c;
      coords[i][0] = pt.x;
      coords[i][1] = pt.y;
      coords[i][2] = pt.z;
    }
  }

  //reset protein
  m.rec_conf.position = vec(0, 0, 0);
  m.rec_conf.orientation = qt(1, 0, 0, 0);

  //calc center for ligand
  current_center = vec(0, 0, 0);
  for (auto coord : m.coordinates()) {
    current_center += coord;
  }
  current_center /= (float) m.coordinates().size();

  //mgrid is shared across threads, so do not set center except when guarded by mutex

  if (cnnopts.verbose) {
    std::cout << "new center: ";
    current_center.print(std::cout);
    std::cout << "\n";
  }
}


//populate score and aff with current network output
void CNNScorer::get_net_output(Dtype& score, Dtype& aff, Dtype& loss) {
  const caffe::shared_ptr<Blob<Dtype> > outblob = net->blob_by_name("output");
  const caffe::shared_ptr<Blob<Dtype> > lossblob = net->blob_by_name("loss");
  const caffe::shared_ptr<Blob<Dtype> > affblob = net->blob_by_name("predaff");

  const Dtype* out = outblob->cpu_data();
  score = out[1];
  aff = 0.0;
  if (affblob) {
    aff = affblob->cpu_data()[0];
  }

  loss = lossblob->cpu_data()[0];
}

// Extract ligand atoms and coordinates
std::pair<atomv, vecv> CNNScorer::get_ligand(model& m) const{
  atomv ligand_atoms = atomv(
    m.get_movable_atoms().cbegin() + m.ligands[0].node.begin,
    m.get_movable_atoms().cbegin() + m.m_num_movable_atoms
  );

  vecv ligand_coords = vecv(
    m.coordinates().cbegin() + m.ligands[0].node.begin,
    m.coordinates().cbegin() + m.m_num_movable_atoms
  );

  if(!cnnopts.flexopt){ // No flexible residues
    CHECK_EQ(ligand_atoms.size(), m.get_movable_atoms().size());
    CHECK_EQ(ligand_coords.size(), m.coordinates().size());
  }

  return std::pair<atomv, vecv>(ligand_atoms, ligand_coords);
}

std::pair<atomv,sz> CNNScorer::get_receptor(model& m, sz num_ligand_atoms) const{
  // Get receptor flexible atoms
  atomv flexible_atoms = atomv(
    m.get_movable_atoms().cbegin(),
    m.get_movable_atoms().cbegin() + m.ligands[0].node.begin
  );

  // Get inflex atoms
  atomv inflex_atoms = atomv(
    m.get_movable_atoms().cbegin() + m.m_num_movable_atoms,
    m.get_movable_atoms().cend()
  );

  // Check that ligand atoms and movable atoms in flexible residues sum up to the 
  // total number of movable atoms
  CHECK_EQ(flexible_atoms.size() + num_ligand_atoms, m.m_num_movable_atoms);

  if(!cnnopts.flexopt){ // No atoms from flexible residues should be present
    CHECK_EQ(flexible_atoms.size(), 0);
    CHECK_EQ(inflex_atoms.size(), 0);
  }

  // rmeli: Use std::move_iterator instead?
  atomv receptor_atoms;
  receptor_atoms.insert( // Insert flexible residues movable atoms
    receptor_atoms.end(), 
    flexible_atoms.cbegin(), 
    flexible_atoms.cend()
  );
  receptor_atoms.insert( // Insert inflex atoms
    receptor_atoms.end(), 
    inflex_atoms.cbegin(), 
    inflex_atoms.cend()
  );
  receptor_atoms.insert( // Insert fixed receptor atoms
    receptor_atoms.end(), 
    m.get_fixed_atoms().cbegin(), 
    m.get_fixed_atoms().cend()
  );

  if(!cnnopts.flexopt){ // No atoms from flexible residues should be present
    CHECK_EQ(receptor_atoms.size(), m.get_fixed_atoms().size());
  }

  return std::pair<atomv,sz>(receptor_atoms, flexible_atoms.size());
}

//return score of model, assumes receptor has not changed from initialization
//also sets affinity (if available) and loss (for use with minimization)
//if compute_gradient is set, also adds cnn atom gradient to m.minus_forces
//if maintain center, it will not reposition the molecule
//ALERT: clears minus forces
float CNNScorer::score(model& m, bool compute_gradient, float& affinity,
    float& loss) {
  boost::lock_guard<boost::recursive_mutex> guard(*mtx);
  if (!initialized()) return -1.0;

  caffe::Caffe::set_random_seed(cnnopts.seed); //same random rotations for each ligand..

  if (!isnan(cnnopts.cnn_center[0])) {
    mgrid->setGridCenter(cnnopts.cnn_center);
    current_center = mgrid->getGridCenter();
  } else if(!isnan(current_center[0])){
    mgrid->setGridCenter(current_center);
  }

  // rmeli: Carve out ligand from movable atoms
  atomv ligand_atoms;
  vecv ligand_coords;
  std::tie(ligand_atoms, ligand_coords) = get_ligand(m);

  mgrid->setLigand(ligand_atoms, ligand_coords, cnnopts.move_minimize_frame);

  // rmeli: Carve out receptor atoms from movable atoms
  atomv receptor_atoms;
  sz num_flexres_atoms;
  std::tie(receptor_atoms, num_flexres_atoms) = get_receptor(m, ligand_atoms.size());

  // rmeli: what about receptor coordinates?
  if (!cnnopts.move_minimize_frame) { //if fixed_receptor, rec_conf will be identify
    mgrid->setReceptor(receptor_atoms, m.rec_conf.position, m.rec_conf.orientation);
  } else { //don't move receptor
    mgrid->setReceptor(receptor_atoms);
    current_center = mgrid->getGridCenter(); //has been recalculated from ligand
    if(cnnopts.verbose) {
      std::cout << "current center: ";
      current_center.print(std::cout);
      std::cout << "\n";
    }
  }

  if(compute_gradient || cnnopts.outputxyz) {
    mgrid->enableLigandGradients();
    if(cnnopts.moving_receptor() || cnnopts.outputxyz){
      mgrid->enableReceptorGradients();
    }
    else if(cnnopts.flexopt){
      mgrid->enableReceptorGradients(); // rmeli: TODO flexres gradients only
    }
  }

  m.clear_minus_forces();
  double score = 0.0;
  affinity = 0.0;
  loss = 0.0;
  Dtype s = 0.0;
  Dtype a = 0.0;
  Dtype l = 0.0;

  unsigned cnt = 0;
  mgrid->setLabels(1); //for now pose optimization only
  for (unsigned r = 0, n = max(cnnopts.cnn_rotations, 1U); r < n; r++) {
    net->Forward(); //do all rotations at once if requested
    get_net_output(s, a, l);
    score += s;
    affinity += a;
    loss += l;

    if (cnnopts.cnn_rotations > 1) {
      std::cout << "RotateScore: " << s << "\n";
      if (a) std::cout << "RotateAff: " << a << "\n";
    }

    if (compute_gradient || cnnopts.outputxyz) {

      net->Backward();

      // Get ligand gradient
      mgrid->getLigandGradient(0, gradient);

      CHECK_EQ(gradient.size(), ligand_coords.size());

      // Get receptor gradient
      std::vector<gfloat3> gradient_rec;
      if (cnnopts.flexopt) { // Optimization of flexible residues
        mgrid->getReceptorGradient(0, gradient_rec);
        
        // rmeli: Select components of flexible residues!
        gradient_rec = std::vector<gfloat3>(
          gradient_rec.cbegin(), 
          gradient_rec.cbegin() + num_flexres_atoms
        );

        // rmeli: receptor coords?
        CHECK_EQ(gradient_rec.size(), num_flexres_atoms);
      }

      // Merge ligand and flexible residues gradient
      // Flexible residues comes first
      gradient.insert(gradient.begin(), gradient_rec.cbegin(), gradient_rec.cend());

      // Ligand plus flexible residues gradient
      m.add_minus_forces(gradient);

      // Gradient for rigid receptor transformation: translation and torque
      if(cnnopts.moving_receptor())
        mgrid->getReceptorTransformationGradient(0, m.rec_change.position, m.rec_change.orientation);
    }
    cnt++;
  }

  if (cnnopts.outputxyz) {
    const string& ligname = cnnopts.xyzprefix + "_lig.xyz";
    const string& recname = cnnopts.xyzprefix + "_rec.xyz";

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
    check_gradient();
  }

  if (cnnopts.outputdx) {
    //DANGER! This modifies the values in the network
    outputDX(m.get_name());
  }

  //if there were rotations, scale appropriately
  if (cnt > 1) {
    m.scale_minus_forces(1.0 / cnt);
  }
  affinity /= cnt;
  loss /= cnt;

  if (cnnopts.verbose)
    std::cout << std::fixed << std::setprecision(10) << "cnnscore "
        << score / cnt << "\n";

  return score / cnt;
}

//return only score
float CNNScorer::score(model& m) {
  float aff = 0;
  float loss = 0;
  return score(m, false, aff, loss); 
}

// To aid in debugging, will compute the gradient at the
// grid level, apply it with different multiples, and evaluate
// the effect. Perhaps may evaluate atom gradients as well?
//
// IMPORTANT: assumes model is already setup
// Prints out the results
void CNNScorer::check_gradient() {
  Dtype origscore = 0;
  Dtype origaff = 0;
  Dtype origloss = 0;
  Dtype newscore = 0.0;
  Dtype newaff = 0.0;
  Dtype newloss = 0.0;

  std::cout << std::scientific;
  Dtype lambda = 1.0;
  for (unsigned i = 0; i < 4; i++) {
    //score pose
    net->Forward();
    get_net_output(origscore, origaff, origloss);

    //backprop
    net->Backward();

    //get grid and diff blobs
    const caffe::shared_ptr<Blob<Dtype> > datablob = net->blob_by_name("data");
    Dtype *data = datablob->mutable_cpu_data();
    Dtype *diff = datablob->mutable_cpu_diff();

    //apply gradient
    caffe_cpu_axpby(datablob->count(), -lambda, diff, 1.0f, data); //sets data

    //propagate forward, starting _after_ molgrid
    net->ForwardFrom(1);

    //compare scores
    get_net_output(newscore, newaff, newloss);

    std::cout << "CHECK   " << origscore - newscore << "    LAMBDA: " << lambda
        << "   OLD: " << origscore << "," << origaff << "   NEW: " << newscore
        << "," << newaff << std::endl;

    //test a single channel
    unsigned channel = 16 + 15; //compile time constant
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
      if (*(data + off + i) > 0) //only modify if density actually exists
        *(data + off + i) += -lambda * *(diff + off + i);
    }
    //caffe_cpu_axpby(datablob->count(2), -lambda, diff+off, 1.0f, data+off);
    net->ForwardFrom(1);
    get_net_output(newscore, newaff, newloss);

    std::cout << "CHECKch " << origscore - newscore << " " << channel
        << " LAMBDA: " << lambda << "   OLD: " << origscore << "," << origaff
        << "   NEW: " << newscore << "," << newaff << std::endl;

    //super expensive - evaluate every grid point for channel
    net->Forward();
    get_net_output(origscore, origaff, origloss);
    for (unsigned i = 0; i < n; i++) {
      data = datablob->mutable_cpu_data();
      diff = datablob->mutable_cpu_diff();
      float gval = *(data + off + i);
      float gdiff = *(diff + off + i);
      if (gval > 0) {
        *(data + off + i) += lambda;
        net->ForwardFrom(1);
        get_net_output(newscore, newaff, newloss);
        std::cout << "GRIDch " << i << " gval " << gval << " gdiff " << gdiff
            << " lambda " << lambda << " change " << origscore - newscore
            << "\n";

        data = datablob->mutable_cpu_data();
        diff = datablob->mutable_cpu_diff();
        *(data + off + i) = max(gval - lambda, 0.0f);
        net->ForwardFrom(1);
        get_net_output(newscore, newaff, newloss);
        std::cout << "GRIDch " << i << " gval " << gval << " gdiff " << gdiff
            << " lambda " << -lambda << " change " << origscore - newscore
            << "\n";
        *(data + off + i) = gval;
      }
    }

    lambda /= 10.0;
  }

}

//dump dx files of the diff
//zero_values: run backward relevance with only dead node values
void CNNScorer::outputDX(const string& prefix, double scale, bool lrp,
    string layer_to_ignore, bool zero_values)

    {
  const caffe::shared_ptr<Blob<Dtype>> datablob = net->blob_by_name("data");

  const vector<caffe::shared_ptr<Layer<Dtype> > >& layers = net->layers();
//    if(datablob) {
  //this is a big more fragile than I would like.. if there is a pooling layer before
  //the first convoluational of fully connected layer and it is a max pooling layer,
  //change it to average before the backward to avoid a discontinuous map
  PoolingLayer<Dtype> *pool = NULL;
  for (unsigned i = 1, nl = layers.size(); i < nl; i++) {
    pool = dynamic_cast<PoolingLayer<Dtype>*>(layers[i].get());
    if (pool)
      break; //found it
    else
      if (layers[i]->type() == string("Convolution"))
        break; //give up
      else
        if (layers[i]->type() == string("InnerProduct")) break;
  }

  if (pool) {
    if (pool->pool() == PoolingParameter_PoolMethod_MAX) {
      pool->set_pool(PoolingParameter_PoolMethod_AVE);
    } else {
      pool = NULL; //no need to reset to max
    }
  }

  //must redo backwards with average pooling
  if (lrp) {
    net->Backward_relevance(layer_to_ignore, zero_values);
  } else
    net->Backward();

  string p = prefix;
  if (p.length() == 0) p = "dx";
  mgrid->dumpDiffDX(p, datablob.get(), scale);

  if (pool) {
    pool->set_pool(PoolingParameter_PoolMethod_MAX);
  }

}

void CNNScorer::outputXYZ(const string& base, const vector<gfloat3>& atoms,
    const vector<short>& whichGrid, const vector<gfloat3>& gradient) {
  const char* sym[] = { "C", "C", "C", "C", "Ca", "Fe", "Mg", "N", "N", "N",
      "N", "O", "O", "P", "S", "Zn", "C", "C", "C", "C", "Br", "Cl", "F", "N",
      "N", "N", "N", "O", "O", "O", "P", "S", "S", "I" };

  ofstream out(base.c_str());
  out.precision(5);

  // Count number of valid (>= 0) whichGrid entries
  size_t n_atoms{0};
  for(size_t i = 0, n = atoms.size(); i < n; ++i){
    if(whichGrid[i] >= 0){
      n_atoms++;
    }
  }
  out << n_atoms << "\n\n"; // xyz
  
  // Print coordinates and gradients
  for (unsigned i = 0, n = atoms.size(); i < n; ++i) {

    // Skip invalid channel
    if(whichGrid[i] < 0){
      continue;
    }

    out << sym[whichGrid[i]] << " ";
    out << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z << " ";
    out << gradient[i].x << " " << gradient[i].y << " " << gradient[i].z;
    if (i + 1 < n) out << "\n";
  }
}

