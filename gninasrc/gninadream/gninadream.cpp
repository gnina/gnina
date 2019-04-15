#include <boost/program_options.hpp>
#include <cmath>
#include <unistd.h>
#include <glob.h>
#include <regex>
#include "../lib/cnn_scorer.h"
#include <boost/filesystem>
#include <boost/algorithm/string.hpp>

#define CUDA_NUM_THREADS 512

typedef caffe::MolGridDataLayer<float> mgridT;

std::vector<std::string> glob(const std::string& pattern) {
    using namespace std;

    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // glob
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if(return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames
    vector<string> filenames;
    for(size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}

bool readDXGrid(istream& in, vec& center, double& res, float* grid, unsigned numgridpoints, 
    std::string& fname) {
  string line;
  vector<string> tokens;

  res = 0;
  getline(in, line);
  split(tokens, line, is_any_of(" \t"), token_compress_on);
  if (tokens.size() != 8) return false;
  unsigned n = lexical_cast<unsigned>(tokens[7]);
  if (lexical_cast<unsigned>(tokens[6]) != n) return false;
  if (lexical_cast<unsigned>(tokens[5]) != n) return false;

  //the center
  getline(in, line);
  split(tokens, line, is_any_of(" \t"), token_compress_on);
  if (tokens.size() != 4) return false;
  double x = lexical_cast<double>(tokens[1]);
  double y = lexical_cast<double>(tokens[2]);
  double z = lexical_cast<double>(tokens[3]);

  //the transformation matrix, which has the resolution
  getline(in, line);
  split(tokens, line, is_any_of(" \t"), token_compress_on);
  if (tokens.size() != 4) return false;
  res = lexical_cast<float>(tokens[1]);

  getline(in, line);
  split(tokens, line, is_any_of(" \t"), token_compress_on);
  if (tokens.size() != 4) return false;
  if (res != lexical_cast<float>(tokens[2])) return false;

  getline(in, line);
  split(tokens, line, is_any_of(" \t"), token_compress_on);
  if (tokens.size() != 4) return false;
  if (res != lexical_cast<float>(tokens[3])) return false;

  //figure out center
  double half = res * n / 2.0;
  center[0] = x + half;
  center[1] = y + half;
  center[2] = z + half;

  //grid connections
  getline(in, line);
  //object 3
  getline(in, line);

  unsigned total = 0;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < n; j++) {
      for (unsigned k = 0; k < n; k++) {
        in >> grid[((i * n) + j) * n + k]; 
        total++;
      }
    }
  }
  if (total != n * n * n) return false;
  if (total != numgridpoints) {
    std::cerr << "Number of grid points in file " << fname << " does not equal the number of grid points expected by trained net.\n" << std::endl;
    exit(1);
  }

  return true;
}

// device functions for warp-based reduction using shufl operations
// TODO: should probably just be factored out into gpu_math or gpu_util
template<class T>
__device__   __forceinline__ T warp_sum(T mySum) {
  for (int offset = warpSize >> 1; offset > 0; offset >>= 1)
    mySum += shuffle_down(mySum, offset);
  return mySum;
}

__device__ __forceinline__
bool isNotDiv32(unsigned int val) {
  return val & 31;
}

/* requires blockDim.x <= 1024, blockDim.y == 1 */
template<class T>
__device__   __forceinline__ T block_sum(T mySum) {
  const unsigned int lane = threadIdx.x & 31;
  const unsigned int wid = threadIdx.x >> 5;

  __shared__ T scratch[32];

  mySum = warp_sum(mySum);
  if (lane == 0) scratch[wid] = mySum;
  __syncthreads();

  if (wid == 0) {
    mySum = (threadIdx.x < blockDim.x >> 5) ? scratch[lane] : 0;
    mySum = warp_sum(mySum);
    if (threadIdx.x == 0 && isNotDiv32(blockDim.x))
      mySum += scratch[blockDim.x >> 5];
  }
  return mySum;
}


void cpu_l2(const float* optgrid, const float* screengrid, float* scoregrid, 
    size_t M, size_t N, size_t gsize) {
  // optimized grids
  for (size_t i=0; i<M; ++i) {
    // conformers to screen against
    for (size_t j=0; j<N; ++j) {
      float sum = 0.;
#pragma omp parallel for reduction(+:sum)
      for (size_t k=0; k<gsize; ++k) {
        float diff = optgrid[i * gsize + k] - screengrid[j * gsize + k];
        float sqdiff = diff * diff;
        sum += sqdiff;
      }
    scoregrid[i * N + j] = std::sqrt(sum);
    }
  }
}

__global__
void gpu_l2(const float* optgrid, const float* screengrid, float* scoregrid, 
    size_t M, size_t N, size_t gsize) {
  unsigned tidx = threadIdx.x;
  // optimized grids
  for (size_t i=0; i<M; ++i) {
    // conformers to screen against
    for (size_t j=0; j<N; ++j) {
      float sum = 0.;
      for (size_t k=0; k<gsize; k+=CUDA_NUM_THREADS) {
        float diff = optgrid[i * gsize + k] - screengrid[j * gsize + k];
        float sqdiff = diff * diff;
        sum += sqdiff;
      }
    float total = block_sum<float>(sum);
    if (tidx == 0)
      scoregrid[i * N + j] = sqrtf(total);
    }
  }
}

void do_exact_vs(std::shared_ptr<mgridT>& opt_mgrid, shared_ptr<Net<float> >& net, 
    std::vector<std::ostream>& out, bool gpu) {
  // use net top blob to do virtual screen against input sdf
  // produce output file for each input from which we started optimization
  // output will just be overlap score, in order of the compounds in the
  // original file
  // right now we assume these are pre-generated poses, although we could dock
  // them internally if desired (but would take forever)
  //
  // reinit MolGrid with params for virtual screening and a vanilla provider,
  // for each example rec is the lig used to set center and lig is one of the
  // vs ligands; we set use_rec_center and ignore_rec.  this is annoying
  // because done the naive way we have to regrid the same ligand many times
  //
  unsigned ncompounds;
  std::ifstream infile(vsfile.c_str());
  CHECK((bool)infile) << "Could not open " << vsfile;
  while (getline(infile, line))
  {
    ++ncompounds;
  }
  if (!ncompounds) throw usage_error("No compounds in virtual screen file");
  unsigned niters = ncompounds / default_batch_size;
  unsigned bs_lastiter = ncompounds % default_batch_size;

  MolGridDataParameter* mparam = opt_mgrid.mutable_molgrid_data_param();
  mparam->set_ignore_rec(true);
  mparam->set_use_rec_center(true);
  mparam->set_source(vsfile); // this is probably going to change
  unsigned batch_size = ncompounds > default_batch_size ? default_batch_size : ncompounds;
  mparam->set_batch_size = batch_size;
  // set up blobs for virtual screen compound grids and set up
  vector<Blob<Dtype>*> bottom; // will always be empty
  vector<Blob<Dtype>*> top(2); // want to use MGrid::Forward so we'll need a dummy labels blob
  opt_mgrid.VSLayerSetUp(bottom, top);
  shared_ptr<Blob<float> > scores(new Blob<float>());
  unsigned nopts = net->top_vecs()[0][0].shape()[0];
  unsigned example_size = opt_mgrid.getExampleSize();
  std::vector<int> score_shape = {nopts};
  scores->Reshape(score_shape);
  float* scoregrid;

  for (size_t i=0; i<niters; ++i) {
    if (i==niters-1 && batch_size > bs_lastiter)  {
      mparam->set_batch_size = bs_lastiter;
      auto data_shape = top[0].shape();
      data_shape[0] = bs_lastiter;
      vector<int> label_shape;
      label_shape.push_back(bs_lastiter);
      top[0].Reshape(data_shape);
      top[1].Reshape(label_shape);
    }
    if (gpu) {
      opt_mgrid.Forward_gpu(bottom, top);
      float* optgrid = net->top_vecs()[0][0]->gpu_data();
      float* screengrid = top[0]->gpu_data();
      scoregrid = scores->mutable_gpu_data();
      gpu_l2<<<1, CUDA_NUM_THREADS>>>(optgrid, screengrid, scoregrid, nopts, batch_size, example_size);
    }
    else {
      opt_mgrid.Forward_cpu(bottom, top);
      float* optgrid = net->top_vecs()[0][0]->cpu_data();
      float* screengrid = top[0]->cpu_data();
      scoregrid = scores->mutable_cpu_data();
      cpu_l2(optgrid, screengrid, scoregrid, nopts, batch_size, example_size);
    }
    
    // write to output
    scoregrid = scores->mutable_cpu_data();
    for (size_t j=0; j<nopts; ++j) {
      float* datastart = scoregrid + j * batch_size;
      for (size_t k=0; k<batch_size; ++k) {
        out[j] << *(datastart + k);
      }
    }
  }
}

void do_approx_vs(mgridT& opt_mgrid, mgridT& vs_mgrid) {
  // TODO?
  assert(0);
}

int main(int argc, char* argv[]) {
  using namespace boost::program_options;

  std::string receptor_name, ligand_names, grid_prefix, vsfile, solverstate, 
    outname, out_prefix, types;

  cnn_options cnnopts;
  cnnopts.cnn_scoring = True;
  int iterations;
  int gpu;
  float base_lr;
  bool dump_all = false;
  bool dump_last = false;
  bool exclude_receptor = false;
  bool exclude_ligand = false;
  bool ignore_ligand = false;
  bool allow_neg = false;
  std::string sigint_effect = "stop";
  std::string sighup_effect = "snapshot";
  unsigned default_batch_size = 50;
  unsigned nopts = 0;
  std::vector<std::string> opt_names;

  options_description inputs("General Input");
  inputs.add_options()("receptor, r",
      value<std::string>(&receptor_name), "receptor to provide optimization context")(
      "ligand, l", value<std::string>(&ligand_names),
      "one or more ligands to provide additional optimization context")(
      "types, t", value<std::string>(&types), 
      "MolGrid .types file, formatted with one example per line, <score> <affinity> <receptor_file> <ligand_file>")(
      "grid, g", value<std::string>(&grid_prefix), 
      "prefix for grid files from which to begin optimization (instead of molecules), filenames assumed to be [prefix]_[Rec/Lig]_[channel0]_[channel1][...].dx")(
      "virtualscreen, vs", value<std::string>(&vsfile), 
      "file of compounds to score according to overlap with optimized grid");

  options_description cnn("CNN Input");
  cnn.add_options()("cnn_model", value<std::string>(&cnnopts.cnn_model),
      "CNN model file (*.model), if not provided current default model will be used")(
      "cnn_model_name", value<std::string>(&cnnopts.cnn_model_name), 
      ("built-in model to use: " + builtin_cnn_models()).c_str())(
      "cnn_weights", value<std::string>(&cnnopts.cnn_weights),
      "CNN weights file (*.caffemodel), if not provided current default weights will be used; if no other network input is provided an attempt will be made to restart optimization from this")(
      "solverstate", value<std::string>(&solverstate),
      "CNN solverstate for restarting optimization in-progress");

  options_description output("Output");
  output.add_options()("vs_output, o",
      value<std::string>(&outname), "output virtual screen compounds with grid overlap score")(
      "out_prefix, p", value<std::string>(&out_prefix),
      "prefix for storing checkpoint files for optimization in progress, default is gninadream.<PID>")(
      "dump_all, a", bool_switch(&dump_all)->default_value(false), 
      "dump all intermediate grids from optimization")(
      "dump_last, dl", bool_switch(&dump_last)->default_value(false),
      "dump last grid from optimization");

  options_description options("Options");
  options.add_options()("iterations, i",
      value<int>(&iterations), "number of iterations to run")(
      "base_lr", value<float>(&base_lr),
      "base learning rate for density updates")(
      "gpu, g", value<int>(&gpu)->default_value(-1), "gpu to run on")(
      "exclude_receptor, er", bool_switch(&exclude_receptor)->default_value(false), 
      "don't update the receptor grids")(
      "exclude_ligand, el",  bool_switch(&exclude_ligand)->default_value(false), 
      "don't update the ligand grids")(
      "ignore_ligand, il", bool_switch(&ignore_ligand)->default_value(false),
      "just use ligand to set center")(
      "allow_negative, an", bool_switch(&allow_neg)->default_value(false),
      "allow optimization to result in negative atom density");

  options_description desc;
  desc.add(inputs).add(cnn).add(output);

  positional_options_description positional; 
  variables_map vm;
  try {
    store(

        command_line_parser(argc, argv).options(desc).style(
            command_line_style::default_style
                ^ command_line_style::allow_guessing).positional(positional).run(),
        vm);
    notify(vm);
  } catch (boost::program_options::error& e) {
    std::cerr << "Command line parse error: " << e.what() << '\n'
        << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  //treat rec + lig + grids as mutually exclusive with types (the former are
  //"in_mem" for MolGrid)
  //an existing solverstate is used only if none of the other options are set
  //if starting from a rec+lig or solverstate, we optimize for requested number
  //of iterations and terminate; if we have a .types file we optimize for that
  //number of iterations for each example in the batch (batching so that we
  //don't run out of memory)
  if (vm.count("receptor") <= 0 && vm.count("ligand") <= 0 && vm.count("grid") <= 0 && 
      vm.count("types") <= 0 && vm.count("solverstate") <= 0) {
    std::cerr << "Missing optimization context.\n" << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  if (!out_prefix.size()) 
    out_prefix = "gninadream." + std::to_string(getpid()) + ".out";

  google::InitGoogleLogging(argv[0]);
  google::SetStderrLogging(2);

  //find out/set up some stuff about the CNN model
  NetParameter net_param;

  if (cnnopts.cnn_model.size() == 0) {
    if(cnn_models.count(cnnopts.cnn_model_name) == 0) {
      throw usage_error("Invalid model name: "+cnnopts.cnn_model_name);
    }

    const char *model = cnn_models[cnnopts.cnn_model_name].model;
    google::protobuf::io::ArrayInputStream modeldata(model, strlen(model));
    bool success = google::protobuf::TextFormat::Parse(&modeldata, &net_param);
    if (!success) throw usage_error("Error with built-in cnn model "+cnnopts.cnn_model_name);
    UpgradeNetAsNeeded("default", &net_param);
  } else {
    ReadNetParamsFromTextFileOrDie(cnnopts.cnn_model, &net_param);
  }

  net_param.mutable_state()->set_phase(TRAIN);

  LayerParameter *first = net_param.mutable_layer(0);
  MolGridDataParameter* mgridparam = first->mutable_molgrid_data_net_param();
  if (mgridparam == NULL) {
    throw usage_error("First layer of model must be MolGridData.");
  }

  if (cnnopts.cnn_model.size() == 0) {
    const char *recmap = cnn_models[cnnopts.cnn_model_name].recmap;
    const char *ligmap = cnn_models[cnnopts.cnn_model_name].ligmap;
    mgridparam->set_mem_recmap(recmap);
    mgridparam->set_mem_ligmap(ligmap);
  }

  // if we have structure file(s) or grids as input, we'll do just one example,
  // inmem=true
  if (receptor.size() || grid_prefix.size()) {
    mgridparam->set_inmemory(true);
    mgridparam.set_batch_size(1);
  }
  else {
    batch_size = mgridparam.batch_size();
    if (types.size())
      mgridparam->set_source(types.c_str());
  }

  net_param.set_force_backward(true);
  mgridparam.set_ignore_ligand(ignore_ligand);

  //set up solver params and then construct solver
  if (gpu > -1) {
    caffe::Caffe::SetDevice(gpu);
    caffe::Caffe::set_mode(caffe::Caffe::GPU);
    //TODO: should be possible to toggle pool in cudnn layers, just requires 
    //different syntax, maybe it doesn't matter though
    caffe::Caffe::set_cudnn(false);
  }

  caffe::SignalHandler signal_handler(
        GetRequestedAction(sigint_effect),
        GetRequestedAction(sighup_effect));

  caffe::SolverParameter solver_param;
  solver_param.net_param = net_param;
  solver_param.base_lr = base_lr;
  solver_param.max_iter = iterations;
  solver_param.lr_policy = fixed;
  solver_param.snapshot_prefix = "inputopt_";
  solver_param.type = "InputOptSGD";
  if (cnnopts.cnn_weights.size())
    solver_param.weights = cnnopts.cnn_weights;

  shared_ptr<caffe::Solver<float> >
      solver(caffe::SolverRegistry<float>::CreateSolver(solver_param));

  solver->SetActionFunction(signal_handler.GetActionFunction());
  shared_ptr<Net<float> > net = solver->net();
  //if there wasn't a weights file, check that we're using one of the provided
  //cnn models and attempt to load the appropriate stored weights
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

  const vector<caffe::shared_ptr<Layer<float> > >& layers = net->layers();
  mgridT mgrid = dynamic_cast<MolGridDataLayer<float>*>(layers[0].get());
  if (mgrid == NULL) {
    throw usage_error("First layer of model must be MolGridDataLayer.");
  }

  // figure out what the initial state should be and run optimization
  if (receptor.size()) {
    // can start from receptor and 0+ ligands
    // expected behavior is starting from rec with a ligand available for
    // setting center at minimum; without ligand things _should_ work (verify
    // there aren't errors about not having lig atoms) and expect that grid
    // center will be set to origin in that case. possibly that should
    // be changed
    MolGetter mols(receptor);
    ++nopts;
    if (!ligand_names.size()) {
      ligand_names.push_back("");
    }
    for (unsigned l = 0, nl = ligand_names.size(); l < nl; l++) {
      const std::string ligand_name = ligand_names[l];
      mols.setInputFile(ligand_name);
      for (;;)  {
        model* m = new model;

        if (!mols.readMoleculeIntoModel(*m)) {
          delete m;
          break;
        }
        mgrid->setLigand(m.get_movable_atoms(), m.coordinates());
        mgrid->setReceptor(m.get_fixed_atoms());
        solver.ResetIter();
        solver.Solve();
        boost::filesystem::path rec(receptor);
        boost::filesystem::path lig(ligand_name);
        opt_names.push_back(rec.stem() + "_" + lig.stem());
        ++nopts;
      }
    }
  }
  else if (types.size()) {
    //molgrid will load in batches from types
    unsigned n_passes = 1;
    std::ifstream infile(types.c_str());
    CHECK((bool)infile) << "Could not open " << types;
    while (getline(infile, line))
    {
      ++nopts;
      std::vector<std::string> contents;
      boost::split(contents, line, boost::is_any_of(" "));
      size_t size = contents.size();
      boost::filesystem::path rec(contents[size-1]);
      boost::filesystem::path lig(contents[size-2]);
      opt_names.push_back(rec.stem() + "_" + lig.stem());
    }
    if (!nopts) throw usage_error("No examples in types file");
    n_passes = std::ceil((float)(nopts) / batch_size);
    for (unsigned i=0; i<n_passes; ++i) {
      net->ForwardFromTo(1,1);
      solver.ResetIter();
      solver.Solve();
      if (vsfile.size()) {
        std::vector<ostream> out(nopts);
        for (size_t j=0; j<nopts; ++j) {
          out.push_back(std::ofstream((opt_names[j] + ".vsout").c_str()));
        }
        //use the resulting grids to perform a virtual screen against these
        //compounds; if there were multiple inputs optimized, return a separate
        //preds file for each of them 
        //right now we only support an L2 loss on the grid; generate grids for the
        //virtual screen compounds centered on whatever the grid center for the
        //inputopt grids was
      }
    }
  }
  else if (grid_prefix.size()) {
    std::vector<std::string> filenames = glob(grid_prefix + "*.dx");
    /* 
     * want a vector of unordered sets for the maps
     */
    std::vector<std::string> rectypes = mgrid->getRecTypes();
    std::vector<std::unordered_set<std::string> > recset;
    for (auto& substr : rectypes) {
      // these look like Rec_CHANNEL1_CHANNEL2_...
      std::string chs = std::regex_replace(substr, std::regex("Rec_"), "");
      std::vector<std::string> subchs;
      boost::split(subchs, chs, boost::is_any_of("_"));
      recset.push_back(std::unordered_set<std::string>(&subchs[0], &subchs[subchs.size()]));
    }

    std::vector<std::string> ligtypes = mgrid->getLigTypes();
    std::vector<std::unordered_set<std::string> > ligset;
    for (auto& substr : ligtypes) {
      // these look like Lig_CHANNEL1_CHANNEL2_...
      std::string chs = std::regex_replace(substr, std::regex("Lig_"), "");
      std::vector<std::string> subchs;
      boost::split(subchs, chs, boost::is_any_of("_"));
      ligset.push_back(std::unordered_set<std::string>(&subchs[0], &subchs[subchs.size()]));
    }
    /* 
     * we need to figure out which channels/groups of channels make up the 
     * maps used to generate these grids and compare with the trained net 
     * we're using. if they don't match, error out, otherwise load in the grids
     * in the right order for the map
     */ 
    float* inputblob;
    if (gpu)  {
      inputblob = net->top_vecs()[0][0]->mutable_gpu_data();
      CUDA_CHECK_GNINA(cudaMemset(inputblob, 0, numgridpoints * sizeof(float)));
    }
    else {
      inputblob = net->top_vecs()[0][0]->mutable_cpu_data();
      memset(inputblob, 0, numgridpoints * sizeof(float));
    }
    double resolution = mgrid->getResolution();
    unsigned npts_allchs = mgrid->getNumGridPoints();
    unsigned nchannels = numgridpoints / (rectypes.size() + ligtypes.size());
    unsigned npts_onech = npts_allchs / nchannels;

    for (auto& fname : filenames) {
      boost::filesystem::path p(fname);
      std::string stem = p.stem();
      std::string channels = std::regex_replace(stem, std::regex(grid_prefix + "_"), "");
      std::vector<std::string> channelvec;
      boost::split(channelvec, channels, boost::is_any_of("_"));
      // the first entry in channelvec is either Rec_ or Lig_
      std::unordered_set<std::string> channelmap(&channelvec[1], &channelvec[channelvec.size()]);
      for (auto& thisch : channelmap) {
        // does this correspond to channels in the net's map, if so what's the index
        // and use that to copy it
        unsigned idx = 0;
        if (channelvec[0] == "Rec_") {
          for (auto& chanset : recset) {
            if (channelmap == chanset) {
              // cool it's a match read it in
		          ifstream gridfile(fname.c_str());
		          if(!gridfile)
		          {
                std::cerr << "Could not open example grid file " << fname << "\n";
		          	abort();
		          }
		          vec gridcenter;
		          double gridres = 0;
              // if we're running on the CPU we can read directly into the
              // right place, otherwise we read into a temporary host buffer
              // and memcpy to the device
              float* g;
              std::vector<float> hostbuf;
              if (gpu) {
                hostbuf.resize(npts_onech);
                g = &hostbuf[0];
              }
              else 
                g = inputblob + idx * numgridpoints;
		          if(!readDXGrid(gridfile, gridcenter, gridres, g, npts_onech, fname))
		          {
                std::cerr << "I couldn't understand the provided dx file " << fname << ". I apologize for not being more informative and possibly being too picky about my file formats.\n";
		          	exit(1);
		          }
              if (gridres - resolution > 0.001) {
                std::cerr << "grid resolution in file " << fname << " does not match grid resolution specified for trained net.\n" << std::endl;
                exit(1);
              }
              if (gpu) 
                CUDA_CHECK_GNINA(cudaMemcpy(g, inputblob + idx * npts_onech, 
                      npts_onech * sizeof(float), cudaMemcpyHostToDevice));
              break;
            }
            ++idx;
          }
        }
        else if (channelvec[0] == "Lig_") {
          idx = rectypes.size();
          for (auto& chanset : ligset) {
            if (channelmap == chanset) {
              // cool it's a match read it in
		          ifstream gridfile(fname.c_str());
		          if(!gridfile)
		          {
                std::cerr << "Could not open example grid file " << fname << "\n";
		          	abort();
		          }
		          vec gridcenter;
		          double gridres = 0;
              // if we're running on the CPU we can read directly into the
              // right place, otherwise we read into a temporary host buffer
              // and memcpy to the device
              float* g;
              std::vector<float> hostbuf;
              if (gpu) {
                hostbuf.resize(npts_onech);
                g = &hostbuf[0];
              }
              else 
                g = inputblob + idx * numgridpoints;
		          if(!readDXGrid(gridfile, gridcenter, gridres, g, npts_onech, fname))
		          {
                std::cerr << "I couldn't understand the provided dx file " << fname << ". I apologize for not being more informative and possibly being too picky about my file formats.\n";
		          	exit(1);
		          }
              if (gridres - resolution > 0.001) {
                std::cerr << "grid resolution in file " << fname << " does not match grid resolution specified for trained net.\n" << std::endl;
                exit(1);
              }
              if (gpu) 
                CUDA_CHECK_GNINA(cudaMemcpy(g, inputblob + idx * npts_onech, 
                      npts_onech * sizeof(float), cudaMemcpyHostToDevice));
              break;
            }
            ++idx;
          }
        }
        else {
          std::cerr << "Grid files don't match pattern PREFIX_(Rec|Lig)_CHANNELS.\n";
          exit(1);
        }
      }
    }
  }
  else if (solverstate.size()) {
    char* ss_cstr = solverstate.c_str();
    //restart from optimization in progress
    // FIXME: do we have atom type info?
    solver->Restore(ss_cstr);
    solver.Solve();
    nopts = 1;
    // just use out_prefix as outname for vs?
  }
  else {
    std::cerr << "No valid initial input for optimization provided.\n";
    exit(1);
  }
}
