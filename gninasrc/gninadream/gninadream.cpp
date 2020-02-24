#include <boost/program_options.hpp>
#include <cmath>
#include <unistd.h>
#include <cstdio>
#include <glob.h>
#include <regex>
#include <unordered_set>
#include <boost/filesystem.hpp>
#include <libmolgrid/grid_io.h>
#include "caffe/util/signal_handler.h"
#include "gninadream.h"

struct atom_info {
	float x, y, z;
	int type;

	atom_info(): x(0), y(0), z(0), type(-1) {}
	atom_info(float X, float Y, float Z, int T): x(X), y(Y), z(Z), type(T) {}
};

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

    globfree(&glob_result);
    return filenames;
}

// Translate the signal effect the user specified on the command-line to the
// corresponding enumeration.
caffe::SolverAction::Enum GetRequestedAction(
    const std::string& flag_value) {
  if (flag_value == "stop") {
    return caffe::SolverAction::STOP;
  }
  if (flag_value == "snapshot") {
    return caffe::SolverAction::SNAPSHOT;
  }
  if (flag_value == "none") {
    return caffe::SolverAction::NONE;
  }
  LOG(FATAL) << "Invalid signal effect \""<< flag_value << "\" was specified";
}

bool readDXGrid(std::istream& in, vec& center, double& res, float* grid, unsigned numgridpoints, 
    std::string& fname) {
  string line;
  vector<string> tokens;

  res = 0;
  getline(in, line);
  boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
  if (tokens.size() != 8) return false;
  unsigned n = boost::lexical_cast<unsigned>(tokens[7]);
  if (boost::lexical_cast<unsigned>(tokens[6]) != n) return false;
  if (boost::lexical_cast<unsigned>(tokens[5]) != n) return false;

  //the center
  getline(in, line);
  boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
  if (tokens.size() != 4) return false;
  double x = boost::lexical_cast<double>(tokens[1]);
  double y = boost::lexical_cast<double>(tokens[2]);
  double z = boost::lexical_cast<double>(tokens[3]);

  //the transformation matrix, which has the resolution
  getline(in, line);
  boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
  if (tokens.size() != 4) return false;
  res = boost::lexical_cast<float>(tokens[1]);

  getline(in, line);
  boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
  if (tokens.size() != 4) return false;
  if (res != boost::lexical_cast<float>(tokens[2])) return false;

  getline(in, line);
  boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);
  if (tokens.size() != 4) return false;
  if (res != boost::lexical_cast<float>(tokens[3])) return false;

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
    std::exit(1);
  }

  return true;
}

void write_dx(const std::string& prefix, mgridT* mgrid, 
    const libmolgrid::Grid<float,4>& grid) {
  std::vector<std::string> recnames = mgrid->getRecTypes();
  gfloat3 center(mgrid->getGridCenter());
  double resolution = mgrid->getResolution();
  for (unsigned a = 0, na = recnames.size(); a < na; ++a) {
    std::string name = recnames[a];
    std::string fname = prefix + "_rec_" + name + ".dx";
    std::ofstream out(fname.c_str());
    libmolgrid::write_dx<float>(out, grid[a], center, resolution);
  }
  std::vector<std::string> lignames = mgrid->getLigTypes();
  unsigned offset = recnames.size();
  for (unsigned a = 0, na = lignames.size(); a < na; ++a) {
    std::string name = lignames[a];
    std::string fname = prefix + "_lig_" + name + ".dx";
    std::ofstream out(fname.c_str());
    libmolgrid::write_dx<float>(out, grid[a+offset], center, resolution);
  }
}

void cpu_constant_fill(float* fillgrid, size_t gsize, float fillval) {
#pragma omp parallel
  for (size_t i=0; i<gsize; ++i) {
    if (fillgrid[i] == 0) 
      fillgrid[i] = fillval;
  }
}

void setLigand(model& m, std::vector<float3>& coords, std::vector<smt>& smtypes) {
  // Ligand atoms start at ligand root node start idx
  sz n = m.m_num_movable_atoms - m.ligands[0].node.begin;

  auto m_atoms = m.get_movable_atoms().cbegin() + m.ligands[0].node.begin;
  auto m_coords = m.coordinates().cbegin() + m.ligands[0].node.begin;

  for(sz i=0;  i<n; ++i) {
    smtypes.push_back(m_atoms[i].sm);
    const vec& coord = m_coords[i];
    coords.push_back(gfloat3(coord[0], coord[1], coord[2]));
  }
}

void setReceptor(model& m, std::vector<float3>& coords, std::vector<smt>& smtypes) {
  auto& at = m.get_fixed_atoms();
  sz n = at.size();
  for (sz i=0; i<n; ++i) {
    const vec& coord = at[i].coords;
    coords.push_back(make_float3(coord[0], coord[1], coord[2]));
    smtypes.push_back(at[i].sm);
  }
}

void do_exact_vs(LayerParameter param, caffe::Net<float>& net, 
    std::string vsfile, std::vector<std::string>& ref_ligs,
    std::vector<caffe::shared_ptr<std::ostream> >& out, 
    bool gpu, std::string dist_method, float positive_threshold, float negative_threshold, 
    bool compute_cost) {
  // use net top blob to do virtual screen against input sdf
  // produce output file for each input from which we started optimization
  // output will just be overlap score, in order of the compounds in the
  // original file
  // right now we assume these are pre-generated poses, although we could dock
  // them internally or generate conformers in theory. 
  //
  // reinit MolGrid with params for virtual screening
  // for each example rec is the lig used to set center (unless there was no
  // lig, in which case we effectively fix center to origin) and lig is one of the
  // vs ligands; we set use_rec_center and ignore_rec if there was an autocenter
  // lig. this is annoying because done the naive way we have to regrid the
  // same ligand many times
  unsigned nopts = net.top_vecs()[0][0]->shape()[0];
  unsigned batch_size = 1; // inmem for now

  tee log(true);
  unsigned subgrid_dim = 4; // right now require grid dimension be divisible by this if using EMD
  FlexInfo finfo(log);
  MolGetter mols(std::string(), std::string(), finfo, false, false, log);
  mols.setInputFile(vsfile);

  MolGridDataParameter* mparam = param.mutable_molgrid_data_param();
  if (!mparam) {
    std::cerr << "Virtual screen passed non-molgrid layer parameter.\n";
    std::exit(1);
  }
  mparam->set_ignore_rec(true);
  mparam->set_ignore_ligand(false);
  mparam->set_has_affinity(false);
  mparam->set_inmemory(true);
  mparam->set_use_rec_center(true);
  mparam->set_batch_size(batch_size);
  // initblobs for virtual screen compound grids and set up
  vector<Blob<float>*> bottom; // will always be empty
  vector<Blob<float>*> top(2); // want to use MGrid::Forward so we'll need a dummy labels blob
  Blob<float> datablob;
  Blob<float> labelsblob;
  top[0] = &datablob;
  top[1] = &labelsblob;

  mgridT opt_mgrid(param);
  opt_mgrid.DataLayerSetUp(bottom, top);
  // actually want size for lig channels only
  unsigned example_size = opt_mgrid.getExampleSize();
  unsigned ligGridSize = opt_mgrid.getNumGridPoints() * opt_mgrid.getNumLigandTypes();
  unsigned recGridSize = opt_mgrid.getNumGridPoints()* opt_mgrid.getNumReceptorTypes();
  float* gpu_score = nullptr;
  // two floats since for the sum method we need two storage locations
  if (gpu) 
    CUDA_CHECK_GNINA(cudaMalloc(&gpu_score, sizeof(float)*2));
  
  // if using EMD for virtual screen, populate cost_matrix once
  flmat cost_matrix;
  unsigned blocks_per_side;
  unsigned ncubes;
  unsigned ntypes;
  unsigned n_matrix_elems;
  unsigned dim;
  if (!std::strcmp(dist_method.c_str(), "emd") && compute_cost) {
    double dimension = opt_mgrid.getDimension();
    dim = dimension * 2 + 1;
    if (dim % subgrid_dim) {
      std::cerr << "Grid dimension not divisible by 4, which is currently required for EMD.\n";
      std::exit(1);
    }
    blocks_per_side = dim / subgrid_dim;
    ncubes = blocks_per_side * blocks_per_side * blocks_per_side;
    ntypes = opt_mgrid.getNumLigandTypes();
    n_matrix_elems = ncubes * ntypes;
    cost_matrix = flmat(n_matrix_elems, 0);
    populate_cost_matrix(ncubes, subgrid_dim, ntypes, blocks_per_side, dimension, cost_matrix);
  }

  //VS compounds are currently done one at a time, inmem
  //out is the vector of output filestreams, one per optimized input to be
  //screened against 
  model m;
  for (size_t i=0; i<out.size(); ++i) {
    std::vector<float> scores;
    unsigned offset = i * example_size + recGridSize;
    std::string& next_ref_lig = ref_ligs[i];
    bool is_gninatypes = boost::algorithm::ends_with(next_ref_lig, ".gninatypes");
    
    if (next_ref_lig != "none" && !is_gninatypes)
      mols.create_init_model(ref_ligs[0], std::string(), finfo, log);
    else
      mols.create_init_model(std::string(), std::string(), finfo, log);
    for (;;) {
      if (!mols.readMoleculeIntoModel(m)) {
        break;
      }
      unsigned natoms = m.num_movable_atoms();
      std::vector<smt> ligand_smtypes;
      std::vector<float3> ligand_coords;
      setLigand(m, ligand_coords, ligand_smtypes);
      opt_mgrid.setLigand(ligand_coords, ligand_smtypes);
      if (next_ref_lig != "none") {
        if (is_gninatypes) {
          struct info {
            float x,y,z;
            int type;
          } gt_atom;
          std::vector<float3> rec_coords;
          std::vector<smt> rec_smtypes;

          std::ifstream in(next_ref_lig.c_str());
          CHECK(in) << "Could not read " << next_ref_lig;

          while(in.read((char*)&gt_atom, sizeof(gt_atom)))
          {
            smt t = (smt)gt_atom.type;
            rec_smtypes.push_back(t);
            rec_coords.push_back(make_float3(gt_atom.x, gt_atom.y, gt_atom.z));

          }
          opt_mgrid.setReceptor(rec_coords, rec_smtypes);
        }
        else {
          std::vector<smt> rec_smtypes;
          std::vector<float3> rec_coords;
          setReceptor(m, rec_coords, rec_smtypes);
          opt_mgrid.setReceptor(rec_coords, rec_smtypes);
        }
      }
      else
        opt_mgrid.setGridCenter(vec(0, 0, 0));
      opt_mgrid.setLabels(1, 10); 
      if (gpu) {
        opt_mgrid.Forward_gpu(bottom, top);
        const float* optgrid = net.top_vecs()[0][0]->gpu_data();
        const float* screengrid = top[0]->gpu_data();
        CUDA_CHECK_GNINA(cudaMemset(gpu_score, 0, sizeof(float)*2));
        if (!std::strcmp(dist_method.c_str(), "l1"))
          do_gpu_l1(optgrid + offset, screengrid + recGridSize, gpu_score, ligGridSize);
        else if (!std::strcmp(dist_method.c_str(), "l2"))
          do_gpu_l2sq(optgrid + offset, screengrid + recGridSize, gpu_score, ligGridSize);
        else if(!std::strcmp(dist_method.c_str(), "mult"))
          do_gpu_mult(optgrid + offset, screengrid + recGridSize, gpu_score, ligGridSize);
        else if(!std::strcmp(dist_method.c_str(), "sum")) {
          float l2;
          float mult;
          do_gpu_l2sq(optgrid + offset, screengrid + recGridSize, gpu_score, ligGridSize);
          CUDA_CHECK_GNINA(cudaMemcpy(&l2, gpu_score, sizeof(float), cudaMemcpyDeviceToHost));
          do_gpu_mult(optgrid + offset, screengrid + recGridSize, gpu_score+1, ligGridSize);
          CUDA_CHECK_GNINA(cudaMemcpy(&mult, gpu_score+1, sizeof(float), cudaMemcpyDeviceToHost));
          l2 = std::sqrt(l2);
          scores.push_back((l2/100.f + mult) / (natoms * ligGridSize));
        }
        else if(!std::strcmp(dist_method.c_str(), "threshold")) {
          do_gpu_thresh(optgrid + offset, screengrid + recGridSize, gpu_score, ligGridSize,
              positive_threshold, negative_threshold);
        }
        else if (!std::strcmp(dist_method.c_str(), "emd")) {
          do_gpu_emd(optgrid + offset, screengrid + recGridSize, gpu_score, dim, subgrid_dim, 
              blocks_per_side, ntypes, ligGridSize, cost_matrix);
        }
        else {
          std::cerr << "Unknown distance method for overlap-based virtual screen\n";
          exit(-1);
        }
        if(std::strcmp(dist_method.c_str(), "sum")) {
          float scoresq;
          CUDA_CHECK_GNINA(cudaMemcpy(&scoresq, gpu_score, sizeof(float), cudaMemcpyDeviceToHost));
          if (!std::strcmp(dist_method.c_str(), "l2"))
            scoresq = std::sqrt(scoresq);
          scores.push_back(scoresq / (natoms * ligGridSize));
        }
      }
      else {
        opt_mgrid.Forward_cpu(bottom, top);
        const float* optgrid = net.top_vecs()[0][0]->cpu_data();
        const float* screengrid = top[0]->cpu_data();
        scores.push_back(float());
        if (!std::strcmp(dist_method.c_str(), "l1"))
          cpu_l1(optgrid + offset, screengrid + recGridSize, &scores.back(), 
              ligGridSize);
        else if (!std::strcmp(dist_method.c_str(), "l2"))
          cpu_l2sq(optgrid + offset, screengrid + recGridSize, &scores.back(), 
              ligGridSize);
        else if(!std::strcmp(dist_method.c_str(), "mult"))
          cpu_mult(optgrid + offset, screengrid + recGridSize, &scores.back(), 
              ligGridSize);
        else if(!std::strcmp(dist_method.c_str(), "sum")) {
          float l2;
          cpu_l2sq(optgrid + offset, screengrid + recGridSize, &l2, 
              ligGridSize);
          float mult;
          cpu_mult(optgrid + offset, screengrid + recGridSize, &mult, 
              ligGridSize);
          scores.back() = (l2/100.f + mult) / (natoms * ligGridSize);
        }
        else if(!std::strcmp(dist_method.c_str(), "threshold")) {
          cpu_thresh(optgrid + offset, screengrid + recGridSize, &scores.back(), 
              ligGridSize, positive_threshold, negative_threshold);
        }
        else if (!std::strcmp(dist_method.c_str(), "emd")) {
          do_cpu_emd(optgrid+offset, screengrid + recGridSize, &scores.back(), dim, subgrid_dim, 
              blocks_per_side, ntypes, ligGridSize, cost_matrix);
        }
        else {
          std::cerr << "Unknown distance method for overlap-based virtual screen\n";
          exit(-1);
        }
        scores.back() = std::sqrt(scores.back()) / (natoms * ligGridSize);
      }
    }
    // write to output
    for (auto& score : scores) {
      *out[i] << score << "\n";
    }
  }
}

void do_approx_vs(mgridT* opt_mgrid, caffe::Net<float>& net, 
    std::string vsfile, std::vector<std::string>& ref_ligs,std::vector<std::ostream>& out, 
    bool gpu) {
  // TODO?
  assert(0);
}

void carve_init(mgridT& mgrid, bool no_density, caffe::shared_ptr<Net<float> >& net, 
    float carve_val) {
  unsigned nlt = mgrid.getNumLigandTypes();
  unsigned nrt = mgrid.getNumReceptorTypes();
  unsigned npts = mgrid.getNumGridPoints();
  unsigned lig_grid_size = nlt * npts;
  unsigned rec_grid_size = nrt * npts;
 if (no_density) {
   // no ligand density, just set the values
   switch (Caffe::mode()) {
   case Caffe::GPU: {
 #ifndef CPU_ONLY
     float* optgrid = net->top_vecs()[0][0]->mutable_gpu_data();
     std::vector<float> fillbuffer(lig_grid_size, carve_val);
     CUDA_CHECK(cudaMemcpy(optgrid + rec_grid_size, &fillbuffer[0], lig_grid_size * sizeof(float), cudaMemcpyHostToDevice));
 #else
     NO_GPU;
 #endif
     break;
   }
   case Caffe::CPU: {
     float* optgrid = net->top_vecs()[0][0]->mutable_cpu_data();
     std::fill(optgrid + rec_grid_size, optgrid + rec_grid_size + lig_grid_size, carve_val);
     break;
   }
   default:
     LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
   }
 }
 else {
   // some ligand density present, set value to carve_val only if it's 0
   switch (Caffe::mode()) {
   case Caffe::GPU: {
 #ifndef CPU_ONLY
     float* optgrid = net->top_vecs()[0][0]->mutable_gpu_data();
     do_constant_fill(optgrid + rec_grid_size, lig_grid_size, carve_val);
 #else
     NO_GPU;
 #endif
     break;
   }
   case Caffe::CPU: {
     float* optgrid = net->top_vecs()[0][0]->mutable_cpu_data();
     cpu_constant_fill(optgrid + rec_grid_size, lig_grid_size, carve_val);
     break;
   }
   default:
     LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
   }
 }
}

int main(int argc, char* argv[]) {
  using namespace boost::program_options;

  std::string receptor_name, grid_prefix, vsfile, solverstate, 
    outname, out_prefix, types;
  std::vector<std::string> ligand_names;

  cnn_options cnnopts;
  cnnopts.cnn_scoring = true;
  int iterations = 1000;
  int gpu;
  float base_lr = 0.1;
  bool dump_all = false;
  bool dump_last = false;
  bool exclude_receptor = false;
  bool exclude_ligand = false;
  bool ignore_ligand = false;
  bool allow_neg = false;
  bool carve = false;
  std::string sigint_effect = "stop";
  std::string sighup_effect = "snapshot";
  std::string dist_method = "l2";
  unsigned default_batch_size = 50;
  unsigned batch_size;
  unsigned nopts = 0;
  float carve_val = 1;
  std::vector<std::string> opt_names;
  float positive_threshold = 0.037526;
  float negative_threshold = 0.033237;

  positional_options_description positional; 
  options_description inputs("General Input");
  inputs.add_options()("receptor,r",
      value<std::string>(&receptor_name), "receptor to provide optimization context")(
      "ligand,l", value<std::vector<std::string> >(&ligand_names),
      "one or more ligands to provide additional optimization context")(
      "types,t", value<std::string>(&types), 
      "MolGrid .types file, formatted with one example per line, <score> <affinity> <receptor_file> <ligand_file>")(
      "grid,g", value<std::string>(&grid_prefix), 
      "prefix for grid files from which to begin optimization (instead of molecules), filenames assumed to be [prefix]_[Rec/Lig]_[channel0]_[channel1][...].dx")(
      "virtualscreen,v", value<std::string>(&vsfile), 
      "file of compounds to score according to overlap with optimized grid")(
      "positive_threshold", value<float>(&positive_threshold), 
      "optionally specify positive threshold value for threshold overlap method")(
      "negative_threshold", value<float>(&negative_threshold), 
      "optionally specify negative threshold value for threshold overlap method");

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
  output.add_options()("vs_output,o",
      value<std::string>(&outname), "name for output file of virtual screen compounds with grid overlap score, default is {receptor_file}.vsout")(
      "out_prefix,p", value<std::string>(&out_prefix),
      "prefix for storing checkpoint files for optimization in progress, default is gninadream.<PID>")(
      "dump_all", bool_switch(&dump_all)->default_value(false), 
      "dump all intermediate grids from optimization")(
      "dump_last", bool_switch(&dump_last)->default_value(false),
      "dump last grid from optimization");

  options_description options("Options");
  options.add_options()("iterations,i",
      value<int>(&iterations), "number of iterations to run, default is 1000")(
      "base_lr", value<float>(&base_lr),
      "base learning rate for density updates, default is 0.1")(
      "gpu", value<int>(&gpu)->default_value(-1), "gpu to run on")(
      "exclude_receptor", bool_switch(&exclude_receptor)->default_value(false), 
      "don't update the receptor grids during optimization")(
      "exclude_ligand",  bool_switch(&exclude_ligand)->default_value(false), 
      "don't update the ligand grids during optimization")(
      "ignore_ligand", bool_switch(&ignore_ligand)->default_value(false),
      "just use ligand to set center, don't use it to initialize grids")(
      "allow_negative", bool_switch(&allow_neg)->default_value(false),
      "allow optimization to result in negative atom density")(
      "distance", value<std::string>(&dist_method), 
      "distance function for virtual screen; options are 'l1', 'l2', 'mult', 'sum' (which is a sum of the scores produced by those two options), 'emd', and 'threshold'; default is l2")(
      "carve", bool_switch(&carve)->default_value(false),
      "do optimization from positive-initialized grids")(
      "carve_val", value<float>(&carve_val),
      "initialization value when carving, default is 1"
        );

  options_description desc;
  desc.add(inputs).add(cnn).add(output).add(options);

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
  FLAGS_alsologtostderr = 1; // for debug unless you wanna get spammed from
  // molgrid about the receptor not being set during virtual screen

  // treat rec + lig + grids as mutually exclusive with types (the former are
  // "in_mem" for MolGrid)
  // an existing solverstate is used only if none of the other options are set
  // if starting from a rec+lig or solverstate, we optimize for requested number
  // of iterations and terminate; if we have a .types file we optimize for that
  // number of iterations for each example in the batch (batching so that we
  // don't run out of memory)
  if (vm.count("receptor") <= 0 && vm.count("ligand") <= 0 && vm.count("grid") <= 0 && 
      vm.count("types") <= 0 && vm.count("solverstate") <= 0) {
    std::cerr << "Missing optimization context.\n" << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  if (exclude_receptor && exclude_ligand) {
    std::cerr << "Need to optimize at least one of ligand or receptor channels\n";
    return 1;
  }

  if (!out_prefix.size()) 
    out_prefix = "gninadream." + std::to_string(getpid()) + ".out";

  google::InitGoogleLogging(argv[0]);
  google::SetStderrLogging(2);

  // find out/set up some stuff about the CNN model
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

  // TODO: the builtin default model has separate TEST/TRAIN phase molgrid layers,
  // but in general how to deal with this? find all molgrid layers and set
  // params for all of them?
  LayerParameter* first = net_param.mutable_layer(1);
  MolGridDataParameter* mgridparam = first->mutable_molgrid_data_param();
  if (mgridparam == NULL) {
    throw usage_error("First layer of model must be MolGridData.");
  }

  mgridparam->set_random_rotation(false);
  mgridparam->set_random_translate(false);
  mgridparam->set_shuffle(false);
  mgridparam->set_balanced(false);
  mgridparam->set_stratify_receptor(false);

  if (cnnopts.cnn_model.size() == 0) {
    const char *recmap = cnn_models[cnnopts.cnn_model_name].recmap;
    const char *ligmap = cnn_models[cnnopts.cnn_model_name].ligmap;
    mgridparam->set_mem_recmap(recmap);
    mgridparam->set_mem_ligmap(ligmap);
  }

  // if we have structure file(s) or grids as input, we'll do just one example,
  // inmem=true
  if (receptor_name.size() || grid_prefix.size()) {
    mgridparam->set_inmemory(true);
    mgridparam->set_batch_size(1);
  }
  else {
    batch_size = mgridparam->batch_size();
    if (types.size())
      mgridparam->set_source(types.c_str());
  }

  net_param.set_force_backward(true);
  mgridparam->set_ignore_ligand(ignore_ligand);

  // set up solver params and then construct solver
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
  solver_param.mutable_net_param()->caffe::NetParameter::MergeFrom(net_param);
  solver_param.set_base_lr(base_lr);
  solver_param.set_max_iter(iterations);
  solver_param.set_lr_policy("fixed");
  solver_param.set_display(1);
  solver_param.set_snapshot_prefix("inputopt_");
  solver_param.set_type("InputOpt");
  if (cnnopts.cnn_weights.size())
    solver_param.add_weights(cnnopts.cnn_weights);

  caffe::shared_ptr<caffe::Solver<float> >
      solver(caffe::SolverRegistry<float>::CreateSolver(solver_param));

  solver->SetActionFunction(signal_handler.GetActionFunction());
  solver->DoThreshold(!allow_neg);
  caffe::shared_ptr<Net<float> > net = solver->net();
  // if there wasn't a weights file and we actually care about network output,
  // check that we're using one of the provided cnn models and attempt to load
  // the appropriate stored weights
  // if we don't do any optimizations of the input we don't need the network
  // for anything except loading data into grids so weights don't matter
  if (iterations != 0) {
    if (cnnopts.cnn_weights.size() == 0) {
      NetParameter wparam;

      if (cnn_models.find(cnnopts.cnn_model_name) == cnn_models.end()) {
        throw usage_error("No weights for CNN model " + cnnopts.cnn_model_name);
      }
      const unsigned char *weights = cnn_models[cnnopts.cnn_model_name].weights;
      unsigned int nbytes = cnn_models[cnnopts.cnn_model_name].num_bytes;

      google::protobuf::io::ArrayInputStream weightdata(weights,nbytes);
      google::protobuf::io::CodedInputStream strm(&weightdata);
      strm.SetTotalBytesLimit(INT_MAX, 536870912);
      bool success = wparam.ParseFromCodedStream(&strm);
      if (!success) throw usage_error("Error with default weights.");

      net->CopyTrainedLayersFrom(wparam);
    } else {
      net->CopyTrainedLayersFrom(cnnopts.cnn_weights);
    }
  }

  const vector<caffe::shared_ptr<Layer<float> > >& layers = net->layers();
  mgridT* mgrid = dynamic_cast<MolGridDataLayer<float>*>(layers[0].get());
  if (mgrid == NULL) {
    throw usage_error("First layer of model must be MolGridDataLayer.");
  }
  
  if (exclude_ligand) {
    solver->SetNligTypes(mgrid->getNumLigandTypes());
    solver->SetNpoints(mgrid->getNumGridPoints());
  }
  if (exclude_receptor) {
    solver->SetNrecTypes(mgrid->getNumReceptorTypes());
    solver->SetNpoints(mgrid->getNumGridPoints());
  }

  // figure out what the initial state should be and run optimization
  if (receptor_name.size()) {
    // can start from receptor and 0+ ligands
    // expected behavior is starting from rec with a ligand available for
    // setting center at minimum; without ligand things _should_ work (verify
    // there aren't errors about not having lig atoms) and expect that grid
    // center will be set to origin in that case. possibly that should
    // be changed (in the past when there wasn't a ligand I would set 
    // use_rec_center, so that you could provide binding site residues without
    // a reference ligand and get a reasonable result)
    tee log(true);
    FlexInfo finfo(log);
    MolGetter mols(receptor_name, std::string(), finfo, true, true, log);
    if (!ligand_names.size()) {
      ligand_names.push_back("");
    }
    boost::filesystem::path rec(receptor_name);
    std::vector<caffe::shared_ptr<std::ostream> > out;
    if (vsfile.size()) {
      if (outname.size())
        out.push_back(boost::make_shared<std::ofstream>((outname + ".vsout").c_str()));
      else
        out.push_back(boost::make_shared<std::ofstream>((rec.stem().string() + ".vsout").c_str()));
    }
    bool compute_cost = true;
    for (unsigned l = 0, nl = ligand_names.size(); l < nl; l++) {
      const std::string ligand_name = ligand_names[l];
      boost::filesystem::path lig(ligand_name);
      mols.setInputFile(ligand_name);
      for (;;)  {
        model m;

        if (!mols.readMoleculeIntoModel(m)) {
          break;
        }
        std::vector<smt> ligand_smtypes;
        std::vector<float3> ligand_coords;
        setLigand(m, ligand_coords, ligand_smtypes);
        mgrid->setLigand(ligand_coords, ligand_smtypes);
        std::vector<smt> rec_smtypes;
        std::vector<float3> rec_coords;
        setReceptor(m, rec_coords, rec_smtypes);
        mgrid->setReceptor(rec_coords, rec_smtypes);
        // with a types file you can target arbitrary pose and affinity values,
        // here we assume you just want a really good active
        mgrid->setLabels(1, 10); 
        net->ForwardFromTo(0,0);
        if (carve) {
          bool no_density = ligand_name.empty() || ignore_ligand;
          carve_init(*mgrid, no_density, net, carve_val);
        }
        std::string prefix = rec.stem().string() + "_" + lig.stem().string();
        float* data = net->top_vecs()[0][0]->mutable_cpu_data();
        unsigned dim = mgrid->getDim();
        libmolgrid::Grid<float,4> grid(data, mgrid->getNumChannels(), dim, dim, dim);
        for (size_t i=0; i<iterations; ++i) {
          solver->Step(1);
          if (dump_all || (i==iterations-1 && dump_last))
            write_dx(prefix + "_iteration" + std::to_string(i) + "_", mgrid, grid);
        }
        if (ligand_names[0] == "")
          ligand_names[0] = "none";
        if (vsfile.size())
          do_exact_vs(*net_param.mutable_layer(1), *net, vsfile, ligand_names, out, gpu+1, 
              dist_method, positive_threshold, negative_threshold, compute_cost);
        compute_cost = false;
      }
    }
  }
  else if (types.size()) {
    //molgrid will load in batches from types
    unsigned n_passes = 1;
    std::ifstream infile(types.c_str());
    if (!(bool)infile) {
      std::cerr << "Could not open " << types;
      std::exit(1);
    } 
    std::string line;
    while (getline(infile, line))
    {
      std::vector<std::string> contents;
      boost::split(contents, line, boost::is_any_of(" "));
      if (mgridparam->max_group_size()-1 || mgridparam->has_rmsd()) {
        std::cerr << "Groups and RMSD not permitted for dream optimization.\n";
        std::exit(1);
      }
      size_t offset = mgridparam->has_affinity() + 1;
      if (contents.size() < offset + 2) {
        std::cerr << ".types file input should be organized as LABEL [AFFINITY] RECFILE LIGFILE with one example per line.\n";
        std::exit(1);
      }
      boost::filesystem::path rec(contents[offset]);
      boost::filesystem::path lig(contents[offset+1]);
      opt_names.push_back(rec.stem().string() + "_" + lig.stem().string());
    }
    if (!nopts) throw usage_error("No examples in types file");
    n_passes = std::ceil((float)(nopts) / batch_size);
    for (unsigned i=0; i<n_passes; ++i) {
      net->ForwardFromTo(0,0);
      solver->ResetIter();
      if (carve) {
        // could figure out if ligand was 'none' too but that's more work 
        bool no_density = ignore_ligand;
        carve_init(*mgrid, no_density, net, carve_val);
      }
      float* data = net->top_vecs()[0][0]->mutable_cpu_data();
      unsigned dim = mgrid->getDim();
      libmolgrid::Grid<float,4> grid(data, mgrid->getNumChannels(), dim, dim, dim);
      for (size_t j=0; j<iterations; ++j) {
        solver->Step(1);
        if (dump_all || (i==iterations-1 && dump_last)) {
          for (size_t k=0; k<opt_names.size(); ++k)
            write_dx(opt_names[k] + "_iteration" + std::to_string(i) + "_", mgrid, grid);
        }
      }
      if (vsfile.size()) {
        std::vector<caffe::shared_ptr<std::ostream> > out;
        for (auto& name : opt_names) {
          out.push_back(boost::make_shared<std::ofstream>((name + ".vsout").c_str()));
        }
        do_exact_vs(*net_param.mutable_layer(1), *net, vsfile, ligand_names, out, gpu+1, 
            dist_method, positive_threshold, negative_threshold);
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
    float* inputblob = nullptr;
    double resolution = mgrid->getResolution();
    unsigned npts_onech = mgrid->getNumGridPoints();
    unsigned nchannels = rectypes.size() + ligtypes.size();
    unsigned npts_allchs = npts_onech * nchannels;
    vec gridcenter(NAN,NAN,NAN);
    if (gpu > -1)  {
      inputblob = net->top_vecs()[0][0]->mutable_gpu_data();
      CUDA_CHECK_GNINA(cudaMemset(inputblob, 0, npts_allchs* sizeof(float)));
    }
    else {
      inputblob = net->top_vecs()[0][0]->mutable_cpu_data();
      memset(inputblob, 0, npts_allchs * sizeof(float));
    }

    for (auto& fname : filenames) {
      std::string base = std::regex_replace(fname, std::regex(".dx"), "");
      std::string channels = std::regex_replace(base, std::regex(grid_prefix), "");
      std::vector<std::string> channelvec;
      boost::split(channelvec, channels, boost::is_any_of("_"));
      // the first entry in channelvec is either Rec or Lig
      std::unordered_set<std::string> channelmap(&channelvec[1], &channelvec[channelvec.size()]);
      for (auto& thisch : channelmap) {
        // does this correspond to channels in the net's map, if so what's the index
        // and use that to copy it
        unsigned idx = 0;
        if (channelvec[0] == "Rec") {
          for (auto& chanset : recset) {
            if (channelmap == chanset) {
              // cool it's a match read it in
              std::ifstream gridfile(fname.c_str());
		          if(!gridfile)
		          {
                std::cerr << "Could not open example grid file " << fname << "\n";
		          	abort();
		          }
		          vec this_gridcenter;
		          double gridres = 0;
              // if we're running on the CPU we can read directly into the
              // right place, otherwise we read into a temporary host buffer
              // and memcpy to the device
              float* g = nullptr;
              std::vector<float> hostbuf;
              if (gpu > -1) {
                hostbuf.resize(npts_onech);
                memset(&hostbuf[0], 0, npts_onech * sizeof(float));
                g = &hostbuf[0];
              }
              else 
                g = inputblob + idx * npts_onech;
		          if(!readDXGrid(gridfile, this_gridcenter, gridres, g, npts_onech, fname))
		          {
                std::cerr << "I couldn't understand the provided dx file " << fname << "\n";
                std::exit(1);
		          }
              if (gridres - resolution > 0.001) {
                std::cerr << "grid resolution in file " << fname << " does not match grid resolution specified for trained net.\n" << std::endl;
                std::exit(1);
              }
              if (gpu > -1) 
                CUDA_CHECK_GNINA(cudaMemcpy(inputblob + idx * npts_onech, g,
                      npts_onech * sizeof(float), cudaMemcpyHostToDevice));
              if (std::isnan(gridcenter[0])) {
                gridcenter = this_gridcenter;
              }
              else {
                assert(gridcenter == this_gridcenter);
              }
              break;
            }
            ++idx;
          }
        }
        else if (channelvec[0] == "Lig") {
          idx = rectypes.size();
          for (auto& chanset : ligset) {
            if (channelmap == chanset) {
              // cool it's a match read it in
              std::ifstream gridfile(fname.c_str());
		          if(!gridfile)
		          {
                std::cerr << "Could not open example grid file " << fname << "\n";
		          	abort();
		          }
		          vec this_gridcenter;
		          double gridres = 0;
              // if we're running on the CPU we can read directly into the
              // right place, otherwise we read into a temporary host buffer
              // and memcpy to the device
              float* g = nullptr;
              std::vector<float> hostbuf;
              if (gpu > -1) {
                hostbuf.resize(npts_onech);
                g = &hostbuf[0];
              }
              else 
                g = inputblob + idx * npts_onech;
		          if(!readDXGrid(gridfile, this_gridcenter, gridres, g, npts_onech, fname))
		          {
                std::cerr << "I couldn't understand the provided dx file " << fname << "\n";
                std::exit(1);
		          }
              if (gridres - resolution > 0.001) {
                std::cerr << "grid resolution in file " << fname << " does not match grid resolution specified for trained net.\n" << std::endl;
                std::exit(1);
              }
              if (gpu > -1) 
                CUDA_CHECK_GNINA(cudaMemcpy(inputblob + idx * npts_onech, g,
                      npts_onech * sizeof(float), cudaMemcpyHostToDevice));
              if (std::isnan(gridcenter[0])) {
                gridcenter = this_gridcenter;
              }
              else {
                if (this_gridcenter != gridcenter) {
                  std::cerr << "Provided grids don't all have the same center.\n";
                  std::exit(1);
                }
              }
              break;
            }
            ++idx;
          }
        }
        else {
          std::cerr << "Grid files don't match pattern PREFIX_(Rec|Lig)_CHANNELS.\n";
          std::exit(1);
        }
      }
    }
    if (carve) {
      bool no_density = false;
      carve_init(*mgrid, no_density, net, carve_val);
    }
    float* data = net->top_vecs()[0][0]->mutable_cpu_data();
    unsigned dim = mgrid->getDim();
    libmolgrid::Grid<float,4> grid(data, mgrid->getNumChannels(), dim, dim, dim);
    for (size_t i=0; i<iterations; ++i) {
      solver->Step(1);
      if (dump_all || (i==iterations-1 && dump_last))
        write_dx(out_prefix + "_iteration" + std::to_string(i) + "_", mgrid, grid);
    }
    if (vsfile.size()) {
      // make temp gninatypes for setting the center for screen ligands
      string gt_name = "tmp_" + std::to_string(getpid()) + ".gninatypes";
	  	std::ofstream gt(gt_name.c_str());
	  	if (!gt)
	  	{
        std::cerr << "Error opening output file " << gt_name << "\n";
	  		exit(1);
	  	}
      atom_info atom;
      atom.x = gridcenter[0];
      atom.y = gridcenter[1];
      atom.z = gridcenter[2];
      atom.type = 2;
	  	gt.write((char*)&atom, sizeof(atom));
	  	gt.close();
      if (!ligand_names.size()) {
        ligand_names.push_back(gt_name);
      }
      std::vector<caffe::shared_ptr<std::ostream> > out;
      out.push_back(boost::make_shared<std::ofstream>((out_prefix + ".vsout").c_str()));
      do_exact_vs(*net_param.mutable_layer(1), *net, vsfile, ligand_names, out, gpu+1, 
          dist_method, positive_threshold, negative_threshold);
      std::remove(gt_name.c_str());
    }
  }
  else if (solverstate.size()) {
    const char* ss_cstr = solverstate.c_str();
    //restart from optimization in progress
    solver->Restore(ss_cstr);
    int remaining_iters = iterations - solver->iter();
    float* data = net->top_vecs()[0][0]->mutable_cpu_data();
    unsigned dim = mgrid->getDim();
    libmolgrid::Grid<float,4> grid(data, mgrid->getNumChannels(), dim, dim, dim);
    for (size_t i=0; i<remaining_iters; ++i) {
      solver->Step(1);
      if (dump_all || (i==iterations-1 && dump_last))
        write_dx(out_prefix + "_iteration" + std::to_string(i) + "_", mgrid, grid);
    }
    if (vsfile.size()) {
      std::vector<caffe::shared_ptr<std::ostream> > out;
      out.push_back(boost::make_shared<std::ofstream>((out_prefix + ".vsout").c_str()));
      do_exact_vs(*net_param.mutable_layer(1), *net, vsfile, ligand_names, out, gpu+1, 
          dist_method, positive_threshold, negative_threshold);
    }
  }
  else {
    std::cerr << "No valid initial input for optimization provided.\n";
    std::exit(1);
  }
}
