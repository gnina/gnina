#include <boost/program_options.hpp>
#include "../lib/cnn_scorer.h"

int main(int argc, char* argv[]) {
  using namespace boost::program_options;

  std::string receptor_name, ligand_name, grid_prefix, vsfile, solverstate, 
    outname, out_prefix, types, box_ligand;

  cnn_options cnnopts;
  cnnopts.cnn_scoring = True;
  int iterations;
  int gpu;
  float base_lr;
  bool dump_all = false;
  bool exclude_receptor = false;
  bool exclude_ligand = false;
  std::string sigint_effect = "stop";
  std::string sighup_effect = "snapshot";

  options_description inputs("General Input");
  inputs.add_options()("receptor, r",
      value<std::string>(&receptor_name), "receptor to provide optimization context")(
      "ligand, l", value<std::string>(&ligand_name),
      "ligand to provide optimization context")(
      "types, t", value<std::string>(&types), 
      "MolGrid .types file, formatted with one example per line, <score> <affinity> <receptor_file> <ligand_file>")(
      "grid, g", value<std::string>(&grid_prefix), 
      "prefix for grid files from which to begin optimization (instead of molecules), filenames assumed to be [prefix]_[channel].dx")(
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
      "dump all intermediate grids from optimization");

  options_description options("Options");
  options.add_options()("iterations, i",
      value<int>(&iterations), "number of iterations to run")(
      "base_lr", value<float>(&base_lr),
      "base learning rate for density updates")(
      "auto_center", value<std::string>(&box_ligand), 
      "provide reference ligand just to set grid center")(
      "gpu, g", value<int>(&gpu)->default_value(-1), "gpu to run on")(
      "exclude_receptor, er", bool_switch(&exclude_receptor)->default_value(false), 
      "don't update the receptor grids")(
      "exclude_ligand, el",  bool_switch(&exclude_ligand)->default_value(false), 
      "don't update the ligand grids");

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

  //right now no "user grids", we assume the grids are also rec/lig and load
  //any/all of rec + lig + grids in at once
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
  mgridnet_param = first->mutable_molgrid_data_net_param();
  if (mgridnet_param == NULL) {
    throw usage_error("First layer of model must be MolGridData.");
  }

  if (cnnopts.cnn_model.size() == 0) {
    const char *recmap = cnn_models[cnnopts.cnn_model_name].recmap;
    const char *ligmap = cnn_models[cnnopts.cnn_model_name].ligmap;
    mgridnet_param->set_mem_recmap(recmap);
    mgridnet_param->set_mem_ligmap(ligmap);
  }

  net_param.set_force_backward(true);
  unsigned batch_size = net_param.batch_size();

  //set up solver net_params and then construct solver
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


  if (!receptor.empty()) {
    //set inmem
    // mgridparam->set_inmemory(true);
    // batch_size = 1;
    std::cerr << "Starting inmem not supported yet\n";
    exit(-1);
  }
  else if (!types.empty()) {
    //molgrid will load in batches from types
    unsigned n_examples;
    unsigned n_iters;
    std::ifstream infile(types.c_str());
    CHECK((bool)infile) << "Could not open " << source;
    while (getline(infile, line))
    {
      ++n_examples;
    }
  }
  else if (!grid_prefix.empty()) {
    std::cerr << "Importing raw grids not supported yet\n";
    exit(-1);
  }
  else if (!solverstate.empty()) {
    //restart from optimization in progress
    std::cerr << "Restarting from previous solverstate not supported yet\n";
    exit(-1);
    // solver.Solve(solverstate.c_str());
  }
  if (!vsfile.empty()) {
    //use the resulting grids to perform a virtual screen against these
    //compounds; if there were multiple inputs optimized, return a separate
    //preds file for each of them 
    //right now we only support an L2 loss on the grid; generate grids for the
    //virtual screen compounds centered on whatever the grid center for the
    //inputopt grids was
  }
}
