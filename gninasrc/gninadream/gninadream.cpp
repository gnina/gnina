#include <boost/program_options.hpp>
#include "../lib/cnn_scorer.h"

int main(int argc, char* argv[]) {
  using namespace boost::program_options;

  std::string receptor_name, ligand_name, grid_prefix, vsfile, solverstate, 
    outname, out_prefix, types, box_ligand;

  int iterations;
  int gpu;
  bool dump_all = false;
  bool exclude_receptor = false;
  bool exclude_ligand = false;

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
      "auto_center", value<std::string>(&box_ligand), 
      "provide reference ligand just to set grid center")(
      "gpu, g", value<int>(&gpu), "gpu to run on")(
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
  //cnn_weights are used only if none of the other options are set
  if (vm.count("receptor") <= 0 && vm.count("ligand") <= 0 && vm.count("grid") <= 0 && 
      vm.count("types") <= 0 && vm.count("cnn_weights") <= 0) {
    std::cerr << "Missing optimization context.\n" << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  google::InitGoogleLogging(argv[0]);
  google::SetStderrLogging(2);

}
