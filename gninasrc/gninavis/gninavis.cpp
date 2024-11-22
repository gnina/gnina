#include "../lib/box.h"
#include "../lib/dl_scorer.h"
#include "../lib/flexinfo.h"
#include "../lib/tee.h"
#include "../lib/torch_models.h"
#include "cnn_visualization.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  bool zero_values;
  vis_options visopts;
  cnn_options cnnopts;

  cnnopts.cnn_rotations = 0; // any reason to make this an option?
  cnnopts.cnn_scoring = CNNall;

  using namespace boost::program_options;

  std::string vis_method;
  std::string target;

  options_description inputs("Input");
  inputs.add_options()("receptor, r", value<std::string>(&visopts.receptor_name), "receptor for coloring")(
      "ligand, l", value<std::string>(&visopts.ligand_name), "ligand for coloring");

  options_description cnn("CNN Input");
  std::string model;
  std::string weights;
  cnn.add_options() //
      ("cnn", value<std::vector<std::string>>(&cnnopts.cnn_model_names)->multitoken(),
       ("built-in model to use, specify PREFIX_ensemble to evaluate an ensemble of models starting with PREFIX: " +
        builtin_torch_models())
           .c_str()) //
      ("cnn_model", value<std::vector<std::string>>(&cnnopts.cnn_models)->multitoken(),
       "cmodel file; if not specified a default model will be used");
  ;

  options_description output("Output");
  output.add_options()("skip_ligand_output", bool_switch(&visopts.skip_ligand_output), "skip ligand visualization") //
      ("skip_receptor_output", bool_switch(&visopts.skip_receptor_output), "skip receptor visualization");

  options_description options("Misc options");
  options.add_options() //
      ("box_size", value<float>(&visopts.box_size)->default_value(23.5),
       "diameter of bounding box for receptor coloring, in angstroms (default 23.5)")                               //
      ("frags_only", bool_switch(&visopts.frags_only)->default_value(false), "only run fragment removal on ligand") //
      ("atoms_only", bool_switch(&visopts.atoms_only)->default_value(false),
       "only run individual removal on ligand") //
      ("verbose", bool_switch(&visopts.verbose)->default_value(false),
       "print full output, including removed atom lists") //
      ("target", value<std::string>(&visopts.target)->default_value("pose"),
       "scoring method for masking (pose or aff)")                                               //
      ("score_scale", value<double>(&visopts.score_scale)->default_value(10.0),
       "Amount to scale score output by (default 10.0)");

  options_description debug("Debug");
  debug.add_options()("output_files", bool_switch(&visopts.output_files)->default_value(false),
                      "write every modified pdbqt file")("additivity", value<std::string>(&visopts.additivity),
                                                         "print additivity data for ligand")(
      "skip_bound_check", bool_switch(&visopts.skip_bound_check)->default_value(false),
      "score all residues, regardless of proximity to ligand");
  options_description desc;
  desc.add(inputs).add(cnn).add(output).add(options).add(debug);

  positional_options_description positional; // remains empty?
  variables_map vm;
  try {
    store(

        command_line_parser(argc, argv)
            .options(desc)
            .style(command_line_style::default_style ^ command_line_style::allow_guessing)
            .positional(positional)
            .run(),
        vm);
    notify(vm);
  } catch (boost::program_options::error &e) {
    std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  // placeholders for center to instantiate
  float center_x = 0, center_y = 0, center_z = 0;
  vec center(center_x, center_y, center_z);

  if (vm.count("receptor") <= 0) {
    std::cerr << "Missing receptor.\n"
              << "\nCorrect usage:\n"
              << desc << '\n';
    return 1;
  }

  if (vm.count("ligand") <= 0) {
    std::cerr << "Missing ligand.\n"
              << "\nCorrect usage:\n"
              << desc << '\n';
    return 1;
  }

  if (visopts.frags_only && visopts.atoms_only) {
    std::cerr << "Cannot use 'frags_only' and 'atoms_only' together.\n"
              << "\nCorrect usage:\n"
              << desc << '\n';
    return 1;
  }

  cnn_visualization vis = cnn_visualization(visopts, cnnopts, center);
  vis.masking(); 
  //TODO: reimplement at least gradient_vis in torch

}
