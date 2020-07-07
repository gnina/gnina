#include <iostream>
#include "cnn_visualization.hpp"
#include "../lib/cnn_scorer.h"
#include "../lib/tee.h"
#include "../lib/flexinfo.h"
#include "../lib/box.h"
#include <vector>
#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {
  bool zero_values;
  vis_options visopts;
  cnn_options cnnopts;

  cnnopts.cnn_rotations = 0; //any reason to make this an option?
  cnnopts.cnn_scoring = CNNall;

  using namespace boost::program_options;

  std::string vis_method;
  std::string target;

  options_description inputs("Input");
  inputs.add_options()("receptor, r",
      value<std::string>(&visopts.receptor_name), "receptor for coloring")(
      "ligand, l", value<std::string>(&visopts.ligand_name),
      "ligand for coloring");

  options_description cnn("CNN Input");
  std::string model;
  std::string weights;
  cnn.add_options()("cnn_model", value<std::string>(&model),
      "CNN model file (*.model)")
      ("cnn_weights", value<std::string>(&weights),
      "CNN weights file (*.caffemodel)")("ignore_layer",
      value<std::string>(&visopts.layer_to_ignore)->default_value(""),
      "zero values in layer with provided name, in the case of split output");

  options_description output("Output");
  output.add_options()("skip_ligand_output",
      bool_switch(&visopts.skip_ligand_output), "skip ligand visualization")(
      "skip_receptor_output", bool_switch(&visopts.skip_receptor_output),
      "skip receptor visualization");

  options_description options("Misc options");
  options.add_options()("box_size",
      value<float>(&visopts.box_size)->default_value(23.5),
      "diameter of bounding box for receptor coloring, in angstroms (default 23.5)")
      ("frags_only", bool_switch(&visopts.frags_only)->default_value(false),
      "only run fragment removal on ligand")("atoms_only",
      bool_switch(&visopts.atoms_only)->default_value(false),
      "only run individual removal on ligand")("verbose",
      bool_switch(&visopts.verbose)->default_value(false),
      "print full output, including removed atom lists")("gpu",
      value<int>(&visopts.gpu)->default_value(-1),
      "gpu id for accelerated scoring")("vis_method",
      value<std::string>(&vis_method)->default_value("masking"),
      "visualization method (lrp, masking, gradient, or all)")("target",
      value<std::string>(&visopts.target)->default_value("pose"),
      "scoring method for masking (pose or aff)")("outputdx",
      bool_switch(&visopts.outputdx)->default_value(false),
      "output DX grid files")("zero_values",
      bool_switch(&visopts.zero_values)->default_value(false),
      "only propagate values from dead nodes")
      ("score_scale", value<double>(&visopts.score_scale)->default_value(10.0),
            "Amount to scale score output by (default 10.0)");

  options_description debug("Debug");
  debug.add_options()("output_files",
      bool_switch(&visopts.output_files)->default_value(false),
      "write every modified pdbqt file")("additivity",
      value<std::string>(&visopts.additivity),
      "print additivity data for ligand")("skip_bound_check",
      bool_switch(&visopts.skip_bound_check)->default_value(false),
      "score all residues, regardless of proximity to ligand");
  options_description desc;
  desc.add(inputs).add(cnn).add(output).add(options).add(debug);

  positional_options_description positional; // remains empty?
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

  //placeholders for center to instantiate
  float center_x = 0, center_y = 0, center_z = 0;
  vec center(center_x, center_y, center_z);

  if (vm.count("receptor") <= 0) {
    std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  if (vm.count("ligand") <= 0) {
    std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  if (visopts.frags_only && visopts.atoms_only) {
    std::cerr << "Cannot use 'frags_only' and 'atoms_only' together.\n"
        << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  cnnopts.cnn_models.push_back(model);
  cnnopts.cnn_weights.push_back(weights);

  google::InitGoogleLogging(argv[0]);
  google::SetStderrLogging(2);

  cnn_visualization vis = cnn_visualization(visopts, cnnopts, center);

  if ("masking" == vis_method) {
    vis.masking();
  } else
    if ("lrp" == vis_method) {
      vis.lrp();
    } else
      if ("gradient" == vis_method) {
        vis.gradient_vis();
      } else
        if ("all" == vis_method) {
          vis.gradient_vis();
          vis.lrp();
          vis.masking();
        } else {
          std::cerr
              << "Specified vis_method not known. Use \"masking\", \"lrp\", or \"gradient\"\n";
          return 1;
        }
}
