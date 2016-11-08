#include <iostream>
#include "cnn_visualization.hpp"
#include "../lib/cnn_scorer.h"
#include "../lib/tee.h"
#include "../lib/flexinfo.h"
#include "../lib/box.h"
#include <vector>
#include <boost/program_options.hpp>

int main(int argc, char* argv[])
{
  vis_options visopts;
  cnn_options cnnopts;

  cnnopts.cnn_rotations = 24; //any reason to make this an option?
  cnnopts.cnn_scoring = true;

  using namespace boost::program_options;

  options_description inputs("Input");
  inputs.add_options()
    ("receptor, r", value<std::string>(&visopts.receptor_name),
                  "receptor for coloring")
    ("ligand, l", value<std::string>(&visopts.ligand_name),
                  "ligand for coloring");
    
  options_description cnn("CNN Input");
  cnn.add_options()
    ("cnn_model", value<std::string>(&cnnopts.cnn_model),
                  "CNN model file (*.model)")
    ("cnn_weights", value<std::string>(&cnnopts.cnn_weights),
                  "CNN weights file (*.caffemodel)");

  options_description outputs("Output");
  outputs.add_options()
    ("receptor_output", value<std::string>(&visopts.receptor_output),
                  "output file for colored receptor (omit to skip receptor coloring)")
    ("ligand_output", value<std::string>(&visopts.ligand_output),
                  "output file for colored ligand (omit to skip ligand coloring)");

  options_description options("Misc options");
  options.add_options()
    ("box_size", value<float>(&visopts.box_size)->default_value(23.5),
                  "diameter of bounding box for receptor coloring, in angstroms (default 23.5)")
    ("frags_only", bool_switch(&visopts.frags_only)->default_value(false),
                  "only run fragment removal on ligand")
    ("atoms_only", bool_switch(&visopts.atoms_only)->default_value(false),
                  "only run individual removal on ligand")
    ("verbose", bool_switch(&visopts.verbose)->default_value(false),
                  "print full output, including removed atom lists");

  options_description debug("Debug");
  debug.add_options()
    ("output_files", bool_switch(&visopts.output_files)->default_value(false),
                    "write every modified pdbqt file")
    ("additivity", bool_switch(&visopts.additivity)->default_value(false),
                    "print additivity data for ligand");

  options_description desc;
  desc.add(inputs).add(cnn).add(outputs).add(options).add(debug);

  positional_options_description positional; // remains empty?
  variables_map vm;
  try
  {
    store(

          command_line_parser(argc,argv).options(desc)
                .style(
                  command_line_style::default_style
                          ^ command_line_style::allow_guessing)
                .positional(positional).run(), vm);
    notify(vm);
  }
  catch (boost::program_options::error& e)
  {
    std::cerr << "Command line parse error: " << e.what() << '\n'
              << "\nCorrect usage:\n" << desc << '\n';
    return 1;
  }

  if(vm.count("receptor") <= 0)
  {
    std::cerr << "Missing receptor.\n" << "\nCorrect usage:\n"
            << desc << '\n';
    return 1;
  }
  
  if(vm.count("ligand") <= 0)
  {
    std::cerr << "Missing ligand.\n" << "\nCorrect usage:\n"
            << desc << '\n';
    return 1;
  }

  if(vm.count("cnn_model") <= 0)
  {
    std::cerr << "Missing cnn_model.\n" << "\nCorrect usage:\n"
            << desc << '\n';
    return 1;
  }

  if(vm.count("cnn_weights") <= 0)
  {
    std::cerr << "Missing cnn_weights.\n" << "\nCorrect usage:\n"
            << desc << '\n';
    return 1;
  }

  if(vm.count("receptor_output") <= 0 && vm.count("ligand_output") <= 0)
  {
    std::cerr << "At least one of 'receptor_output' and 'ligand_output' required.\n" << "\nCorrect usage:\n"
            << desc << '\n';
    return 1;
  }

  if(visopts.frags_only && visopts.atoms_only)
  {
    std::cerr << "Cannot use 'frags_only' and 'atoms_only' together.\n" << "\nCorrect usage:\n"
            << desc << '\n';
    return 1;
  }
  
  google::InitGoogleLogging(argv[0]);
  google::SetStderrLogging(2);

  //placeholders for center to instantiate
  float center_x = 0, center_y = 0, center_z = 0;
  vec center(center_x,center_y,center_z);

  cnn_visualization vis = cnn_visualization(visopts, cnnopts, center);
  vis.color();
}
