/*
 * gninagrid.cpp
 *
 *  Created on: Nov 4, 2015
 *      Author: dkoes
 *
 * Output a voxelation of a provided receptor and ligand.
 * For every (heavy) atom type and grid point compute an occupancy value.
 */

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "atom_type.h"
#include "box.h"
#include "molgetter.h"

using namespace std;
using namespace boost;

/* options that can be set on the command line */
struct cmdoptions
{
	string receptorfile;
	string ligandfile;
	string outname;
	double dim;
	double res;
	double x;
	double y;
	double z;
	int verbosity;
	bool help;
	bool version;

	cmdoptions() :
			dim(16), res(0.5), x(HUGE_VAL), y(HUGE_VAL), z(HUGE_VAL),
			verbosity(1), help(false), version(false)
	{
	}
};

//parse commandline options using boost::program_options and put the values in opts
//return true if successfull and ready to compute
//will exit on error
static bool parse_options(int argc, char *argv[], cmdoptions& o)
{
	using namespace boost::program_options;
	positional_options_description positional; // remains empty

	options_description inputs("Input");
	inputs.add_options()
	("receptor,r", value<std::string>(&o.receptorfile),
			"receptor file")
	("ligand,l", value<std::string>(&o.ligandfile), "ligand(s)");

	options_description outputs("Output");
	outputs.add_options()
	("out,o", value<std::string>(&o.outname),
			"output file name, format taken from file extension");

	options_description options("Options");
	options.add_options()
	("dimension,d", value<double>(&o.dim), "Cubic grid dimension (Angstroms)")
	("resolution,r", value<double>(&o.res), "Cubic grid resolution (Angstroms)")
	("center_x", value<double>(&o.x), "X coordinate of the center, if unspecified use first ligand")
	("center_y", value<double>(&o.y), "Y coordinate of the center, if unspecified use first ligand")
	("center_z", value<double>(&o.z), "Z coordinate of the center, if unspecified use first ligand");

	options_description info("Information (optional)");
	info.add_options()
	("help", bool_switch(&o.help), "display usage summary")
	("version", bool_switch(&o.version), "display program version")
	("verbosity", value<int>(&o.verbosity)->default_value(1),
			"Adjust the verbosity of the output, default: 1");
	options_description desc;
	desc.add(inputs).add(options).add(outputs).add(info);
	variables_map vm;
	try
	{
		store(
			command_line_parser(argc, argv).options(desc)
						.style(
						command_line_style::default_style
								^ command_line_style::allow_guessing)
						.positional(positional).run(), vm);
		notify(vm);
	} catch (boost::program_options::error& e)
	{
		std::cerr << "Command line parse error: " << e.what() << '\n'
				<< "\nCorrect usage:\n" << desc << '\n';
		exit(-1);
	}

	//process informational
	if (o.help)
	{
		cout << desc << '\n';
		return false;
	}
	if (o.version)
	{
		cout << "gnina "  __DATE__ << '\n';
		return false;
	}

	return true;
}

/* Maintains atom grid information.  Stores a model of the receptor/ligand with
 * MolGetter, but also numerical grids for every protein/ligand atom type.
 */
class NNGridder
{
	MolGetter mols; //this stores the models
	grid_dims dims;
public:
	NNGridder(const cmdoptions& opt) {
		tee log(true);
		FlexInfo finfo(log); //dummy
		mols.create_init_model(opt.receptorfile, "", finfo, log);

		int numpts = round(opt.dim/opt.res);
		double half = opt.dim/2.0;
		dims[0].begin = opt.x - half;
		dims[0].end = opt.x + half;
		dims[0].n = numpts;

		dims[1].begin = opt.y - half;
		dims[1].end = opt.y + half;
		dims[1].n = numpts;

		dims[2].begin = opt.z - half;
		dims[2].end = opt.z + half;
		dims[2].n = numpts;
	}
};

int main(int argc, char *argv[])
{
	//setup commandline options
	cmdoptions opt;
	if(!parse_options(argc, argv, opt))
		exit(0);

	//figure out grid center
	if(!isfinite(opt.x + opt.y + opt.z))
	{
		fl dummy; //we wil set the size
		setup_autobox(opt.ligandfile, 0, opt.x,opt.y, opt.z, dummy, dummy, dummy);
	}

	//setup receptor grid
	NNGridder gridder(opt);

	//for each ligand..

		//compute ligand grid

		//and output

}
