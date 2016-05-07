#ifndef GRIDDER_OPTIONS_H
#define GRIDDER_OPTIONS_H

#include <string>

using namespace std;

/* options that can be set on the command line */
struct gridoptions
{
	string receptorfile;
	string ligandfile;
	string centerfile; //ligand to use to compute center
	string outname;
	string ligoutname; //first separate maps
	string recmap;
	string ligmap;
	double dim;
	double res;
	fl x;
	fl y;
	fl z;
	fl randtranslate;
	int verbosity;
	int seed;
	bool randrotate;
	bool help;
	bool version;
	bool outmap;
	bool binary;
	gridoptions() :
			dim(24), res(0.5), x(HUGE_VAL), y(HUGE_VAL), z(HUGE_VAL),
			verbosity(1), seed(0), randrotate(false), randtranslate(0.0), help(false), version(false), outmap(false), binary(false)
	{
	}
};

#endif
