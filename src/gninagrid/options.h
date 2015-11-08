#ifndef GRIDDER_OPTIONS_H
#define GRIDDER_OPTIONS_H

#include <string>

using namespace std;

/* options that can be set on the command line */
struct cmdoptions
{
	string receptorfile;
	string ligandfile;
	string outname;
	string recmap;
	string ligmap;
	double dim;
	double res;
	double x;
	double y;
	double z;
	int verbosity;
	bool help;
	bool version;
	bool outmap;
	bool binary;

	cmdoptions() :
			dim(24), res(0.5), x(HUGE_VAL), y(HUGE_VAL), z(HUGE_VAL),
			verbosity(1), help(false), version(false), outmap(false), binary(false)
	{
	}
};

#endif
