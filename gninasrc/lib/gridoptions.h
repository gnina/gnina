#ifndef GRIDDER_OPTIONS_H
#define GRIDDER_OPTIONS_H

#include <string>
#include <time.h>

using namespace std;

/* options that can be set on the command line */
struct gridoptions
{
	string receptorfile;
	string ligandfile;
	string outname;
	string recmap;
	string ligmap;
	double dim;
	double res;
	fl randtranslate;
	int verbosity;
	int seed;
	bool randrotate;
	bool help;
	bool version;
	bool timeit;
	bool outmap;
	bool outdx; //unfortunately mutually exclusive with outmap
	bool binary;
	bool gpu;
	gridoptions() :
		//a default dimension of 23.5 yields 48x48x48 gridpoints
			dim(23.5), res(0.5), verbosity(1), seed((int)time(NULL)),
			randrotate(false), randtranslate(0.0),
			help(false), version(false), timeit(false), outmap(false), binary(false), gpu(false)
	{
	}
};

#endif
