//============================================================================
// Name        : fromsmina.cpp
// Author      : David Koes
// Created     : Jul 14, 2015
// Copyright   : 2015 University of Pittsburgh
// Description : ?Convert smina format to sdf, for debugging.
//============================================================================

#include <iostream>
#include <fstream>
#include "CommandLine2/CommandLine.h"
#include "GninaConverter.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace std;
using namespace OpenBabel;

cl::opt<string> infile("in", cl::desc("input file"), cl::Required,
    cl::Positional);

int main(int argc, char *argv[]) {
  cl::ParseCommandLineOptions(argc, argv);

  ifstream ifile(infile.c_str());

  model initm;
  size_t position = 0;
  while (ifile) {
    unsigned sz;
    ifile.read((char*) &sz, sizeof(sz));
    char buffer[sz + 1];
    ifile.read(buffer, sz);

    if (!ifile) break;

    stringstream data;
    for (unsigned i = 0; i < sz; i++) {
      data.put(buffer[i]);
    }

    boost::iostreams::filtering_stream<boost::iostreams::input> io_strm;
    io_strm.push(boost::iostreams::gzip_decompressor());

    io_strm.push(data);

    boost::archive::binary_iarchive serialin(io_strm,
        boost::archive::no_header | boost::archive::no_tracking);

    unsigned numtors;
    parsing_struct p;
    context c;

    serialin >> numtors;
    serialin >> p;
    serialin >> c;

    model m = initm;

    non_rigid_parsed nr;
    pdbqt_initializer tmp;

    postprocess_ligand(nr, p, c, numtors);
    tmp.initialize_from_nrp(nr, c, true);
    tmp.initialize(nr.mobility_matrix());
    m.set_name(c.sdftext.name);

    m.append(tmp.m);

    stringstream str;
    m.write_sdf(str);
    cout << str.str() << "$$$$\n";
  }
}
