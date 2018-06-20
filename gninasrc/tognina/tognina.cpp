//============================================================================
// Name        : tosmina.cpp
// Author      : David Koes
// Created     : Jun 4, 2014
// Copyright   : 2014 University of Pittsburgh
// Description : Convert a ligand, read in using OpenBabel, to smina format
//============================================================================

#include <iostream>
#include <fstream>
#include "obmolopener.h" //copied smina
#include "CommandLine2/CommandLine.h"
#include <openbabel/mol.h>
#include "GninaConverter.h"

using namespace std;
using namespace OpenBabel;

cl::opt<string> infile("in", cl::desc("input file"), cl::Required,
    cl::Positional);
cl::opt<string> outfile("out", cl::desc("output file"), cl::Required,
    cl::Positional);
cl::opt<bool> textOutput("text", cl::desc("produce text output"));

int main(int argc, char *argv[]) {
  cl::ParseCommandLineOptions(argc, argv);

  OBConversion conv;

  obmol_opener opener;
  opener.openForInput(conv, infile);

  ostream *out = NULL;
  ofstream outf;
  string outname(outfile);
  if (outname != "-") {
    outf.open(outfile.c_str());
    out = &outf;
  } else //stdout
  {
    out = &cout;
  }

  OBMol mol;
  while (conv.Read(&mol)) {
    if (textOutput)
      GninaConverter::convertText(mol, *out);
    else
      GninaConverter::convertBinary(mol, *out);
  }
  return 0;
}
