/*
 * obmolopener.cpp
 *
 *  Created on: Jan 10, 2013
 *      Author: dkoes
 */

#include "obmolopener.h"
#include "file.h"
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/timer/timer.hpp>
#include <boost/iostreams/device/null.hpp>

using namespace OpenBabel;
using namespace boost;
using namespace boost::iostreams;

void obmol_opener::clear() {
  for (unsigned i = 0, n = streams.size(); i < n; i++) {
    delete streams[i];
  }
  streams.clear();
}

obmol_opener::~obmol_opener() {
  clear();
}

void obmol_opener::openForInput(OBConversion& conv, const std::string& name) {
  OBFormat *format = conv.FormatFromExt(name);

  if (!format || !conv.SetInFormat(format)) {
    throw file_error(path(name), true);
  }

  //dkoes - annoyingly, although openbabel is smart enough to ignore
  //the .gz at the end of a file when determining the file format, it
  //does not actually open the file as a gzip stream
  std::ifstream *uncompressed_inmol = new std::ifstream(name.c_str());
  streams.push_back(uncompressed_inmol);
  filtering_stream<input> *inmol = new filtering_stream<input>();
  streams.push_back((std::istream*) inmol);

  std::string::size_type pos = name.rfind(".gz");
  if (pos != std::string::npos) {
    inmol->push(gzip_decompressor());
  }
  inmol->push(*uncompressed_inmol);

  if (!*uncompressed_inmol || !*inmol) {
    throw file_error(path(name), true);
  }
  conv.SetInStream((std::istream*) inmol);

}

void obmol_opener::obmol_opener::openForOutput(OBConversion& outconv,
    const std::string& outname) {
  OBFormat *outformat = outconv.FormatFromExt(outname);
  if (!outformat || !outconv.SetOutFormat(outformat)) {
    throw file_error(outname, false);
  }

  std::ofstream *uncompressed_outfile = new std::ofstream(outname.c_str());
  filtering_stream<output>* outfile = new filtering_stream<output>();
  streams.push_back((std::ostream*) outfile);
  streams.push_back(uncompressed_outfile); //has to be deleted after filter

  std::string::size_type pos = outname.rfind(".gz");
  if (pos != std::string::npos) {
    outfile->push(gzip_compressor());
  }
  outfile->push(*uncompressed_outfile);
  if (!*outfile || !*uncompressed_outfile) {
    throw file_error(path(outname), false);
  }
  outconv.SetOutStream((std::ostream*) outfile);

}
