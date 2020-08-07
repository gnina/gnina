/*

 Copyright (c) 2006-2010, The Scripps Research Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 Author: Dr. Oleg Trott <ot14@columbia.edu>, 
 The Olson Lab, 
 The Scripps Research Institute

 */

#ifndef VINA_FILE_H
#define VINA_FILE_H

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
#include "common.h"

struct file_error {
    path name;
    bool in;
    file_error(const path& name_, bool in_)
        : name(name_), in(in_) {
    }
};

struct ifile : public boost::filesystem::ifstream { // never use ifstream pointer to destroy ifile - no virtual destructors, possibly
    ifile(const path& name)
        : boost::filesystem::ifstream(name) {
      if (!(*this)) throw file_error(name, true);
    }
    ifile(const path& name, std::ios_base::openmode mode)
        : boost::filesystem::ifstream(name, mode) {
      if (!(*this)) throw file_error(name, true);
    }
};

struct ofile : public boost::filesystem::ofstream { // never use ofstream pointer to destroy ofile - no virtual destructors, possibly
    ofile(const path& name)
        : boost::filesystem::ofstream(name) {
      if (!(*this)) throw file_error(name, false);
    }
    ofile(const path& name, std::ios_base::openmode mode)
        : boost::filesystem::ofstream(name, mode) {
      if (!(*this)) throw file_error(name, false);
    }
};

//dkoes - wrapper for an input file that is optionally gzipped
//name ends in .gz
class izfile : public boost::iostreams::filtering_stream<boost::iostreams::input> {

    std::ifstream uncompressed_infile;
    bool iszipped;
  public:

    izfile() {
    }

    izfile(const path& name, const std::string& ext, bool is_compressed = false)
        : iszipped(is_compressed) {
      open(name, ext);
    }

    //opens file name, but only if it has the appropriate ext
    //otherwise returns false
    bool open(const path& name, const std::string& ext, bool is_compressed =
        false) {
      using namespace boost::filesystem;
      iszipped = is_compressed;
      //clean up if we are already open
      while (!empty())
        pop();
      uncompressed_infile.close();

      std::string fileext = extension(name);
      if (fileext == ".gz") {
        fileext = extension(basename(name));
        iszipped = true;
      }

      if (fileext != ext) return false; //wrong type of file

      uncompressed_infile.open(name.c_str());

      if (iszipped) {
        push(boost::iostreams::gzip_decompressor());
      }
      push(uncompressed_infile);

      if (!uncompressed_infile || !*this) {
        throw file_error(path(name), true);
      }

      return true;
    }

    virtual ~izfile() {
      //must remove streams before deallocating
      while (!empty())
        pop();
    }
};

//dkoes - wrapper for an output file that will be gzipped if the file
//name ends in .gz
class ozfile : public boost::iostreams::filtering_stream<
    boost::iostreams::output> {

    std::ofstream uncompressed_outfile;
  public:

    ozfile() {
    }

    ozfile(const path& name) {
      open(name);
    }

    //opens file name, with gzip filter if name ends with .gz
    //return non-gz extension
    std::string open(const path& name) {
      using namespace boost::filesystem;
      uncompressed_outfile.open(name.c_str());
      if (!uncompressed_outfile) throw file_error(name, false);

      std::string ext = boost::filesystem::extension(name);
      //should we gzip?
      if (ext == ".gz") {
        ext = extension(basename(name));
        push(boost::iostreams::gzip_compressor());
      }
      push(uncompressed_outfile);
      if (!(*this)) throw file_error(name, false);
      return ext;
    }

    virtual ~ozfile() {
      //must remove streams before deallocating
      while (!empty())
        pop();
    }
};

#endif
