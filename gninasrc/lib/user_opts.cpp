/*
 * user_opts.cpp
 *
 *  Created on: Jul 3, 2020
 *      Author: dkoes
 */

#include "user_opts.h"
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

using namespace boost::algorithm;
using namespace boost::program_options;

//enum options and their parsers
std::istream& operator>>(std::istream& in, ApproxType& type) {
  using namespace boost::program_options;

  std::string token;
  in >> token;
  if (token == "spline")
    type = SplineApprox;
  else if (token == "linear")
    type = LinearApprox;
  else if (token == "exact")
    type = Exact;
  else if (token == "gpu")
    type = GPU;
  else
    throw validation_error(validation_error::invalid_option_value);
  return in;
}

//for reading in as a commandline option
std::istream& operator>>(std::istream &in, pose_sort_order &sort_order)
{
  std::string token;
  in >> token;
  to_lower(token);
  if (token == "cnnscore")
    sort_order = CNNscore;
  else if (token == "cnnaffinity")
    sort_order = CNNaffinity;
  else if (token == "energy" || token == "vina")
    sort_order = Energy;
  else
    throw validation_error(validation_error::invalid_option_value);
  return in;
}

std::istream& operator>>(std::istream &in, cnn_scoring_level &cnn_level)
{
  std::string token;
  in >> token;
  to_lower(token);
  if (token == "all" || token == "docking")
    cnn_level = CNNall;
  else if (token == "rescore")
    cnn_level = CNNrescore;
  else if (starts_with(token, "refine") || starts_with(token,"min"))
    cnn_level = CNNrefinement;
  else if(token == "none" || token == "no")
    cnn_level = CNNnone;
  else
    throw validation_error(validation_error::invalid_option_value);
  return in;
}



