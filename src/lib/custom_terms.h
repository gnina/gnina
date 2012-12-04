/*
 * Copyright 2012 David Koes and University of Pittsburgh
 *
 * This implements a customizable scoring function class that can be
 * initialized using a user provided file.
 */
#ifndef VINA_CUSTOM_TERMS_H
#define VINA_CUSTOM_TERMS_H


#include "terms.h"
#include "int_pow.h"
#include <string>

struct custom_terms : public terms {
	custom_terms(); //creates empty set
	add(const std::string& name, double weight);
};


#endif
