/*
 * builtinscoring.h
 *
 *  Created on: Oct 15, 2015
 *      Author: dkoes
 *
 *  Singleton class for setting up builtin scoring terms
 */

#include <vector>
#include <boost/unordered_map.hpp>
#include "custom_terms.h"

class builtin_scoring {
    struct singleterm {
        const char* term;
        double weight;

        singleterm()
            : term(NULL), weight(0) {
        }
        singleterm(const char* t, double w)
            : term(t), weight(w) {
        }
    };

    //map from names to sets of terms
    typedef boost::unordered_map<std::string, std::vector<singleterm> > funcmap;
    //optional - map from names to different atom type data
    typedef boost::unordered_map<std::string, const smina_atom_type::info* > datamap;
    funcmap functions;
    datamap atomdata;

    void add(const std::string& name, const char *term, double w) {
      functions[name].push_back(singleterm(term, w));
    }

    void addparams(const std::string& name, const smina_atom_type::info* data)
    {
      atomdata[name] = data;
    }

  public:
    builtin_scoring();

    //set t to the terms specified by builtin name
    //return true if successful
    bool set(custom_terms& t, const std::string& name) {
      if (functions.count(name) == 0) return false;

      std::vector<singleterm>& terms = functions[name];
      for (unsigned i = 0, n = terms.size(); i < n; i++) {
        t.add(terms[i].term, terms[i].weight);
      }

      if(atomdata.count(name))
      {
        //not an error if doesn't exist - just use default atom parameters
        //but if it does exist, set new params
        const smina_atom_type::info *data = atomdata[name];
        assert(data);
        for(size_t i = 0u; i < smina_atom_type::NumTypes; ++i)
          smina_atom_type::data[i] = data[i];
      }
      return true;
    }

    void print_functions(std::ostream& out);
};

extern builtin_scoring builtin_scoring_functions;
