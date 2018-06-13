// simple struct for storing info on two atoms
#ifndef INTERACTING_PAIR_H
#define INTERACTING_PAIR_H

#include "atom_constants.h"

struct interacting_pair {
    smt t1;
    smt t2;
    sz a;
    sz b;
    interacting_pair()
        : t1(smina_atom_type::Hydrogen), t2(smina_atom_type::Hydrogen), a(0),
            b(0) {
    }
    interacting_pair(smt t1_, smt t2_, sz a_, sz b_)
        : t1(t1_), t2(t2_), a(a_), b(b_) {
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned version) {
      ar & t1;
      ar & t2;
      ar & a;
      ar & b;
    }
};

#endif
