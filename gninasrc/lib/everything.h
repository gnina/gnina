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

/*
 * SMINA NOTICE
 * dkoes - if you add a term here, it must also define a unique regular expression
 * that matches a string description for use in custom scoring
 * createFrom should take a string description and return a newly created and
 * allocated term that matches the string (or null if no match, or throw an error
 * if match but poorly parameterized).
 * As a convenience, provide reasonable default initializers.
 *
 * I'm not a big fan of having objects being their own factories, but want to keep
 * as much implementation as possible contained within the new term's class definition.
 * Suggestions for a better design are welcome.
 *
 * You must also register your term with the term_creators constructor
 *
 */

#ifndef VINA_EVERYTHING_H
#define VINA_EVERYTHING_H

#include "terms.h"
#include "int_pow.h"
#include <boost/lexical_cast.hpp>
#include <iostream>

inline fl gaussian(fl x, fl width) {
  return std::exp(-sqr(x / width));
}

inline fl smooth_div(fl x, fl y) {
  if (std::abs(x) < epsilon_fl) return 0;
  if (std::abs(y) < epsilon_fl) return ((x * y > 0) ? max_fl : -max_fl); // FIXME I hope -max_fl does not become NaN
  return x / y;
}

// distance_additive terms

template<unsigned power>
struct electrostatic : public charge_dependent {
    fl cap;

    electrostatic(fl cap_ = 100, fl cutoff_ = 8)
        : charge_dependent(cutoff_), cap(cap_) {
      name = std::string("electrostatic(i=") + to_string(power) + ",_^="
          + to_string(cap) + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("electrostatic\\(i=(\\S+),_\\^=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);

    }

    result_components eval_components(smt t1, smt t2, fl r) const {
      result_components comp;
      fl tmp = int_pow<power>(r);
      if (tmp < epsilon_fl) //avoid divide by zero
        comp[result_components::ABChargeDependent] = cap;
      else
        comp[result_components::ABChargeDependent] = (std::min)(cap, 1 / tmp);
      return comp;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl i = boost::lexical_cast<fl>(match[1]);
      fl cap = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      if (i == 1)
        return new electrostatic<1>(cap, c);
      else
        if (i == 2)
          return new electrostatic<2>(cap, c);
        else
          throw scoring_function_error(desc, "Invalid exponent: 1 or 2 only");
    }

};

struct ad4_solvation : public charge_dependent {
    fl solvation_q;
    fl desolvation_sigma;
    ad4_solvation(fl desolvation_sigma_ = 3.6, fl solvation_q_ = 0.01097,
        fl cutoff_ = 8)
        : charge_dependent(cutoff_), solvation_q(solvation_q_),
            desolvation_sigma(desolvation_sigma_) {
      name = std::string("ad4_solvation(d-sigma=")
          + to_string(desolvation_sigma) + ",_s/q=" + to_string(solvation_q)
          + +",_c=" + to_string(cutoff) + ")";
      rexpr.assign("ad4_solvation\\(d-sigma=(\\S+),_s/q=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);

    }

    result_components eval_components(smt t1, smt t2, fl r) const {
      result_components ret;
      fl solv1 = solvation_parameter(t1);
      fl solv2 = solvation_parameter(t2);

      fl volume1 = ad_volume(t1);
      fl volume2 = ad_volume(t2);

      fl mysolv = solvation_q;
      fl distfactor = std::exp(-sqr(r / (2 * desolvation_sigma)));

      ret[result_components::TypeDependentOnly] += solv1 * volume2 * distfactor;
      ret[result_components::AbsAChargeDependent] += mysolv * volume2
          * distfactor;

      ret[result_components::TypeDependentOnly] += solv2 * volume1 * distfactor;
      ret[result_components::AbsBChargeDependent] += mysolv * volume1
          * distfactor;

      return ret;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl sigma = boost::lexical_cast<fl>(match[1]);
      fl w = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      return new ad4_solvation(sigma, w, c);
    }
};

inline fl optimal_distance(smt xs_t1, smt xs_t2) {
  return xs_radius(xs_t1) + xs_radius(xs_t2);
}

struct gauss : public charge_independent {
    fl offset; // added to optimal distance
    fl width;
    gauss(fl offset_ = 0, fl width_ = 0.5, fl cutoff_ = 8)
        : charge_independent(cutoff_), offset(offset_), width(width_) {
      name = std::string("gauss(o=") + to_string(offset) + ",_w="
          + to_string(width) + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("gauss\\(o=(\\S+),_w=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);

    }

    using charge_independent::eval;
    virtual fl eval(smt t1, smt t2, fl r) const {
      return gaussian(r - (optimal_distance(t1, t2) + offset), width);
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl o = boost::lexical_cast<fl>(match[1]);
      fl w = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      return new gauss(o, w, c);
    }
};

struct repulsion : public charge_independent {
    fl offset; // added to vdw
    repulsion(fl offset_ = 0, fl cutoff_ = 8)
        : charge_independent(cutoff_), offset(offset_) {
      name = std::string("repulsion(o=") + to_string(offset) + ",_c="
          + to_string(cutoff) + ")";
      rexpr.assign("repulsion\\(o=(\\S+),_c=(\\S+)\\)", boost::regex::perl);

    }
    using charge_independent::eval;
    virtual fl eval(smt t1, smt t2, fl r) const {
      fl d = r - (optimal_distance(t1, t2) + offset);
      if (d > 0) return 0;
      return d * d;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl o = boost::lexical_cast<fl>(match[1]);
      fl c = boost::lexical_cast<fl>(match[2]);
      return new repulsion(o, c);
    }
};

inline fl slope_step(fl x_bad, fl x_good, fl x) {
  if (x_bad < x_good) {
    if (x <= x_bad) return 0;
    if (x >= x_good) return 1;
  } else {
    if (x >= x_bad) return 0;
    if (x <= x_good) return 1;
  }
  return (x - x_bad) / (x_good - x_bad);
}

struct hydrophobic : public charge_independent {
    fl good;
    fl bad;
    hydrophobic(fl good_ = 0.5, fl bad_ = 1.5, fl cutoff_ = 8)
        : charge_independent(cutoff_), good(good_), bad(bad_) {
      name = "hydrophobic(g=" + to_string(good) + ",_b=" + to_string(bad)
          + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("hydrophobic\\(g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }

    using charge_independent::eval;
    virtual fl eval(smt t1, smt t2, fl r) const {
      //std::cout << "HYDRO " << t1 << " " << t2 << " " << r << " " << slope_step(bad, good, r - optimal_distance(t1, t2)) << "\n";
      if (xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2))
        return slope_step(bad, good, r - optimal_distance(t1, t2));
      else
        return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl g = boost::lexical_cast<fl>(match[1]);
      fl b = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      return new hydrophobic(g, b, c);
    }
};

struct non_hydrophobic : public charge_independent {
    fl good;
    fl bad;
    non_hydrophobic(fl good_ = 0.5, fl bad_ = 1.5, fl cutoff_ = 8)
        : charge_independent(cutoff_), good(good_), bad(bad_) {
      name = "non_hydrophobic(g=" + to_string(good) + ",_b=" + to_string(bad)
          + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("non_hydrophobic\\(g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);

    }

    using charge_independent::eval;
    virtual fl eval(smt t1, smt t2, fl r) const {
      if (!xs_is_hydrophobic(t1) && !xs_is_hydrophobic(t2))
        return slope_step(bad, good, r - optimal_distance(t1, t2));
      else
        return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl g = boost::lexical_cast<fl>(match[1]);
      fl b = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      return new non_hydrophobic(g, b, c);
    }
};

template<unsigned n, unsigned m>
void find_vdw_coefficients(fl position, fl depth, fl& c_n, fl& c_m) {
  BOOST_STATIC_ASSERT(n != m);
  c_n = int_pow<n>(position) * depth * m / (fl(n) - fl(m));
  c_m = int_pow<m>(position) * depth * n / (fl(m) - fl(n));
}

template<unsigned i, unsigned j>
struct vdw : public charge_independent {
    fl smoothing;
    fl cap;
    vdw(fl smoothing_ = 1, fl cap_ = 100, fl cutoff_ = 8)
        : charge_independent(cutoff_), smoothing(smoothing_), cap(cap_) {
      name = "vdw(i=" + to_string(i) + ",_j=" + to_string(j) + ",_s="
          + to_string(smoothing) + ",_^=" + to_string(cap) + ",_c="
          + to_string(cutoff) + ")";
      rexpr.assign(
          "vdw\\(i=(\\S+),_j=(\\S+),_s=(\\S+),_\\^=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);

    }

    using charge_independent::eval;
    virtual fl eval(smt t1, smt t2, fl r) const {
      fl d0 = optimal_distance(t1, t2);
      fl depth = 1;
      fl c_i = 0;
      fl c_j = 0;
      find_vdw_coefficients<i, j>(d0, depth, c_i, c_j);
      if (r > d0 + smoothing)
        r -= smoothing;
      else
        if (r < d0 - smoothing)
          r += smoothing;
        else
          r = d0;

      fl r_i = int_pow<i>(r);
      fl r_j = int_pow<j>(r);
      if (r_i > epsilon_fl && r_j > epsilon_fl)
        return (std::min)(cap, c_i / r_i + c_j / r_j);
      else
        return cap;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl vi = boost::lexical_cast<fl>(match[1]);
      fl vj = boost::lexical_cast<fl>(match[2]);
      fl s = boost::lexical_cast<fl>(match[3]);
      fl cap = boost::lexical_cast<fl>(match[4]);
      fl c = boost::lexical_cast<fl>(match[5]);
      if (vi == 4.0 && vj == 8)
        return new vdw<4, 8>(s, cap, c);
      else
        if (vi == 6 && vj == 12)
          return new vdw<6, 12>(s, cap, c);
        else
          throw scoring_function_error(desc,
              "Unsupported LJ exponents: try <4,8> or <6,12>.");
    }
};

/* A 10-12 LJ potential */
struct non_dir_h_bond_lj : public charge_independent {
    fl offset;
    fl cap;
    non_dir_h_bond_lj(fl offset_ = -0.7, fl cap_ = 100, fl cutoff_ = 8)
        : charge_independent(cutoff_), offset(offset_), cap(cap_) {
      name = std::string("non_dir_h_bond_lj(o=") + to_string(offset) + ",_^="
          + to_string(cap) + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("non_dir_h_bond_lj\\(o=(\\S+),_\\^=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }

    using charge_independent::eval;
    fl eval(smt t1, smt t2, fl r) const {
      if (xs_h_bond_possible(t1, t2)) {
        fl d0 = optimal_distance(t1, t2) + offset;
        fl depth = 5;
        fl c_i = 0;
        fl c_j = 0;
        find_vdw_coefficients<10, 12>(d0, depth, c_i, c_j);

        fl r_i = int_pow<10>(r);
        fl r_j = int_pow<12>(r);
        if (r_i > epsilon_fl && r_j > epsilon_fl)
          return (std::min)(cap, c_i / r_i + c_j / r_j);
        else
          return cap;
      }
      return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl o = boost::lexical_cast<fl>(match[1]);
      fl cap = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      return new non_dir_h_bond_lj(o, cap, c);
    }
};

/* This mimics repulsion, but only between polar atoms that can't possibly hydrogen bond */
struct non_dir_anti_h_bond_quadratic : public charge_independent {
    fl offset;
    non_dir_anti_h_bond_quadratic(fl offset_ = 0, fl cutoff_ = 8)
        : charge_independent(cutoff_), offset(offset_) {
      name = std::string("non_dir_anti_h_bond_quadratic(o=") + to_string(offset)
          + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("non_dir_anti_h_bond_quadratic\\(o=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }

    using charge_independent::eval;
    fl eval(smt t1, smt t2, fl r) const {
      if (xs_anti_h_bond(t1, t2)) {
        fl d = r - (optimal_distance(t1, t2) + offset);
        if (d > 0) return 0;
        return d * d;
      }
      return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl o = boost::lexical_cast<fl>(match[1]);
      fl c = boost::lexical_cast<fl>(match[2]);
      return new non_dir_anti_h_bond_quadratic(o, c);
    }
};

/* Quadratic potential (see repulsion) between donor atoms */
struct donor_donor_quadratic : public charge_independent {
    fl offset;
    donor_donor_quadratic(fl offset_ = 0, fl cutoff_ = 8)
        : charge_independent(cutoff_), offset(offset_) {
      name = std::string("donor_donor_quadratic(o=") + to_string(offset)
          + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("donor_donor_quadratic\\(o=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }

    using charge_independent::eval;
    fl eval(smt t1, smt t2, fl r) const {
      if (xs_is_donor(t1) && xs_is_donor(t2)) {
        fl d = r - (optimal_distance(t1, t2) + offset);
        if (d > 0) return 0;
        return d * d;
      }
      return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl o = boost::lexical_cast<fl>(match[1]);
      fl c = boost::lexical_cast<fl>(match[2]);
      return new donor_donor_quadratic(o, c);
    }
};

/* Quadratic potential (see repulsion) between acceptor atoms */
struct acceptor_acceptor_quadratic : public charge_independent {
    fl offset;
    acceptor_acceptor_quadratic(fl offset_ = 0, fl cutoff_ = 8)
        : charge_independent(cutoff_), offset(offset_) {
      name = std::string("acceptor_acceptor_quadratic(o=") + to_string(offset)
          + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("acceptor_acceptor_quadratic\\(o=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }

    using charge_independent::eval;
    fl eval(smt t1, smt t2, fl r) const {
      if (xs_is_acceptor(t1) && xs_is_acceptor(t2)) {
        fl d = r - (optimal_distance(t1, t2) + offset);
        if (d > 0) return 0;
        return d * d;
      }
      return 0;
    }
    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl o = boost::lexical_cast<fl>(match[1]);
      fl c = boost::lexical_cast<fl>(match[2]);
      return new acceptor_acceptor_quadratic(o, c);
    }
};

//classic Vina hbond term
struct non_dir_h_bond : public charge_independent {
    fl good;
    fl bad;
    non_dir_h_bond(fl good_ = -0.7, fl bad_ = 0, fl cutoff_ = 8)
        : charge_independent(cutoff_), good(good_), bad(bad_) {
      name = std::string("non_dir_h_bond(g=") + to_string(good) + ",_b="
          + to_string(bad) + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign("non_dir_h_bond\\(g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }
    using charge_independent::eval;
    fl eval(smt t1, smt t2, fl r) const {
      if (xs_h_bond_possible(t1, t2))
        return slope_step(bad, good, r - optimal_distance(t1, t2));
      return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;

      fl g = boost::lexical_cast<fl>(match[1]);
      fl b = boost::lexical_cast<fl>(match[2]);
      fl c = boost::lexical_cast<fl>(match[3]);
      return new non_dir_h_bond(g, b, c);
    }
};

/* dkoes - atom type functions.  These are generic functions that only apply to
 * the specified atom type pair.  A variety of distance dependent functions are
 * provided. Atom types are specifed by strings so we can pretty print them.
 */
struct atom_type_base : public charge_independent {
    //base type for atom type functions
    std::string name1, name2;
    smt t1, t2; //atom type pair

    using charge_independent::eval;

    atom_type_base(const std::string& n1, const std::string& n2, fl cutoff_ = 8)
        : charge_independent(cutoff_), name1(n1), name2(n2),
            t1(string_to_smina_type(n1)), t2(string_to_smina_type(n2)) {
      if (name1.length() > 0 && t1 == smina_atom_type::NumTypes) //ignore default empty
      throw scoring_function_error(name1, "Invalid atom type: ");
      if (name2.length() > 0 && t2 == smina_atom_type::NumTypes)
        throw scoring_function_error(name2, "Invalid atom type: ");
    }

  protected:

    bool types_match(smt t1_, smt t2_) const {
      //match any order
      return (t1_ == t1 && t2_ == t2) || (t1_ == t2 && t2_ == t1);
    }
};

/* inverse power potential (see electrostatic) between atom types */
template<unsigned power>
struct atom_type_inverse_power : public atom_type_base {
    fl cap;
    atom_type_inverse_power(const std::string& n1 = "", const std::string& n2 =
        "", fl cap_ = 100, fl cutoff_ = 8)
        : atom_type_base(n1, n2, cutoff_), cap(cap_) {
      name = std::string("atom_type_inverse_power(t1=") + name1 + ",t2=" + name2
          + ",i=" + to_string(power) + ",_^=" + to_string(cap) + ",_c="
          + to_string(cutoff) + ")";
      rexpr.assign(
          "atom_type_inverse_power\\(t1=(\\S+),t2=(\\S+),i=(\\S+),_\\^=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }
    using charge_independent::eval;

    fl eval(smt T1, smt T2, fl r) const {
      if (types_match(T1, T2)) {
        fl tmp = int_pow<power>(r);
        if (tmp < epsilon_fl) //avoid divide by zero
          return cap;
        else
          return (std::min)(cap, 1 / tmp);
      }
      return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;
      std::string n1 = match[1];
      std::string n2 = match[2];

      fl i = boost::lexical_cast<fl>(match[3]);
      fl cap = boost::lexical_cast<fl>(match[4]);
      fl c = boost::lexical_cast<fl>(match[5]);
      if (i == 1)
        return new atom_type_inverse_power<1>(n1, n2, cap, c);
      else {
        if (i == 2)
          return new atom_type_inverse_power<2>(n1, n2, cap, c);
        else
          throw scoring_function_error(desc, "Invalid exponent: 1 or 2 only");
      }
    }
};

/* Gaussian potential (see gauss) between atom types */
struct atom_type_gaussian : public atom_type_base {
    fl width, offset;
    atom_type_gaussian(const std::string& n1 = "", const std::string& n2 = "",
        fl o = 0, fl w = 0, fl cutoff_ = 8)
        : atom_type_base(n1, n2, cutoff_), width(w), offset(o) {
      name = std::string(
          "atom_type_gaussian(t1=" + name1 + ",t2=" + name2 + ",o=")
          + to_string(offset) + ",_w=" + to_string(width) + ",_c="
          + to_string(cutoff) + ")";
      rexpr.assign(
          "atom_type_gaussian\\(t1=(\\S+),t2=(\\S+),o=(\\S+),_w=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }
    using charge_independent::eval;
    virtual fl eval(smt T1, smt T2, fl r) const {
      if (types_match(T1, T2)) {
        return gaussian(r - (optimal_distance(t1, t2) + offset), width);
      }
      return 0;
    }

    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;
      std::string n1 = match[1];
      std::string n2 = match[2];
      fl o = boost::lexical_cast<fl>(match[3]);
      fl w = boost::lexical_cast<fl>(match[4]);
      fl c = boost::lexical_cast<fl>(match[5]);
      return new atom_type_gaussian(n1, n2, o, w, c);
    }
};

/* 6-12 LJ potential (see vdw) between atom types */
struct atom_type_lennard_jones: public atom_type_base
{
    fl optimal_distance;
    fl cap;
    atom_type_lennard_jones(const std::string& n1="", const std::string& n2="", fl o=0, fl c=100, fl cutoff_=8) :
        atom_type_base(n1, n2, cutoff_), optimal_distance(o), cap(c)
    {
        name = std::string("atom_type_lennard_jones(t1="+name1+",t2="+name2+",o=")
                + to_string(optimal_distance)+",_^="+to_string(cap)+",_c=" + to_string(cutoff) + ")";
        rexpr.assign("atom_type_lennard_jones\\(t1=(\\S+),t2=(\\S+),o=(\\S+),_\\^=(\\S+),_c=(\\S+)\\)",boost::regex::perl);
    }

    virtual fl eval(smt T1, smt T2, fl r) const
    {
        fl d0 = optimal_distance;
        fl depth = 1;
        fl c_i = 0;
        fl c_j = 0;
        find_vdw_coefficients<6, 12>(d0, depth, c_i, c_j);

        fl r_i = int_pow<6>(r);
        fl r_j = int_pow<12>(r);
        if (r_i > epsilon_fl && r_j > epsilon_fl)
            return (std::min)(cap, c_i / r_i + c_j / r_j);
        else
            return cap;
    }

    virtual term* createFrom(const std::string& desc) const {
        boost::smatch match;
        if(!regex_match(desc, match, rexpr))
            return NULL;
        std::string n1 = match[1];
        std::string n2 = match[2];
        fl o = boost::lexical_cast<fl>(match[3]);
        fl cap = boost::lexical_cast<fl>(match[4]);
        fl c = boost::lexical_cast<fl>(match[5]);
        return new atom_type_lennard_jones(n1,n2,o,cap,c);
    }
};

/* Linear potential (see hbond) between atom types */
struct atom_type_linear : public atom_type_base {
    fl good, bad;
    atom_type_linear(const std::string& n1 = "", const std::string& n2 = "",
        fl good_ = 0, fl bad_ = 0, fl cutoff_ = 8)
        : atom_type_base(n1, n2, cutoff_), good(good_), bad(bad_) {
      name = std::string(
          "atom_type_linear(t1=" + name1 + ",t2=" + name2 + ",g=")
          + to_string(good) + ",_b=" + to_string(bad) + ",_c="
          + to_string(cutoff) + ")";
      rexpr.assign(
          "atom_type_linear\\(t1=(\\S+),t2=(\\S+),g=(\\S+),_b=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }
    using charge_independent::eval;
    virtual fl eval(smt T1, smt T2, fl r) const {
      if (types_match(T1, T2)) {
        return slope_step(bad, good, r - optimal_distance(t1, t2));
      }
      return 0;
    }
    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;
      std::string n1 = match[1];
      std::string n2 = match[2];
      fl g = boost::lexical_cast<fl>(match[3]);
      fl b = boost::lexical_cast<fl>(match[4]);
      fl c = boost::lexical_cast<fl>(match[5]);
      return new atom_type_linear(n1, n2, g, b, c);
    }
};

/* Quadratic potential (see repulsion) between atom types */
struct atom_type_quadratic : public atom_type_base {
    fl offset;
    atom_type_quadratic(const std::string& n1 = "", const std::string& n2 = "",
        fl offset_ = 0, fl cutoff_ = 8)
        : atom_type_base(n1, n2, cutoff_), offset(offset_) {
      name = std::string(
          "atom_type_quadratic(t1=" + name1 + ",t2=" + name2 + ",o=")
          + to_string(offset) + ",_c=" + to_string(cutoff) + ")";
      rexpr.assign(
          "atom_type_quadratic\\(t1=(\\S+),t2=(\\S+),o=(\\S+),_c=(\\S+)\\)",
          boost::regex::perl);
    }
    using charge_independent::eval;
    virtual fl eval(smt T1, smt T2, fl r) const {
      if (types_match(T1, T2)) {
        fl d = r - (optimal_distance(t1, t2) + offset);
        if (d > 0) return 0;
        return d * d;
      }
      return 0;
    }
    virtual term* createFrom(const std::string& desc) const {
      boost::smatch match;
      if (!regex_match(desc, match, rexpr)) return NULL;
      std::string n1 = match[1];
      std::string n2 = match[2];
      fl o = boost::lexical_cast<fl>(match[3]);
      fl c = boost::lexical_cast<fl>(match[4]);
      return new atom_type_quadratic(n1, n2, o, c);
    }
};

inline fl read_iterator(flv::const_iterator& i) {
  fl x = *i;
  ++i;
  return x;
}

fl smooth_div(fl x, fl y);

struct num_tors_add : public conf_independent {
    num_tors_add() {
      name = "num_tors_add";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      //fl w = 0.1 * read_iterator(i); // [-1 .. 1]
      fl w = read_iterator(i); // FIXME?
      return x + w * in.num_tors;
    }

    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;

      return new num_tors_add();
    }
};

struct num_tors_sqr : public conf_independent {
    num_tors_sqr() {
      name = "num_tors_sqr";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 0.1 * read_iterator(i); // [-1 .. 1]
      fl add = w * sqr(fl(in.num_tors)) / 5;
      fl ret = x + add;
      return ret;
    }

    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_tors_sqr();
    }
};

struct num_tors_sqrt : public conf_independent {
    num_tors_sqrt() {
      name = "num_tors_sqrt";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 0.1 * read_iterator(i); // [-1 .. 1]
      return x + w * std::sqrt(fl(in.num_tors)) / sqrt(5.0);
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_tors_sqrt();
    }
};

struct num_tors_div : public conf_independent {
    num_tors_div() {
      name = "num_tors_div";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 0.1 * (read_iterator(i) + 1); // w is in [0..0.2]
      //std::cout << "Num_tors_factor: " << test << std::endl;
      return smooth_div(x, 1 + w * in.num_tors / 5.0);
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_tors_div();
    }
};

//just divide energy by 1+w*num_tors as opposed to the more complicated vina formula
struct num_tors_div_simple : public conf_independent {
    num_tors_div_simple() {
      name = "num_tors_div_simple";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = read_iterator(i);
      //std::cout << "Num_tors_factor: " << test << std::endl;
      return smooth_div(x, 1 + w * in.num_tors);
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_tors_div_simple();
    }
};

struct ligand_length : public conf_independent {
    ligand_length() {
      name = "ligand_length";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = read_iterator(i);
      return x + w * in.ligand_lengths_sum;
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new ligand_length();
    }
};

struct num_ligands : public conf_independent {
    num_ligands() {
      name = "num_ligands";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 1 * read_iterator(i); // w is in [-1.. 1]
      return x + w * in.num_ligands;
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_ligands();
    }
};

struct num_heavy_atoms_div : public conf_independent {
    num_heavy_atoms_div() {
      name = "num_heavy_atoms_div";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 0.05 * read_iterator(i);
      return smooth_div(x, 1 + w * in.num_heavy_atoms);
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_heavy_atoms_div();
    }
};

struct num_heavy_atoms : public conf_independent {
    num_heavy_atoms() {
      name = "num_heavy_atoms";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 0.05 * read_iterator(i);
      return x + w * in.num_heavy_atoms;
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_heavy_atoms();
    }
};

struct num_hydrophobic_atoms : public conf_independent {
    num_hydrophobic_atoms() {
      name = "num_hydrophobic_atoms";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = 0.05 * read_iterator(i);
      return x + w * in.num_hydrophobic_atoms;
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new num_hydrophobic_atoms();
    }
};

struct constant_term : public conf_independent {
    constant_term() {
      name = "constant_term";
      rexpr.assign(name);
    }
    sz size() const {
      return 1;
    }
    fl eval(const conf_independent_inputs& in, fl x,
        flv::const_iterator& i) const {
      fl w = read_iterator(i);
      return x + w;
    }
    virtual term* createFrom(const std::string& desc) const {
      if (!regex_match(desc, rexpr)) return NULL;
      return new constant_term();
    }
};

//vector of terms
//ADD ALL TERMS HERE
struct term_creators : public std::vector<term*> {
    term_creators() {
      push_back(new electrostatic<2>());
      push_back(new ad4_solvation());
      push_back(new gauss());
      push_back(new repulsion());
      push_back(new hydrophobic());
      push_back(new non_hydrophobic());
      push_back(new vdw<6, 12>());
      push_back(new non_dir_h_bond_lj());
      push_back(new non_dir_anti_h_bond_quadratic());
      push_back(new non_dir_h_bond());
      push_back(new acceptor_acceptor_quadratic());
      push_back(new donor_donor_quadratic());

      push_back(new atom_type_gaussian());
      push_back(new atom_type_linear());
      push_back(new atom_type_quadratic());
      push_back(new atom_type_inverse_power<0>());
      push_back(new atom_type_lennard_jones());

      push_back(new num_tors_add());
      push_back(new num_tors_sqr());
      push_back(new num_tors_sqrt());
      push_back(new num_tors_div());
      push_back(new num_tors_div_simple());
      push_back(new ligand_length());
      push_back(new num_ligands());
      push_back(new num_heavy_atoms_div());
      push_back(new num_heavy_atoms());
      push_back(new num_hydrophobic_atoms());
      push_back(new constant_term());
    }

    virtual ~term_creators() {
      for (unsigned i = 0, n = size(); i < n; i++) {
        delete (*this)[i];
      }
      clear();
    }
};

#endif
