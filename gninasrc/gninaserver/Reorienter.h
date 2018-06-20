/*
 * ReOrienter.h
 *
 *  Created on: Jun 11, 2014
 *      Author: dkoes
 *
 *  A class for transforming cartesian coordinates
 */

#ifndef REORIENTER_H_
#define REORIENTER_H_

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include "common.h"
#include "parsing.h"

using namespace Eigen;
using namespace std;

class Reorienter {
    Matrix3d rotation;
    Vector3d translation;

  public:
    Reorienter()
        : rotation(Matrix3d::Identity()), translation(Vector3d::Zero()) {

    }

    void reorient(vec& pt) const {
      Vector3d v(pt[0], pt[1], pt[2]);
      v = rotation * v + translation;
      pt = vec(v[0], v[1], v[2]);
    }

    //modify points by rot/trans
    void reorient(vecv& pts) const {
      Vector3d pt;
      for (unsigned i = 0, n = pts.size(); i < n; i++) {
        reorient(pts[i]);
      }
    }

    void reorient(parsing_struct::node& node) const {
      reorient(node.a.coords);
      for (unsigned i = 0, n = node.ps.size(); i < n; i++) {
        reorient(node.ps[i]);
      }
    }

    //recursively modify p
    void reorient(parsing_struct& p) const {
      for (unsigned i = 0, n = p.atoms.size(); i < n; i++) {
        reorient(p.atoms[i]);
      }
    }

    void read(istream& in) {
      in.read((char*) rotation.data(), 9 * sizeof(double));
      //eigen uses column major ordering!
      rotation.transposeInPlace();
      in.read((char*) translation.data(), 3 * sizeof(double));
    }
};

#endif /* REORIENTER_H_ */
