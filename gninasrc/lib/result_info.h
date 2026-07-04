/*
 * result_info.h
 *
 *  Created on: Jun 10, 2014
 *      Author: dkoes
 */

#ifndef RESULT_INFO_H_
#define RESULT_INFO_H_

#include <iostream>
#include <string>
#include "model.h"
#include "weighted_terms.h"

// this class holds the contents of the result of a minization/docking
// it handles outputing the molecular data in the appropriate format
class result_info {
    fl energy = 0;
    fl cnnscore = -1;
    fl cnnaffinity = 0;
    fl cnnvariance = 0;         // affinity (regression) variance across ensemble × rotations
    fl cnnscore_variance = 0;   // pose-score (classification) variance across ensemble × rotations
    fl rmsd = -1;
    std::string molstr;
    std::string flexstr;
    std::string atominfo;
    std::string name;
    bool sdfvalid = false;;
    bool flexsdfvalid = false;

  public:
    result_info() { }
    result_info(fl e, fl c, fl ca, fl cv, fl csv, fl r, const model& m)
        : energy(e), cnnscore(c), cnnaffinity(ca), cnnvariance(cv), cnnscore_variance(csv), rmsd(r), sdfvalid(false) {
      setMolecule(m);
    }

    //set the molecular data using the current conformation of model m
    void setMolecule(const model& m);

    //write a table (w/header) of per atom values to out
    void writeAtomValues(std::ostream& out, const weighted_terms *wt) const;

    //computes per-atom term values and formats them into the atominfo string
    void setAtomValues(const model& m, const weighted_terms *wt);

    void write(std::ostream& out, const std::string& ext, bool include_atom_terms,
        const weighted_terms *wt = NULL, int modelnum = 0) ;

    void writeFlex(std::ostream& out, const std::string& ext, int modelnum = 0);

    fl getEnergy() const { return energy;}
    fl getCNNScore() const { return cnnscore; }
    fl getCNNAffinity() const { return cnnaffinity;}
};

#endif /* RESULT_INFO_H_ */
