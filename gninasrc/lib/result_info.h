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
#include "cnn_scorer.h"

// this class holds the contents of the result of a minization/docking
// it handles outputing the molecular data in the appropriate format
class result_info {
    fl energy;
    fl cnnscore;
    fl cnnaffinity;
    fl rmsd;
    std::string molstr;
    std::string flexstr;
    std::string atominfo;
    std::string name;
    bool sdfvalid;

  public:
    result_info()
        : energy(0), cnnscore(-1), cnnaffinity(0), rmsd(-1), sdfvalid(false) {
    }
    result_info(fl e, fl c, fl ca, fl r, const model& m)
        : energy(e), cnnscore(c), cnnaffinity(ca), rmsd(r), sdfvalid(false) {
      setMolecule(m);
    }

    //set the molecular data using the current conformation of model m
    void setMolecule(const model& m);

    //write a table (w/header) of per atom values to out
    void writeAtomValues(std::ostream& out, const weighted_terms *wt) const;

    //computes per-atom term values and formats them into the atominfo string
    void setAtomValues(const model& m, const weighted_terms *wt);

    void write(std::ostream& out, std::string& ext, bool include_atom_terms,
        const weighted_terms *wt = NULL, int modelnum = 0);

    void writeFlex(std::ostream& out, std::string& ext, int modelnum = 0);
};

#endif /* RESULT_INFO_H_ */
