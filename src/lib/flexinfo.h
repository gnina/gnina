#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include "obmolopener.h"
#include "box.h"
#include "tee.h"

#ifndef SMINA_FLEXINFO_H
#define SMINA_FLEXINFO_H

/* Store information for identifying flexible residues in receptor and
 * provide routines for extracting these residues as needed.
 */
class FlexInfo
{
	double flex_dist;
	boost::unordered_set< std::pair<char, int> > residues;
	OpenBabel::OBMol distligand;
	tee& log;
public:
	FlexInfo(const std::string& flexres, double flexdist, const std::string& ligand, tee& l);
	bool hasContent() const
	{
		return residues.size() >0 || flex_dist > 0;
	}

	void extractFlex(OpenBabel::OBMol& receptor, OpenBabel::OBMol& rigid, std::string& flexpdbqt);

};

#endif /* SMINA_FLEXINFO_H */
