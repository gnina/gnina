#include <openbabel/mol.h>
#include <cmath> // for ceila
#include "common.h"

#ifndef SMINA_BOX_H
#define SMINA_BOX_H

/* store a bounding box */
struct Box
{
	fl min_x, max_x;
	fl min_y, max_y;
	fl min_z, max_z;


	Box(): min_x(HUGE_VAL), max_x(-HUGE_VAL), min_y(HUGE_VAL), max_y(-HUGE_VAL),
			min_z(HUGE_VAL), max_z(-HUGE_VAL) {}

	void add_ligand_box(OpenBabel::OBMol& mol)
	{
		using namespace OpenBabel;

		FOR_ATOMS_OF_MOL(a, mol)
		{
			min_x = std::min(min_x, (fl) a->x());
			min_y = std::min(min_y, (fl) a->y());
			min_z = std::min(min_z, (fl) a->z());
			max_x = std::max(max_x, (fl) a->x());
			max_y = std::max(max_y, (fl) a->y());
			max_z = std::max(max_z, (fl) a->z());
		}
	}

	//grow box in all dimensions by d
	void expand(fl d)
	{
		min_x -= d;
		max_x += d;
		min_y -= d;
		max_y += d;
		min_z -= d;
		max_z += d;
	}

	//return true if x,y,z is within box
	bool ptIn(fl x, fl y, fl z) const
	{
		return (x >= min_x && x <= max_x) &&
				(y >= min_y && y <= max_y) &&
				(z >= min_z && z <= max_z);
	}
};

#endif /* SMINA_BOX_H */
