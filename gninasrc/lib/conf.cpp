#include "conf.h"

bool change::operator==(const change& other) const {
  sz n = num_floats();
  assert(other.num_floats() == n);
  for (sz i = 0; i < num_floats(); i++) {
    sz ignore;
    if (get_with_node_idx(i, &ignore, &ignore)
        != other.get_with_node_idx(i, &ignore, &ignore)) return false;
  }
  return true;
}

fl change::get_with_node_idx(sz index, sz* node_idx, sz* offset_in_node) const {
  *node_idx = 0;
  *offset_in_node = 0;
  VINA_FOR_IN(i, ligands) {
    const ligand_change& lig = ligands[i];
    if (index < 6) *offset_in_node = index;
    if (index < 3) return lig.rigid.position[index];
    index -= 3;
    if (index < 3) return lig.rigid.orientation[index];
    index -= 3;

    // Count the number of nodes we have/will skip over on this
    // iteration.
    *node_idx += 1 + std::min(index, lig.torsions.size());

    if (index < lig.torsions.size()) return lig.torsions[index];
    index -= lig.torsions.size();
  }

  VINA_FOR_IN(i, flex) {
    const residue_change& res = flex[i];

    *node_idx += std::min(index, res.torsions.size());
    if (index < res.torsions.size()) return res.torsions[index];
    index -= res.torsions.size();
  }

  if (index < 3)
    return receptor.position[index];
  else
    if (index < 6)
      return receptor.orientation[index - 3];
    else
      abort();
}

bool conf::operator==(const conf& other) const {
  sz n = num_floats();
  assert(other.num_floats() == n);
  for (sz i = 0; i < num_floats(); i++) {
    sz ignore;
    if (get_with_node_idx(i, &ignore, &ignore)
        != other.get_with_node_idx(i, &ignore, &ignore)) return false;
  }
  return true;
}

fl conf::get_with_node_idx(sz index, sz* node_idx, sz* offset_in_node) const {
  *node_idx = 0;
  *offset_in_node = 0;
  VINA_FOR_IN(i, ligands) {
    const ligand_conf& lig = ligands[i];
    if (index < 7) *offset_in_node = index;
    if (index < 3) return lig.rigid.position[index];
    index -= 3;
    // TODO CPU version has this
    // if(index < 3)
    // {
    //     vec ang = quaternion_to_angle(lig.rigid.orientation);
    //     return ang[index];
    // }
    // TODO instead of this:
    if (index < 4) return ((float*) &lig.rigid.orientation)[index];
    index -= 4;
    // Count the number of nodes we have/will skip over on this
    // iteration.
    *node_idx += 1 + std::min(index, lig.torsions.size());

    if (index < lig.torsions.size()) return lig.torsions[index];
    index -= lig.torsions.size();
  }

  VINA_FOR_IN(i, flex) {
    const residue_conf& res = flex[i];

    *node_idx += std::min(index, res.torsions.size());
    if (index < res.torsions.size()) return res.torsions[index];
    index -= res.torsions.size();
  }
  __builtin_unreachable();
}

fl& conf::flat_index(sz index) { // returns by value
  VINA_FOR_IN(i, ligands) {
    ligand_conf& lig = ligands[i];
    if (index < 3) return lig.rigid.position[index];
    index -= 3;
    //TODO
    // if(index < 3)
    // {
    //     vec ang = quaternion_to_angle(lig.rigid.orientation);
    //     return ang[index];
    // }
    if (index < 4) return ((float*) &lig.rigid.orientation)[index];

    index -= 4;
    if (index < lig.torsions.size()) return lig.torsions[index];
    index -= lig.torsions.size();
  }
  VINA_FOR_IN(i, flex) {
    residue_conf& res = flex[i];
    if (index < res.torsions.size()) return res.torsions[index];
    index -= res.torsions.size();
  }
  //VINA_CHECK(false); //avoid compiler warnings
  __builtin_unreachable();
}
