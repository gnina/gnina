#include "conf.h"

fl change::get_with_node_idx(sz index, sz* node_idx) const {
    VINA_FOR_IN(i, ligands) {
        *node_idx = i;

        const ligand_change& lig = ligands[i];
        if(index < 3) return lig.rigid.position[index];
        index -= 3;
        if(index < 3) return lig.rigid.orientation[index];
        index -= 3;
        if(index < lig.torsions.size()) return lig.torsions[index];
        index -= lig.torsions.size();
    }

    VINA_FOR_IN(i, flex) {
        *node_idx = i + ligands.size();

        const residue_change& res = flex[i];
        if(index < res.torsions.size()) return res.torsions[index];
        index -= res.torsions.size();
    }

    __builtin_unreachable();
}

fl conf::get_with_node_idx(sz index, sz* node_idx) const {
    VINA_FOR_IN(i, ligands) {
        *node_idx = i;

        const ligand_conf& lig = ligands[i];
        if(index < 3) return lig.rigid.position[index];
        index -= 3;
        // TODO CPU version has this
        // if(index < 3)
        // {
        //     vec ang = quaternion_to_angle(lig.rigid.orientation);
        //     return ang[index];
        // }
        // TODO instead of this:
        if(index < 4) return ((float*)&lig.rigid.orientation)[index];

        index -= 4;
        if(index < lig.torsions.size()) return lig.torsions[index];
        index -= lig.torsions.size();
    }

    VINA_FOR_IN(i, flex) {
        *node_idx = i + ligands.size();

        const residue_conf& res = flex[i];
        if(index < res.torsions.size()) return res.torsions[index];
        index -= res.torsions.size();
    }
    __builtin_unreachable();
}

fl& conf::operator()(sz index) { // returns by value
    VINA_FOR_IN(i, ligands) {
        ligand_conf& lig = ligands[i];
        if(index < 3) return lig.rigid.position[index];
        index -= 3;
        //TODO
        // if(index < 3)
        // {
        //     vec ang = quaternion_to_angle(lig.rigid.orientation);
        //     return ang[index];
        // }
        if(index < 4) return ((float*)&lig.rigid.orientation)[index];

        index -= 4;
        if(index < lig.torsions.size()) return lig.torsions[index];
        index -= lig.torsions.size();
    }
    VINA_FOR_IN(i, flex) {
        residue_conf& res = flex[i];
        if(index < res.torsions.size()) return res.torsions[index];
        index -= res.torsions.size();
    }
    //VINA_CHECK(false); //avoid compiler warnings
    __builtin_unreachable();
}
