#include "box.h"

#include "obmolopener.h"
#include "model.h"

//generate a box around the provided ligand padded by autobox_add
//initial model is provided to include flexible residues
//if centers are always overwritten, but if sizes are non zero they are preserved
void setup_autobox(const model& m, const std::string& autobox_ligand, fl autobox_add,
    fl& center_x, fl& center_y, fl& center_z, fl& size_x, fl& size_y,
    fl& size_z) {
  using namespace OpenBabel;
  obmol_opener opener;

  //a ligand file can be provided from which to autobox
  OBConversion conv;
  opener.openForInput(conv, autobox_ligand);

  OBMol mol;
  center_x = center_y = center_z = 0;
  Box b;
  while (conv.Read(&mol)) //openbabel divides separate residues into multiple molecules
  {
    b.add_ligand_box(mol);
  }

  //flexible residues
  for(const auto& a : m.get_movable_atoms()) {
    b.add_coord(a.coords[0], a.coords[1], a.coords[2]);
  }

  //set to center of bounding box (as opposed to center of mass
  center_x = (b.max_x + b.min_x) / 2.0;
  center_y = (b.max_y + b.min_y) / 2.0;
  center_z = (b.max_z + b.min_z) / 2.0;

  b.expand(autobox_add);
  if (size_x == 0) size_x = (b.max_x - b.min_x);
  if (size_y == 0) size_y = (b.max_y - b.min_y);
  if (size_z == 0) size_z = (b.max_z - b.min_z);

  if (!std::isfinite(b.max_x) || !std::isfinite(b.max_y)
      || !std::isfinite(b.max_z)) {
    std::stringstream msg;
    msg << "Unable to read  " << autobox_ligand;
    throw usage_error(msg.str());
  }
}
