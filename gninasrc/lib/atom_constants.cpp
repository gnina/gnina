#include "atom_constants.h"

//auto-initailize smina atom types to the default with a global constructor

namespace smina_atom_type {
info data[NumTypes] = { { }, };

struct atom_data_initializer {
    atom_data_initializer() {
      for (size_t i = 0u; i < smina_atom_type::NumTypes; ++i)
        smina_atom_type::data[i] = smina_atom_type::default_data[i];
    }
};

atom_data_initializer initialize_defaults;
}
