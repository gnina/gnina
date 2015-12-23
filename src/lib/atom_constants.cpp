#include "atom_constants.h"

//auto-initailize smina atom types to the default with a global constructor
struct atomconstants_initializer
{
public:
	atomconstants_initializer()
	{
		for(size_t i = 0u; i < smina_atom_type::NumTypes; ++i)
			smina_atom_type::data[i] = smina_atom_type::default_data[i];
	}
};

atomconstants_initializer aci;
