#!/bin/env python3

'''Given list of model file names, write out cpp that puts them in a dictionary.'''

import sys
from pathlib import Path

models = sys.argv[1] # ; separate list
outname = sys.argv[2]
# On Linux, GNU ld's `-r -b binary` makes _end a true linker section-boundary
# symbol, so it must be declared as an array (its "value" is its address).
# On macOS, embed_model.py has no such mechanism and instead defines _end as
# a real pointer variable (start + size), so it must be declared as a pointer.
# Mixing these up compiles fine but silently corrupts every embedded model's
# computed size, since the two symbol kinds are read completely differently.
is_apple = len(sys.argv) > 3 and sys.argv[3] == 'apple'

out = open(outname,'wt')

out.write('''/*
 * GENERATED - DO NOT EDIT
 */

#include <boost/assign.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <vector>
#include <boost/unordered_map.hpp>

''')

names = []
for model in models.split(';'):
    name = Path(model).stem
    name = name.replace('.','_')
    names.append(name)
    end_decl = f'extern char* _binary_lib_models_{name}_pt_end;' if is_apple \
        else f'extern char _binary_lib_models_{name}_pt_end[];'
    out.write(f'''extern char _binary_lib_models_{name}_pt_start[];
{end_decl}

''')

out.write('\n\nboost::unordered_map<std::string, std::pair<char*, char*> > torch_models = boost::assign::map_list_of')

for name in names:
    out.write(f'("{name}",std::make_pair(_binary_lib_models_{name}_pt_start,_binary_lib_models_{name}_pt_end))\n')
out.write(';\n\n')

out.write('''
std::string builtin_torch_models()
{
  std::vector<std::string> names;
  names.reserve(torch_models.size());
  for (auto kv : torch_models) {
    names.push_back(kv.first);
  }
  names.push_back("fast");
  names.push_back("default1.0"); 

  sort(names.begin(), names.end());
  return boost::algorithm::join(names, " ");
}
''')

