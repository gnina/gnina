#!/bin/env python3

'''Given list of model file names, write out cpp that puts them in a dictionary.'''

import sys
from pathlib import Path

models = sys.argv[1] # ; separate list
outname = sys.argv[2]

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
    names.append(name)
    out.write(f'''extern char _binary_lib_models_{name}_pt_start[];
extern char _binary_lib_models_{name}_pt_end[];

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
  sort(names.begin(), names.end());
  return boost::algorithm::join(names, " ");
}
''')

