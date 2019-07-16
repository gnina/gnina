#!/usr/bin/env python3

import sys
import subprocess

gnina = sys.argv[1]  # path to gnina executable

subprocess.check_call('%s  -r data/rec.pdb -l data/lig.mol2 \
    --cnn_scoring --minimize \
    --flexdist 3 --flexdist_ligand data/lig.mol2 \
    -o lig-min.pdb --out_flex flex.pdb \
    --gpu'
    %(gnina), shell=True)