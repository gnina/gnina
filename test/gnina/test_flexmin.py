#!/usr/bin/env python3

import sys
import subprocess

gnina = sys.argv[1]  # path to gnina executable

subprocess.check_call('%s  -r data/1w4o_rec.pdb -l data/1w4o_lig.mol2 \
    --cnn_scoring --minimize \
    --flexdist 3 --flexdist_ligand data/1w4o_lig.mol2 \
    -o lig-min.pdb --out_flex flex-min.pdb \
    --gpu'
    %(gnina), shell=True)

subprocess.check_call('%s  -r data/10gs_rec.pdb -l data/10gs_lig.mol2 \
    --cnn_scoring --minimize \
    --flexdist 3 --flexdist_ligand data/10gs_lig.mol2 \
    -o lig-min.pdb --out_flex flex-min.pdb \
    --gpu'
    %(gnina), shell=True)