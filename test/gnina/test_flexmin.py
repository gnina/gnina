#!/usr/bin/env python3

import sys, os
import subprocess

gnina = sys.argv[1]  # path to gnina executable

def rmout(*files):
    for file in files:
        try:
            os.remove(file)
        except OSError:
            pass

outlig="lig-min.pdb"
outflex="flex-min.pdb"

subprocess.check_call("%s  -r data/1w4o_rec.pdb -l data/1w4o_lig.mol2 \
    --cnn_scoring --minimize \
    --flexdist 3 --flexdist_ligand data/1w4o_lig.mol2 \
    -o %s --out_flex %s \
    --gpu"
    %(gnina, outlig, outflex), shell=True)

rmout(outlig, outflex)

subprocess.check_call("%s  -r data/10gs_rec.pdb -l data/10gs_lig.mol2 \
    --cnn_scoring --minimize \
    --flexdist 3 --flexdist_ligand data/10gs_lig.mol2 \
    -o %s --out_flex %s \
    --gpu"
    %(gnina, outlig, outflex), shell=True)

rmout(outlig, outflex)