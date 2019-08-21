#!/usr/bin/env python3

import sys, os
import subprocess

import numpy as np

from openbabel import pybel

gnina = sys.argv[1]  # path to gnina executable

def rmout(*files):
    for file in files:
        try:
            os.remove(file)
        except OSError:
            pass

def moved(output, reference, refformat):
    mola = next(pybel.readfile('pdb', output))
    molb = next(pybel.readfile(refformat, reference))

    for a, b in zip(mola.atoms, molb.atoms):           
            dist = np.linalg.norm(np.array(a.coords)-np.array(b.coords))
            if dist > 1e-3:
                return True

    print(output)

    return False

outlig="lig-min.pdb"
outflex="flex-min.pdb"

subprocess.check_call("%s  -r data/10gs_rec.pdb -l data/10gs_lig.mol2 \
    --cnn_scoring --minimize \
    --flexdist 3 --flexdist_ligand data/10gs_lig.mol2 \
    -o %s --out_flex %s \
    --gpu"
    %(gnina, outlig, outflex), shell=True)

assert moved(outlig, "data/10gs_lig.mol2", "mol2")
assert moved(outflex, "data/10gs_rec_ref.pdb", "pdb")
rmout(outlig, outflex)