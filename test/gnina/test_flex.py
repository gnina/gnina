#!/usr/bin/env python3

# Test docking/minimsation with flexible residues

# Reference selections are checked against the following selection in PyMol
# (after stripping away all hydrogen atoms):
#
# select byres 
#   (RECEPTOR and not backbone 
#       and not resname ALA and not resname GLY and not resname PRO
#   ) within DISTANCE of LIGAND

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
    """
    Check if the output moved.
    If the reference is a PDB file, check that the residue names correspond.
    """
    mola = next(pybel.readfile('pdb', output))
    molb = next(pybel.readfile(refformat, reference))

    if refformat == "pdb":
        assert len(mola.residues) == len(molb.residues)

        for resa, resb in zip(mola.residues, molb.residues):
            assert resa.name.strip() == resb.name.strip()
            assert resa.OBResidue.GetChain() == resb.OBResidue.GetChain()

    for a, b in zip(mola.atoms, molb.atoms):
            dist = np.linalg.norm(np.array(a.coords)-np.array(b.coords))
            if dist > 1e-3:
                return True

    print(output)

    return False

# Test if side chains moved after minimisation with flexible residues

for system, flexdist in [("10gs", 3.0), ("184l", 3.9)]:

    outlig="lig-min-{system}.pdb".format(system=system)
    outflex="flex-min-{system}.pdb".format(system=system)

    subprocess.check_call("{gnina} -r data/{system}_rec.pdb -l data/{system}_lig.sdf \
        --cnn_scoring=all --minimize \
        --flexdist {flexdist} --flexdist_ligand data/{system}_lig.sdf \
        -o {outlig} --out_flex {outflex} \
        --gpu".format(gnina=gnina, outlig=outlig, outflex=outflex, system=system, flexdist=flexdist),
        shell=True)

    assert moved(outlig, "data/{system}_lig.sdf".format(system=system), "sdf")
    assert moved(outflex, "data/{system}_rec_ref.pdb".format(system=system), "pdb")
    rmout(outlig, outflex)

# Test if side chains moved after docking with flexible residues

for system, flexdist in [("10gs", 3.0), ("184l", 3.9)]:

    outlig="lig-dock-{system}.pdb".format(system=system)
    outflex="flex-dock-{system}.pdb".format(system=system)

    subprocess.check_call("{gnina} -r data/{system}_rec.pdb -l data/{system}_lig.sdf \
        --autobox_ligand data/{system}_lig.sdf --autobox_add 6 \
        --flexdist {flexdist} --flexdist_ligand data/{system}_lig.sdf \
        --num_mc_steps 50 \
        -o {outlig} --out_flex {outflex}".format(gnina=gnina, outlig=outlig, outflex=outflex, system=system, flexdist=flexdist), 
        shell=True)

    assert moved(outlig, "data/{system}_lig.sdf".format(system=system), "sdf")
    assert moved(outflex, "data/{system}_rec_ref.pdb".format(system=system), "pdb")
    rmout(outlig, outflex)

# Test --flex_max selection of closest residues
# --flexdist 3.9 selects 4 residues for 184l
for closest in [1, 2, 3]:

    outlig="lig-dock-184l-c{closest}.pdb".format(closest=closest)
    outflex="flex-dock-184l-c{closest}.pdb".format(closest=closest)

    subprocess.check_call("{gnina} -r data/184l_rec.pdb -l data/184l_lig.sdf \
        --autobox_ligand data/184l_lig.sdf --autobox_add 6 \
        --flexdist 3.9 --flexdist_ligand data/184l_lig.sdf \
        --flex_max {closest} \
        --num_mc_steps 50 \
        -o {outlig} --out_flex {outflex}".format(gnina=gnina,
            outlig=outlig,
            outflex=outflex,
            closest=closest),
        shell=True)

    assert moved(outlig, "data/184l_lig.sdf", "sdf")
    assert moved(outflex, "data/184l_rec_ref_c{closest}.pdb".format(closest=closest), "pdb")
    rmout(outlig, outflex)
