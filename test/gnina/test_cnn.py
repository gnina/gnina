#!/usr/bin/env python3

'''Evaluate CNN scoring'''

import sys, os, subprocess, re
import numpy as np
import pytest


gnina = sys.argv[1]  # take path to gnina executable as only argument

def getscores(out):
    '''Read scores for docked poses and return list of tuples'''
    scores = re.findall(r'^\s*\d+\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)', out.decode(), re.MULTILINE)
    return [(float(a),float(b),float(c)) for (a,b,c) in scores]

def issorted(scores, index):
    mult = 1
    if index > 0:
        mult = -1;
    return all(mult*scores[i][index] <= mult*scores[i+1][index] for i in range(len(scores)-1))


#need to set num_modes=10 since the default is to generate 10 and throw out one afte rsorting
defaultout = subprocess.check_output('%s  --no_gpu -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10'%gnina,shell=True)
#should be sorted by CNNscore
defaultout = getscores(defaultout)
assert issorted(defaultout,1)

energyout = subprocess.check_output('%s --no_gpu -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --pose_sort_order=energy'%gnina,shell=True)
energyout = getscores(energyout)
assert issorted(energyout,0)

affout = subprocess.check_output('%s --no_gpu -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --pose_sort_order=CNNaffinity'%gnina,shell=True)
affout = getscores(affout)
assert issorted(affout,2)

defaultoutg = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10'%gnina,shell=True)
#should be sorted by CNNscore
defaultoutg = getscores(defaultoutg)
assert issorted(defaultoutg,1)

np.testing.assert_array_almost_equal(defaultout,defaultoutg,3)

energyoutg = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --pose_sort_order=energy'%gnina,shell=True)
energyoutg = getscores(energyoutg)
assert issorted(energyoutg,0)
np.testing.assert_array_almost_equal(energyout,energyoutg,3)


affoutg = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --pose_sort_order=CNNaffinity'%gnina,shell=True)
affoutg = getscores(affoutg)
assert issorted(affoutg,2)
np.testing.assert_array_almost_equal(affout,affoutg,3)

refine = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn_scoring=refinement'%gnina,shell=True)
#should be sorted by CNNscore
refine = getscores(refine)
assert refine[0][1] > defaultout[0][1]

singleout = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn=general_default2018'%gnina,shell=True)
#should be sorted by CNNscore
singleout = getscores(singleout)
assert issorted(singleout,1)

#make sure lack of gpu still works
nogpuout = subprocess.check_output('CUDA_VISIBLE_DEVICES= %s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn=general_default2018'%gnina,shell=True)
#should be sorted by CNNscore
nogpuout = getscores(nogpuout)
assert issorted(nogpuout,1)
np.testing.assert_array_almost_equal(nogpuout,singleout,decimal=2)


ensembleout = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn=general_default2018_ensemble -o testingout.sdf'%gnina,shell=True)
#should be sorted by CNNscore
ensembleout = getscores(ensembleout)
assert issorted(ensembleout,1)

assert not np.array_equal(ensembleout,singleout)

#fewer steps for metropolis check
ensembleout = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn=general_default2018_ensemble --num_mc_steps 1000 '%gnina,shell=True)
#should be sorted by CNNscore
ensembleout = getscores(ensembleout)
assert issorted(ensembleout,1)

#metropolis sampling
metropolisrescoreout = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn=general_default2018_ensemble --num_mc_steps 1000 --cnn_scoring metrorescore '%gnina,shell=True)
#should be sorted by CNNscore
metropolisrescoreout = getscores(metropolisrescoreout)
assert issorted(metropolisrescoreout,1)

#metropolis sampling
metropolisrefineout = subprocess.check_output('%s  -r data/184l_rec.pdb -l data/184l_lig.sdf --autobox_ligand data/184l_lig.sdf --seed 2 --num_modes=10 --cnn=general_default2018_ensemble --num_mc_steps 1000 --cnn_scoring metrorefine '%gnina,shell=True)
#should be sorted by CNNscore
metropolisrefineout = getscores(metropolisrefineout)
assert issorted(metropolisrefineout,1)

assert not np.array_equal(metropolisrefineout,metropolisrescoreout)
assert not np.array_equal(ensembleout,metropolisrescoreout)

assert 'CNNaffinity_variance' in open('testingout.sdf').read()

os.remove('testingout.sdf')
