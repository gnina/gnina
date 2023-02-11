#!/usr/bin/env python3

'''Regression and other general tests of gnina'''

from openbabel import pybel
import sys, os, subprocess
import numpy as np
import pytest
import re


gnina = sys.argv[1]  # take path to gnina executable as only argument

#single atoms
outfile = 'test.sdf'
def rmout():
    try:
        os.remove(outfile)
    except OSError:
        pass
        
def are_similar(xyz, sdf):
    mola = next(pybel.readfile('xyz',xyz))
    molb = next(pybel.readfile('sdf',sdf))
    #do an n^s comparison, ensure bijection
    atommap = dict()
    bseen = set()
    for a in mola.atoms:
        for b in molb.atoms:            
            dist = np.linalg.norm(np.array(a.coords)-np.array(b.coords))
            if dist < 0.1:
                assert a.idx not in atommap
                assert b.idx not in bseen
                atommap[a.idx] = b.idx
                bseen.add(b.idx)
                break
        else: #did not break, nothing matched a
            return False
    return True


def get_scores(output):
    output = output.decode();
    aff = None
    cnn = None
    m = re.search('Affinity: (\S+)',output)
    if m: aff = float(m.group(1))
    m = re.search('CNNaffinity: (\S+)',output)
    if m: cnn = float(m.group(1))
    return aff,cnn


rmout()
output = subprocess.check_output('%s  -r data/noelem_rec.pdb -l data/noelem.sdf --score_only '%(gnina),shell=True)
aff,cnn = get_scores(output)
assert aff < -8
assert cnn > 5
rmout()

output = subprocess.check_output('%s  -r data/noelem_rec.pdb -l data/noelem.sdf --score_only --scoring vinardo'%(gnina),shell=True)
aff,cnn = get_scores(output)
assert aff < -8
assert cnn > 5
rmout()
