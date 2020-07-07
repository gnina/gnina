#!/usr/bin/env python3

'''Evaluate CNN minimization using simple overlay model'''

from openbabel import pybel
import sys, os, subprocess
import numpy as np
import pytest


gnina = sys.argv[1]  # take path to gnina executable as only argument

#single atoms
outfile = '/tmp/test.sdf'
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
        
         
rmout()
subprocess.check_call('%s  -r data/C.xyz -l data/C1.xyz --cnn_scoring=refinement --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s'%(gnina,outfile),shell=True)
#output should be similar to receptor
assert are_similar('data/C.xyz',outfile)

#gpu
rmout()
subprocess.check_call('%s  -r data/C.xyz -l data/C1.xyz --cnn_scoring=refinement --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s --gpu'%(gnina,outfile),shell=True)
#output should be similar to receptor
assert are_similar('data/C.xyz',outfile)
    
rmout()
subprocess.check_call('%s  -r data/CC.xyz -l data/CC2.xyz --cnn_scoring=all --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/CC.xyz',outfile)

#gpu
rmout()
subprocess.check_call('%s  -r data/CC.xyz -l data/CC2.xyz --cnn_scoring=all --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s --gpu'%(gnina,outfile),shell=True)
assert are_similar('data/CC.xyz',outfile)
    
    
#fully flexible - need smaller radius to avoid local minima
rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_update_min_frame --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_freeze_receptor --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

#gpu
rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --minimize -o %s --gpu'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_update_min_frame --minimize -o %s --gpu'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_freeze_receptor --minimize -o %s --gpu'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)
