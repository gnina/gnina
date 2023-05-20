#!/usr/bin/env python3

'''Evaluate CNN minimization using simple overlay model'''

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

def get_CNNscore(outfile):
    mol = next(pybel.readfile('sdf',outfile))
    return float(mol.data['CNNscore'])

def get_energies(output):
    for line in output.decode().split('\n'):
        if 'Total energy after refinement' in line:
            total_e=float(line.split()[-1])
        elif 'Empirical energy after refinement' in line:
            emp_e=float(line.split()[-1])
    return total_e,emp_e

def validate_energies(outfile,output,weight):
    CNNscore=get_CNNscore(outfile)
    total_e,emp_e=get_energies(output)
    # (CNN loss + weight * Vina)/( 1 + weight)
    calc_e= ( -np.log(CNNscore) + weight * emp_e) / ( 1 + weight )
    if np.abs(total_e - calc_e) < 0.001:
        return True
    else:
        return False

rmout()
subprocess.check_call('%s  -r data/C.xyz -l data/C1.xyz --cnn_scoring=refinement --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --no_gpu --minimize -o %s'%(gnina,outfile),shell=True)
#output should be similar to receptor
assert are_similar('data/C.xyz',outfile)

#gpu
rmout()
subprocess.check_call('%s  -r data/C.xyz -l data/C1.xyz --cnn_scoring=refinement --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s '%(gnina,outfile),shell=True)
#output should be similar to receptor
assert are_similar('data/C.xyz',outfile)
    
rmout()
subprocess.check_call('%s  -r data/CC.xyz -l data/CC2.xyz --cnn_scoring=all --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s --no_gpu'%(gnina,outfile),shell=True)
assert are_similar('data/CC.xyz',outfile)

#gpu
rmout()
subprocess.check_call('%s  -r data/CC.xyz -l data/CC2.xyz --cnn_scoring=all --cnn_model data/overlap.model \
    --cnn_weights data/overlay.caffemodel --minimize -o %s '%(gnina,outfile),shell=True)
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
  --cnn_update_min_frame=true --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)
assert re.search(r'<molecular weight>',open(outfile).read()) #should retain sddata tag

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_freeze_receptor --minimize -o %s --no_gpu'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

#gpu
rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_update_min_frame=true --minimize -o %s'%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

rmout()
subprocess.check_call('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_freeze_receptor --minimize -o %s '%(gnina,outfile),shell=True)
assert are_similar('data/C8flat.xyz',outfile)

#Vina + CNN

rmout()
output=subprocess.check_output('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_mix_emp_force --cnn_mix_emp_energy --cnn_empirical_weight 1 --verbosity 2 --minimize --no_gpu -o %s'%(gnina,outfile),shell=True)
#output should not be similar to receptor
assert not are_similar('data/C8flat.xyz',outfile)
assert validate_energies(outfile,output,1)

rmout()
output=subprocess.check_output('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_mix_emp_force --cnn_mix_emp_energy --cnn_empirical_weight 1 --verbosity 2 --minimize -o %s'%(gnina,outfile),shell=True)
#output should not be similar to receptor
assert not are_similar('data/C8flat.xyz',outfile)
assert validate_energies(outfile,output,1)

rmout()
output=subprocess.check_output('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_mix_emp_force --cnn_mix_emp_energy --cnn_empirical_weight 0.7 --verbosity 2 --minimize --no_gpu -o %s'%(gnina,outfile),shell=True)
#output should not be similar to receptor
assert not are_similar('data/C8flat.xyz',outfile)
assert validate_energies(outfile,output,0.7)

rmout()
output=subprocess.check_output('%s  -r data/C8flat.xyz -l data/C8bent.sdf --cnn_scoring=all \
  --cnn_model data/overlap_smallr.model --cnn_weights data/overlay.caffemodel \
  --cnn_mix_emp_force --cnn_mix_emp_energy --cnn_empirical_weight 0.7 --verbosity 2 --minimize -o %s'%(gnina,outfile),shell=True)
#output should not be similar to receptor
assert not are_similar('data/C8flat.xyz',outfile)
assert validate_energies(outfile,output,0.7)

rmout()
