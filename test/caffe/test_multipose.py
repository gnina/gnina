#!/usr/bin/env python2

'''Test multiple pose mol gridder by checking outputs'''

import caffe
import numpy as np
import os,sys,re

#these might change
LIGCHANNELS = 14
RECCHANNELS = 14

def test_multipose():

    #create input file in subdir of current directory
    try:
        os.mkdir('multifiles')
    except OSError:
        pass
        
    
    xyz = open('multifiles/C.xyz','w')
    xyz.write('''1
    
    C          0.00000        0.00000        0.00000
    ''')
    xyz.close()
    
    xyz = open('multifiles/N.xyz','w')
    xyz.write('''1
    
    N          0.00000        0.00000        0.00000
    ''')
    xyz.close()
    
    xyz = open('multifiles/O.xyz','w')
    xyz.write('''1
    
    O          0.00000        0.00000        0.00000
    ''')
    xyz.close()
    
    xyz = open('multifiles/CC.xyz','w')
    xyz.write('''8
    
    C          1.06088        0.05280        0.06303
    C          2.57294        0.05281        0.06302
    H          0.67651       -0.71498       -0.61505
    H          0.67651        1.02393       -0.26286
    H          0.67651       -0.15053        1.06698
    H          2.95731        0.25614       -0.94093
    H          2.95731       -0.91832        0.38891
    H          2.95731        0.82059        0.74110
    ''')
    xyz.close()
    
    model = '''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      molgrid_data_param {
        source: "multifiles/input.types"
        batch_size: 2
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: false
        root_folder: "multifiles"
        random_rotation: true
        random_translate: 6.0
        num_poses: 3
        duplicate_poses: DUPLICATE
      }
    }'''
    
    modelf = open('multifiles/duplicate.model','w')
    modelf.write(model.replace('DUPLICATE','true'))
    modelf.close()
    
    modelf = open('multifiles/noduplicate.model','w')
    modelf.write(model.replace('DUPLICATE','false'))
    modelf.close()
    
    types = '''1 CC.xyz C.xyz N.xyz O.xyz
    0 CC.xyz CC.xyz CC.xyz CC.xyz
    '''
    
    typef = open('multifiles/input.types','w')
    typef.write(types)
    typef.close()
    
    net = caffe.Net('multifiles/duplicate.model',caffe.TRAIN)
    res = net.forward()
    
    if not res['label'][:3].all() or not np.all(res['label'][3:] == 0):
        print("Label not correct:",res['label'])
        sys.exit(1)
        
    data = res['data']
    if data.shape[0] != 6 or data.shape[1] != (LIGCHANNELS+RECCHANNELS):
        print("Data incorrect shape:",data.shape)
        sys.exit(1)
    
    #every batch, ligand and receptor, should be nonzero
    for i in range(6):
        if np.sum(data[i]) == 0:
            print("Zero batch",i)
            sys.exit(1)
        if np.sum(data[i][:RECCHANNELS]) == 0:
            print("Zero receptor",i)
            sys.exit(1)
        if np.sum(data[i][RECCHANNELS:]) == 0:
            print("Zero ligand",i)
            sys.exit(1)
            
    #in the first example, each ligand should populate a different atom channel
    seenligch = set()
    for i in range(3):
        ch = np.apply_over_axes(np.sum,data[i][RECCHANNELS:],[1,2,3]).flatten().argmax()
        seenligch.add(ch)
    if len(seenligch) != 3:
        print("Incorrect lig channels")
        sys.exit(1)
        
    #the receptor total density should be roughly the same
    recdensity = np.apply_over_axes(np.sum,data[:,:RECCHANNELS],[1,2,3,4]).flatten()
    averec = np.mean(recdensity)
    recdiff = np.abs(recdensity-averec)
    if recdiff.max() > 2.0:
        print("receptor density too different",recdiff.max())
        sys.exit(1)
    
    #the second example, the receptor and ligand should be rotated the same
    for i in range(3,6):
        diff = np.sum((data[i][:4] - data[i][RECCHANNELS:RECCHANNELS+4])**2) #4 is carbon offset
        if diff > 0.01:
            print("ligand/receptor transformation mismatch",i)
            sys.exit(1)
    
    
    #for nonduplicate, should be a single receptor and multiple ligands
    net = caffe.Net('multifiles/noduplicate.model',caffe.TRAIN)
    res = net.forward()
    data = res['data']
    labels = res['label']
    
    if len(labels) != 2 or labels[0] != 1.0 or labels[1] != 0:
        print("Incorrect labels",labels)
        sys.exit(1)
        
    #check size
    if data.shape[0] != 2:
        print("Incorrect batchsize",data.shape)
        sys.exit(1)
        
    if data.shape[1] != RECCHANNELS+LIGCHANNELS*3:
        print("Incorrect channel number",data.shape)
        sys.exit(1)
    
    #check receptor
    for i in range(2):
        rd = data[i][:RECCHANNELS].sum()
        if np.abs(rd-averec) > 2:
            print("Incorrect receptor density",rd,averec,i)
            sys.exit(1)
            
    #check channels
    for i in range(3):
        start = RECCHANNELS+LIGCHANNELS*i
        end = start+LIGCHANNELS
        if data[0][start:end].sum() == 0:
            print("Zero dup ligand density",i)
            sys.exit(1)
        #second example should be identical - same rotation applied
        diff = np.abs(data[1][start:end]-data[1][RECCHANNELS:RECCHANNELS+LIGCHANNELS]).sum()
        if diff != 0:
            print("Dup not identical",diff,i)
            sys.exit(1)
            
