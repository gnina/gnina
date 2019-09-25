#!/usr/bin/env python2

'''test ligand peturbation feature of molgrid layer
  Ligand peturbation moves the ligand relative to the receptor and outputs
  the _inverse_ transformation as x,y,z and angles (i.e., the transformation
  to apply to get the original pose back)
  
  optional bool peturb_ligand = 33 [default = false]; //apply additional random displacement/rotation to ligand, generate output of inverse transform
  optional float peturb_ligand_translate = 34 [default = 2.0]; //max amount to translate ligand by
  optional bool peturb_ligand_rotate = 35 [default = true]; //do random rotation
  
  We need to verify that the output transformation is correct and that
  peturbation does the right thing with random rotate/translate
  '''

import caffe
import numpy as np
import os,sys,re
from pyquaternion import Quaternion as qt
import molgrid

#these might change - are also assuming the first channel is aliphatic carbon
LIGCHANNELS = 19
RECCHANNELS = 16

import math

def test_perturb():

    def isclose(a, b, rel_tol=1e-06, abs_tol=1e-05):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
        
    def toEuler(q):
        '''wikipedia, I am forever in your debt - convert quaternion q to Euler angles'''
        (w,x,y,z) = q.elements
        sinr_cosp = +2.0 * (w*x + y*z);
        cosr_cosp = +1.0 - 2.0 * (x*x + y*y);
        roll = math.atan2(sinr_cosp, cosr_cosp)
    
        # pitch (y-axis rotation)
        sinp = +2.0 * (w*y - z*x);
        if math.fabs(sinp) >= 1:
            pitch = math.copysign(math.pi / 2, sinp) # use 90 degrees if out of range
        else:
            pitch = math.asin(sinp)
    
        #yaw (z-axis rotation)
        siny_cosp = +2.0 * (w*z + x*y);
        cosy_cosp = +1.0 - 2.0 * (y*y + z*z);  
        yaw = math.atan2(siny_cosp, cosy_cosp);
        
        return (roll,pitch,yaw)
        
    #create input file in subdir of current directory
    try:
        os.mkdir('peturbfiles')
    except OSError:
        pass
        
    
    coords = np.array([[-0.756,  0.   ,  0.   ],[0.756,  0.   ,  0.   ]])
    
    moltemplate='''2
    
    C          %f        %f        %f
    C          %f        %f        %f
    '''
    
    def setmol(coords,fname='CC.xyz'):    
        xyz = open('peturbfiles/'+fname,'w')
        xyz.write(moltemplate%tuple(coords.flatten()))
        xyz.close()
        
    setmol(coords)
    setmol(coords,'rec.xyz')
    
    
    types = '''1 rec.xyz CC.xyz
    '''
    
    typef = open('peturbfiles/input.types','w')
    typef.write(types)
    typef.close()


    model = '''layer {{
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      {topline}
      molgrid_data_param {{
        source: "peturbfiles/input.types"
        batch_size: 1
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: false
        root_folder: "peturbfiles"
        random_rotation: {rotate}
        random_translate: {translate}
        peturb_ligand: {peturb}
        peturb_ligand_rotate: {protate}
        peturb_ligand_translate: {ptranslate}
        cache_structs: false
        fix_center_to_origin: {fix} 
    
      }}
    }}'''
    
    
    modelf = open('peturbfiles/plain.model','w')
    modelf.write(model.format(rotate="false",translate=0.0,peturb='false',protate='false',ptranslate=0.0,topline='',fix='true'))
    modelf.close()
      
    modelf = open('peturbfiles/plainpeturb.model','w')
    modelf.write(model.format(rotate="false",translate=0.0,peturb='true',protate='true',ptranslate=4.0,topline='top: "peturb"',fix='true'))
    modelf.close()
        
    modelf = open('peturbfiles/rand.model','w')
    modelf.write(model.format(rotate="true",translate=5.0,peturb='false',protate='false',ptranslate=4.0,topline='',fix='true'))
    modelf.close()    
    
    modelf = open('peturbfiles/randpeturb.model','w')
    modelf.write(model.format(rotate="true",translate=5.0,peturb='true',protate='true',ptranslate=4.0,topline='top: "peturb"',fix='true'))
    modelf.close()   
    
    modelf = open('peturbfiles/fullrandpeturb.model','w')
    modelf.write(model.format(rotate="true",translate=5.0,peturb='true',protate='true',ptranslate=4.0,topline='top: "peturb"',fix='false'))
    modelf.close()   
    
    caffe.set_mode_gpu()
    
    #to check that a peturbation is correct, we generate a peturbed ligand grid,
    #then use the inverse of the output peturbation to transform the original ligand, use
    #the plain model to generate a grid from that and compare the generated densities
    
    center = coords.mean(axis=0)
    setmol(coords)
    net = caffe.Net('peturbfiles/plainpeturb.model',caffe.TEST)
    res = net.forward()
    peturb = res['peturb'][0].copy() #need to copy because returned values are references to blobs
    data = res['data'][0].copy()
    
    res = net.forward()
    peturb2 = res['peturb'][0].copy()
    data2 = res['data'][0].copy()
    
    if len(peturb) != 10:
        print("Wrong number of rotational numbers")
        sys.exit(1)
        
    trans = peturb[:3]
    q = qt(peturb[3:7])    
    euler = toEuler(q)
    if not isclose(np.square(euler-peturb[-3:]).sum(),0):
        print("Euler angles don't match quaternion")
        sys.exit(1)
                
    #first make sure different peturbations genreate different grids
    if np.square(data-data2).sum() < 0.25:
        print("Not generating distinct grids with peturbation")
        sys.exit(1)
        
    
    newcoords = np.array(list(map(q.inverse.rotate,coords-center)))+center-trans
    
    setmol(newcoords)
    plainnet = caffe.Net('peturbfiles/plain.model',caffe.TEST)
    plainres = plainnet.forward()
    plaindata = plainres['data'][0].copy()
    
    diff = np.square(plaindata - data).sum()
    if not isclose(diff,0):
        print(("Peturbed ligand not as expected %f"%diff))
        sys.exit(1)
    
    
    #now have to check the rotated case
    #first make sure random transformation is working as expected
    setmol(coords)
    setmol(coords,'rec.xyz')
    net = caffe.Net('peturbfiles/rand.model',caffe.TEST)
    res = net.forward()
    data = res['data'][0].copy()
    l =net.layers[0]
    mt = l.get_moltransform(0)
    trans = np.array(tuple(mt.get_translation())) 
    Q = mt.get_quaternion()
    q = qt(Q.R_component_1(),Q.R_component_2(),Q.R_component_3(),Q.R_component_4())
    #rotating
    newcoords = np.array(list(map(q.rotate,coords-center)))+center+trans
    
    setmol(newcoords)
    setmol(newcoords,'rec.xyz')
    plainnet = caffe.Net('peturbfiles/plain.model',caffe.TEST)
    plainres = plainnet.forward()
    plaindata = plainres['data'][0].copy()
    
    diff = np.square(plaindata - data).sum()
    if not isclose(diff,0):
        print(("Rotated ligand not as expected %f"%diff))
        sys.exit(1)
        
    #check peturbed ligand with rotation
    setmol(coords)
    setmol(coords,'rec.xyz')
    net = caffe.Net('peturbfiles/randpeturb.model',caffe.TEST)
    res = net.forward()
    data = res['data'][0].copy()
    peturb = res['peturb'][0].copy() 
    trans = peturb[:3]
    q = qt(peturb[3:7])  
    
    l =net.layers[0]
    mt = l.get_moltransform(0)
    Rtrans = np.array(tuple(mt.get_translation())) 
    Q = mt.get_quaternion()
    Rq = qt(Q.R_component_1(),Q.R_component_2(),Q.R_component_3(),Q.R_component_4())    
    
    #first randomize, then peturb
    center = np.mean(coords,axis=0)
    newreccoords = np.array(list(map(Rq.rotate,coords-center)))+center+Rtrans
    newligcoords = np.array(list(map(Rq.rotate,coords-center)))+center+Rtrans
    center = np.mean(newligcoords,axis=0)
    newligcoords = np.array(list(map(q.inverse.rotate,newligcoords-center)))+center-trans
    
    setmol(newligcoords)
    setmol(newreccoords,'rec.xyz')
    plainnet = caffe.Net('peturbfiles/plain.model',caffe.TEST)
    plainres = plainnet.forward()
    plaindata = plainres['data'][0].copy()
    
    diff = np.square(plaindata - data).sum()
    if not isclose(diff,0):
        print(("Rotated and peturbed ligand not as expected %f"%diff))
        sys.exit(1)
    
    #don't fix to the origin - result will be centered around ligand, and then translated
    #rotation center is grid center, not molecule center
    setmol(coords)
    setmol(coords,'rec.xyz')
    net = caffe.Net('peturbfiles/fullrandpeturb.model',caffe.TEST)
    res = net.forward()
    data = res['data'][0].copy()
    peturb = res['peturb'][0].copy() 
    trans = peturb[:3]
    q = qt(peturb[3:7])  
    
    l = net.layers[0]
    mt = l.get_moltransform(0)
    Rtrans = np.array(tuple(mt.get_translation())) 
    Q = mt.get_quaternion()
    Rq = qt(Q.R_component_1(),Q.R_component_2(),Q.R_component_3(),Q.R_component_4())    
    
    center = np.mean(coords,axis=0)
    
    #first rotate
    newligcoords = np.array(list(map(Rq.rotate,coords-center)))+center+Rtrans
    newreccoords = np.array(list(map(Rq.rotate,coords-center)))+center+Rtrans
    center = np.mean(newligcoords,axis=0)
    #then peturb, using updated ligand center
    newligcoords = np.array(list(map(q.inverse.rotate,newligcoords-center)))+center-trans
    
    setmol(newligcoords)
    setmol(newreccoords,'rec.xyz')
    plainnet = caffe.Net('peturbfiles/plain.model',caffe.TEST)
    plainres = plainnet.forward()
    plaindata = plainres['data'][0].copy()
    
    diff = np.square(plaindata - data).sum()
    if not isclose(diff,0):
        print(("Fully Rotated and peturbed ligand not as expected %f"%diff))
        sys.exit(1)

