#!/usr/bin/env python3

'''Test molgrid grouping'''

import caffe
import numpy as np
import os,sys,re

import pytest
from pytest import approx

def test_groups():

    caffe.set_random_seed(800)
    caffe.set_mode_gpu()
   
    
    def runmodel(model):
        m = open('tmp.model','w')
        m.write(model)
        m.close()
        net = caffe.Net('tmp.model',caffe.TRAIN)
        res = net.forward()
        os.remove('tmp.model')
        return res,net
    
    #plain
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      top: "seqcont"
      molgrid_data_param {
        source: "typesfiles/grouped.types"
        batch_size: 2
        max_group_size: 4
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: false
        root_folder: "typesfiles"
      }
    }''')
    labels = res['label']
    seqcont = res['seqcont']
    data = res['data']
    assert labels.shape == (4,2)
    assert seqcont.shape == (4,2)
    assert list(seqcont[0]) == [0,0]
    assert list(seqcont[1]) == [1,1]
    assert data[3][0].sum() == 0
    assert data[3][1].sum() == 0
    
    res = net.forward()
    labels = res['label']
    seqcont = res['seqcont']
    data = res['data']
    assert labels.shape == (4,2)
    assert seqcont.shape == (4,2)
    assert list(seqcont[0]) == [0,0]
    assert list(seqcont[1]) == [1,1]
    assert data[3][0].sum() > 0
    assert data[3][1].sum() > 0
    
    #with chunks
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      top: "seqcont"
      molgrid_data_param {
        source: "typesfiles/grouped.types"
        batch_size: 2
        max_group_size: 4
        max_group_chunk_size: 2
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: false
        root_folder: "typesfiles"
      }
    }''')
    labels = res['label']
    seqcont = res['seqcont']
    data = res['data']
    assert labels.shape == (2,2)
    assert seqcont.shape == (2,2)
    assert list(seqcont[0]) == [0,0]
    assert list(seqcont[1]) == [1,1]
    assert data[0][0].sum() > 0
    assert data[0][1].sum() > 0
    assert data[1][0].sum() > 0
    assert data[1][1].sum() > 0
    
    res = net.forward()
    labels = res['label']
    seqcont = res['seqcont']
    data = res['data']
    assert labels.shape == (2,2)
    assert seqcont.shape == (2,2)
    assert list(seqcont[0]) == [1,1]
    assert list(seqcont[1]) == [1,1]
    assert data[0][0].sum() > 0
    assert data[0][1].sum() > 0
    assert data[1][0].sum() == 0
    assert data[1][1].sum() == 0
 
