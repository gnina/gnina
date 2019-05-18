#!/usr/bin/env python2

'''Test basic molgrid functionality.  More concerned with how sampling is done than gridding'''

import caffe
import numpy as np
import os,sys,re

import pytest
from pytest import approx

def test_molgrid():

    caffe.set_random_seed(800)
    caffe.set_mode_gpu()
    #these might change
    LIGCHANNELS = 19
    RECCHANNELS = 16
    
    
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
      molgrid_data_param {
        source: "typesfiles/small.types"
        batch_size: 10
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: false
        root_folder: "typesfiles"
      }
    }''')
    labels = res['label']
    
    assert list(labels) == [1., 1., 1., 1., 0., 0., 0., 0., 0., 0.]
    
    labels2 = net.forward()['label']
    assert list(labels2) == [1., 1., 1., 1., 0., 0., 0., 0., 0., 0.]
    
    #balanced
    labels = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      molgrid_data_param {
        source: "typesfiles/small.types"
        batch_size: 10
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: true
        root_folder: "typesfiles"
      }
    }''')[0]['label']
    
    assert list(labels) == [1., 0., 1., 0., 1., 0., 1., 0., 1., 0.]
    
    
    #shuffled
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      molgrid_data_param {
        source: "typesfiles/small.types"
        batch_size: 10
        dimension: 23.5
        resolution: 0.5
        shuffle: true
        balanced: false
        root_folder: "typesfiles"
      }
    }''')
    labels = list(res['label'])
    labels2 = list(net.forward()['label'])
    assert labels != labels2
    
    #stratify by receptor
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      top: "affinity"
      molgrid_data_param {
        source: "typesfiles/smallaff.types"
        batch_size: 8
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: false
        has_affinity: true
        stratify_receptor: true
        root_folder: "typesfiles"
      }
    }''')
    
    assert list(res['affinity']) == approx([ 3.2, -5. ,  6. , -4. ,  4.2,  2. , -8.2, -1. ])
    
    #balanced stratified - will drop some receptors without both true/false examples
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      top: "affinity"
      molgrid_data_param {
        source: "typesfiles/smallaff.types"
        batch_size: 8
        dimension: 23.5
        resolution: 0.5
        shuffle: false
        balanced: true
        has_affinity: true
        stratify_receptor: true
        root_folder: "typesfiles"
      }
    }''')
    assert list(res['affinity']) == approx([-5. ,  2. ,  6. , -8.2, -5. , -3. ,  6. , -8.2])
    assert list(res['label']) == approx([1., 0., 1., 0., 1., 0., 1., 0.])
    
    #stratify by affinity, also shuffle
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      top: "affinity"
      molgrid_data_param {
        source: "typesfiles/smallaff.types"
        batch_size: 8
        dimension: 23.5
        resolution: 0.5
        shuffle: true
        balanced: false
        has_affinity: true
        stratify_affinity_min: 0
        stratify_affinity_max: 8
        stratify_affinity_step: 4
        root_folder: "typesfiles"
      }
    }''')
    
    assert np.sum(np.abs(res['affinity']) < 4) == 4
    
    #stratify by affinity and balance, also shuffle
    res, net = runmodel('''layer {
      name: "data"
      type: "MolGridData"
      top: "data"
      top: "label"
      top: "affinity"
      molgrid_data_param {
        source: "typesfiles/smallaff.types"
        batch_size: 8
        dimension: 23.5
        resolution: 0.5
        shuffle: true
        balanced: true
        has_affinity: true
        stratify_affinity_min: 0
        stratify_affinity_max: 8
        stratify_affinity_step: 4
        root_folder: "typesfiles"
      }
    }''')
    
    assert np.sum(np.abs(res['affinity']) < 4) == 4
    assert np.sum(res['label']) == 4
