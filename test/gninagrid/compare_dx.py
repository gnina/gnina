#!/usr/bin/env python3

'''Compare two dx files'''

import sys
import pytest

from pytest import approx

f1 = open(sys.argv[1])
f2 = open(sys.argv[2])

#7 lines of Header
for i in range(7):
    l1 = f1.readline()
    l2 = f2.readline()
    assert l1 == l2
    
#triplets of floats
for l1 in f1:
    l2 = f2.readline()
    vals1 = list(map(float,l1.split()))
    vals2 = list(map(float,l2.split()))
    assert vals1 == approx(vals2,abs=1e-4)
