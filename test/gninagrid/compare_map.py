#!/usr/bin/env python3

'''Compare two map files'''

import sys
import pytest

from pytest import approx

f1 = open(sys.argv[1])
f2 = open(sys.argv[2])

#6 lines of Header
for i in range(6):
    l1 = f1.readline()
    l2 = f2.readline()
    assert l1 == l2
    
#single flots
for l1 in f1:
    l2 = f2.readline()
    val1 = float(l1)
    val2 = float(l2)
    assert val1 == approx(val2,abs=1e-4)
