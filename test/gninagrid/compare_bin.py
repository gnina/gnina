#!/usr/bin/env python3

'''Compare two binmap files which are nothing but floating point numbers'''

import sys,struct
import pytest

from pytest import approx

f1 = open(sys.argv[1],'rb')
f2 = open(sys.argv[2],'rb')

buf1 = f1.read()
buf2 = f2.read()

assert len(buf1) == len(buf2)
n = int(len(buf1)/4)

vals1 = struct.unpack('f'*n,buf1)
vals2 = struct.unpack('f'*n,buf2)

assert vals1 == approx(vals2,abs=1e-4)


