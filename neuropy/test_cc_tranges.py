"""Test util.cc_tranges"""

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

import numpy as np
import time

#a = np.arange(20, dtype=np.float64)
#a.shape = 4, 5
m = 5
n = 10
c = np.zeros((m, n), dtype=np.int8)
c[0,0] = 1
c[4,7] = 1
c[4,4] = 1

t = np.arange(n, dtype=np.int64)
nids = np.arange(m, dtype=np.int64)
ntranges = 100
tranges = np.arange(ntranges*2, dtype=np.int64)
tranges.shape = ntranges, 2

util.cc_tranges(c, t, nids, tranges)
