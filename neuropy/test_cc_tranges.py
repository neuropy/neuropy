"""Test util.cc_tranges"""

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

import numpy as np
import time

nn = 5 # number of neurons
nt = 100
c = np.int8(np.random.randint(2, size=(nn, nt)))
print(c)

nids = np.arange(nn, dtype=np.int64)
t = np.arange(nt, dtype=np.int64)
#ntranges = 10
tranges = np.array([[0, 10], [10,20], [20,30], [30,40], [40,50],
                    [50,60], [60,70], [70,80], [80,90], [90,100]], dtype=np.int64)
ntranges = len(tranges)
tis = t.searchsorted(tranges) # ntranges x 2 array
means = np.zeros((ntranges, nn))
for i, (lo, hi) in enumerate(tis):
    means[i] = c[:, lo:hi].mean(axis=1)
print('python means:')
print(means)

#print(c.mean(axis=1))

util.cc_tranges(c, t, nids, tranges)
