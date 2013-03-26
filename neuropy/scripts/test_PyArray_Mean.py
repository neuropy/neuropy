"""Test PyArray_Mean.pyx"""

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import PyArray_Mean # .pyx file

import numpy as np
import time

#a = np.arange(20, dtype=np.float64)
#a.shape = 4, 5
a = np.random.random((25, 5000000))

t0 = time.time()
print(a.mean(axis=1))
print('a.mean(axis=1) took %.3f sec' % (time.time()-t0))
t0 = time.time()
print(PyArray_Mean.mean(a, axis=1))
print('PyArray_Mean.mean(a, axis=1) took %.3f sec' % (time.time()-t0))

# both methods are equally fast in this trivial kind of a call
