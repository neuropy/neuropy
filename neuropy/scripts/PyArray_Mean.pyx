# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: profile=False

"""Show how to call numpy functions from within Cython via their C API"""

cimport cython
#from cython.parallel import prange#, parallel
import numpy as np
cimport numpy as np
# import_array() is required for access to NumPy's C API, otherwise calls to something
# like `np.PyArray_EMPTY` segfault. See:
# http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
np.import_array()


def mean(np.ndarray[np.float64_t, ndim=2, mode='c'] a, int axis=0):
    cdef int nrows = a.shape[0]
    #cdef int ncols = a.shape[1]
    cdef np.npy_intp *dims = [nrows]
    #cdef np.ndarray[np.float64_t, ndim=1] out = np.PyArray_EMPTY(1, dims, np.NPY_FLOAT64, 0)
    cdef np.ndarray[np.float64_t, ndim=1] out = np.PyArray_SimpleNew(axis, dims, np.NPY_FLOAT64)
    #return np.PyArray_Mean(a, 1, np.NPY_FLOAT64, out) # this works too
    np.PyArray_Mean(a, 1, np.NPY_FLOAT64, out)
    return out
