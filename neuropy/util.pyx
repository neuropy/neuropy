"""Some functions written in Cython for max performance"""

cimport cython
from cython.parallel import prange#, parallel
import numpy as np
cimport numpy as np

import time
'''
cdef extern from "math.h":
    int abs(int x)
    float fabs(float x)
    double ceil(double x) nogil

cdef extern from "limits.h":
    int INT_MAX

cdef extern from "float.h":
    double DBL_MAX
'''
cdef extern from "stdio.h":
    int printf(char *, ...)
'''
cdef extern from "string.h":
    cdef void *memset(void *, int, size_t) nogil # sets n bytes in memory to constant
'''

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # might be necessary to release the GIL?
@cython.profile(False)
def xcorr(np.ndarray[np.int64_t, ndim=1, mode='c'] x,
          np.ndarray[np.int64_t, ndim=1, mode='c'] y,
          np.ndarray[np.int64_t, ndim=1, mode='c'] trange):
    """Calculate cross-correlation of timepoints in x with y, constrained to lower
    and upper bounds in trange. Assume timepoints in x and y are sorted"""
    # should assert contig of x and y, this seems to happen automatically though
    cdef long long ntx, nty, loti, dtsi, xti, yti, maxxti, maxyti, t, dt
    cdef long long low = trange[0]
    cdef long long high = trange[1]
    cdef long long DTSALLOCSIZE = 1000000
    ntx = x.shape[0]
    nty = y.shape[0]
    maxxti = ntx - 1
    maxyti = nty - 1
    cdef np.ndarray[np.int64_t, ndim=1] dts = np.zeros(DTSALLOCSIZE, dtype=np.int64)
    cdef long long maxdtsi = dts.shape[0] - 1

    loti = 0
    dtsi = 0
    for xti in range(ntx):
        # t is current timepoint in x to compare to all timepoints in y:
        t = x[xti]
        while y[loti] - t < low: # keep checking lower trange bound
            loti += 1
            if loti > maxyti: # no y timepoints fall within trange of t
                break
        # start collecting dt values:
        if loti > maxyti: # no y timepoints fall within trange of t
            continue # to next xti
        yti = loti
        dt = y[yti] - t
        while dt < high: # keep checking upper trange bound
            if dtsi > maxdtsi:
                # when growing an array, pretty much need to allocate a new one,
                # can't very often do it in place:
                dts = np.resize(dts, (dts.shape[0] + DTSALLOCSIZE,))
                maxdtsi = dts.shape[0] - 1
                printf('resized dts array to %d entries\n', dts.shape[0])
            dts[dtsi] = dt
            #printf('%d ', dtsi)
            dtsi += 1 # inc for next loop iter
            yti += 1
            if yti > maxyti: # don't exceed maxyti when indexing into y
                break
            dt = y[yti] - t # update for next loop iter
    dts = dts[:dtsi] # trim it down
    return dts

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # might be necessary to release the GIL?
@cython.profile(False)
def intersect1d_uint8(arrs):
    """Find the intersection of any number of 1D arrays in arrs list.
    Return the sorted, unique values that are in all of the input arrays.
    This is a much faster (at least for many arrs) but type-specific version of
    core.intersect1d()"""
    cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] arr
    cdef np.ndarray[np.uint8_t, ndim=1, mode='c'] common = np.unique(arrs[0])
    cdef int ncommon=common.shape[0]
    cdef int lenarr
    cdef int i, j, k
    cdef bint continuei=False
    #print('common: %s' % common)
    #print('ncommon: %d' % ncommon)
    for arr in arrs[1:]:
        lenarr = arr.shape[0]
        #print('arr: %s' % arr)
        i = 0 # reset
        while i < ncommon:
            #print('i = %d' % i)
            ## TODO: this naive search could be sped up by assuming sorted arrs, or by
            ## first sorting each arr with the quicksort algorithm, although that would
            ## modify each array in place, and so the caller would need to send a copy
            ## to prevent modifying the originals:
            for j in range(lenarr):
                #print('j = %d' % j)
                if common[i] == arr[j]:
                    #print('breaking out of j')
                    continuei = True
                    break # out of j loop
            if continuei:
                continuei = False # reset
                i += 1
                continue # to next i
            # never broke out of j loop, didn't find a value in current arr that matches
            # common[i], common[i] is no longer common:
            ncommon -= 1
            for k in range(i, ncommon):
                common[k] = common[k+1] # shift values above i down by 1
            #print('new common: %s' % common)
            #print('new ncommon: %d' % ncommon)
            # don't inc i, new value at common[i] has just shifted into view
    return common[:ncommon]
