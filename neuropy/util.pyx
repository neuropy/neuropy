# cython: boundscheck=True
# cython: wraparound=True
# cython: cdivision=True
# cython: profile=False

"""Some functions written in Cython for max performance"""

cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np
from libc.math cimport sqrt
# import_array() is required for access to NumPy's C API, otherwise calls to something
# like `np.PyArray_EMPTY` segfault. See:
# http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
#np.import_array()

cdef extern from "string.h" nogil:
    cdef void *memset(void *, int, size_t) # sets n bytes in memory to constant
    cdef void *malloc(size_t) # allocates without clearing to 0
    cdef void *calloc(size_t, size_t) # allocates with clearing to 0


#cdef extern from "Python.h":
#    ctypedef int Py_intptr_t

#import time
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
cdef extern from "stdio.h" nogil:
    int printf(char *, ...)
'''
cdef extern from "string.h":
    cdef void *memset(void *, int, size_t) nogil # sets n bytes in memory to constant
'''

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


def cc_tranges(np.int8_t[:, ::1] c,
               np.int64_t[::1] t,
               np.int64_t[::1] nids,
               np.int64_t[:, ::1] tranges):
    """Calculate all pairwise correlations of codes in 2D array c for every trange
    in tranges. Rows in c are neurons, columns are bins. nids are the row labels
    (neuron ids), t are the column labels (bin times)"""
    cdef int nn = c.shape[0] # number of neurons
    cdef int nt = c.shape[1] # number of time bins
    cdef int ntranges = tranges.shape[0]
    cdef int i, trangei
    cdef double m
    cdef np.int64_t[:, ::1] tis = np.searchsorted(t, tranges) # ntranges x 2 array
    np_tis = np.asarray(tis)
    # calc max widdth of a slice to save memory:
    cdef int maxslice = (np_tis[:, 1] - np_tis[:, 0]).max()
    cdef np.int8_t[:, :, ::1] cslices = np.empty((ntranges, nn, maxslice), dtype=np.int8)
    cdef double[:, ::1] means = np.zeros((ntranges, nn))
    cdef double[:, ::1] stds = np.zeros((ntranges, nn))
    cdef np.int64_t lo, hi

    for trangei in prange(ntranges, nogil=True, schedule='dynamic'):
        lo, hi = tis[trangei, 0], tis[trangei, 1]
        cslices[trangei, :] = c[:, lo:hi]
        mean_int8_axis1(cslices[trangei], means[trangei])
        std_int8_axis1(cslices[trangei], means[trangei], stds[trangei])
    print('cython means:')
    print(np.asarray(means))
    print('cython stds:')
    print(np.asarray(stds))

    '''
    #with nogil, parallel(): # need for setting up thread local buffers for prange

        means = <double *>calloc(nn, sizeof(double))
        #cdef np.float64_t[::1] stds
        #cdef np.int64_t[::1] nhigh
        #cdef int *dims = <int *>malloc(ndims*sizeof(int)) # dimension sizes
        #memset(means
        #means = np.zeros(1000)
        #stds = np.zeros(1000)
        #nhigh = np.zeros(1000, dtype=np.int64)
        for trangei in prange(ntranges, schedule='dynamic'):

            # precalculate mean and std of each cell's codetrain, rows correspond to nids:
            #means = PyArray_SimpleNew(axis, dims, np.NPY_FLOAT64)
            #stds = PyArray_SimpleNew(axis, dims, np.NPY_FLOAT64)
            #PyArray_Mean(c, 1, np.NPY_FLOAT64, means)
            #PyArray_Mean(c, 1, np.NPY_FLOAT64, stds)
            lo, hi = tis[trangei, 0], tis[trangei, 1]
            #lo = tis[trangei, 0]
            #hi = tis[trangei, 1]
            mean_int8_axis1(c[lo:hi], means)
        #if trangei == 9:
        #for i in range(nn):
        #    printf('%.3f, ', means[i])
        #printf('\n')
        #free(means)
    '''







            
    '''

            # precalculate number of high states in each neuron's code:
            #nhigh = PyArray_SimpleNew(axis, dims, np.NPY_INT64)

            uns = get_ipython().user_ns
            if uns['CODEVALS'] != [0, 1]:
                raise RuntimeError("counting of high states assumes CODEVALS = [0, 1]")
            for nii0 in range(nneurons):
                nhigh[nii0] = c[nii0].sum()
            
            #shift, shiftcorrect = self.shift, self.shiftcorrect
            #if shift and shiftcorrect:
            #    raise ValueError("only one of shift or shiftcorrect can be nonzero")

            # iterate over all pairs:
            n = self.r.n
            corrs = []
            counts = []
            pairis = []
            for nii0 in range(nneurons):
                ni0 = nids[nii0]
                for nii1 in range(nii0+1, nneurons):
                    ni1 = nids[nii1]
                    c0 = c[nii0]
                    c1 = c[nii1]
                    # (mean of product - product of means) / product of stds:
                    #numer = (c0 * c1 * binw).mean() - means[nii0] * means[nii1] * meanw
                    numer = np.dot(c0, c1) / nbins - means[nii0] * means[nii1]
                    denom = stds[nii0] * stds[nii1]
                    if numer == 0.0:
                        cc = 0.0 # even if denom is also 0
                    elif denom == 0.0: # numer is not 0, but denom is 0, prevent div by 0
                        print('skipped pair (%d, %d) in r%s' % (ni0, ni1, self.r.id))
                        continue # skip to next pair
                    else:
                        cc = numer / denom
                    # potentially shift correct using only the second spike train of each pair:
                    #if shiftcorrect:
                    #    c1sc = self.r.n[ni1].code(tranges=tranges, shift=shiftcorrect).c
                    #    ccsc = ((c0 * c1sc).mean() - means[ni0] * means[ni1]) / denom
                    #    ## TODO: might also want to try subtracting abs(ccsc)?
                    #    cc -= ccsc
                    corrs.append(cc)
                    pairis.append([nii0, nii1])
                    # take sum of high code counts of pair. Note that taking the mean wouldn't
                    # change results in self.cct(), because it would end up simply normalizing
                    # by half the value
                    counts.append(nhigh[nii0] + nhigh[nii1])
            return corrs, counts, pairis
    '''
    #print('cython mean for last trange:')
    #print(np.asarray(means))
    #print(np.asarray(means).dtype)

'''
cdef double mean_int8(np.int8_t[::1] x) nogil:
    """Return mean of 1D int8 array"""
    cdef int i
    cdef int n = x.shape[0]
    cdef double sum = 0.0
    for i in range(n):
        sum += x[i]
    return sum / n
'''
cdef void mean_int8_axis1(np.int8_t[:, ::1] x, double[::1] means) nogil:
    """Store in `means` the mean of 2D int8 array x along axis 1.
    Assume `means` is initialized to zeros"""
    cdef int i, j, m, n
    m, n = x.shape[0], x.shape[1]
    for i in range(m):
        #means[i] = 0.0 # shouldn't be any need to clear
        for j in range(n):
            means[i] += x[i, j]
        means[i] /= n

cdef void std_int8_axis1(np.int8_t[:, ::1] x, double[::1] means, double[::1] stds) nogil:
    """Store in `stds` the standard deviation of 2D int8 array x along axis 1.
    `means` holds mean of each row in x.
    Assume `means` and `stds` are initialized to zeros"""
    cdef int i, j, m, n
    cdef double d
    m, n = x.shape[0], x.shape[1]
    for i in range(m):
        #stds[i] = 0.0 # shouldn't be any need to clear
        for j in range(n):
            d = x[i, j] - means[i]
            stds[i] += d * d
        stds[i] /= n
        stds[i] = sqrt(stds[i])
