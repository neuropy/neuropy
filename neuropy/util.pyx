# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: profile=False

"""Some functions written in Cython for max performance"""

cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np
from numpy cimport int8_t, int64_t, float64_t
from libc.math cimport sqrt
# import_array() is required for access to NumPy's C API, otherwise calls to something
# like `np.PyArray_EMPTY` segfault. See:
# http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
#np.import_array()
'''
cdef extern from "string.h" nogil:
    cdef void *memset(void *, int, size_t) # sets n bytes in memory to constant
    cdef void *malloc(size_t) # allocates without clearing to 0
    cdef void *calloc(size_t, size_t) # allocates with clearing to 0
'''
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

cdef float64_t listmean(x):
    """This demonstrates that you can't access a Python object within a prange
    without gil"""
    cdef int i, n = len(x)
    cdef float64_t mean = 0.0
    for i in prange(n, nogil=True, schedule='dynamic'):
        mean += x[i] # Cythonization error
    return mean
'''
def list_int_arrs(arrs):
    """This demonstrates how to deal with a (ragged) list of arrays without gil by
    concatenating them into a single big one, and then indexing into that big one
    appropriately from within a prange loop"""
    cdef int64_t i, start, narrs = len(arrs)
    cdef int64_t result = 0
    cdef int64_t[::1] flat = np.concatenate(arrs) # all arrays flattened into one long one
    cdef int64_t[::1] n = np.int64([ len(arr) for arr in arrs ])
    cdef int64_t[::1] offsets = np.concatenate([[0], np.cumsum(n)[:-1]])
    for i in prange(narrs, nogil=True, schedule='dynamic'):
        start = offsets[i]
        result += int_sum(flat[start:start+n[i]])
    return result

cdef int64_t int_sum(int64_t[::1] x) nogil:
    """Sum an int array"""
    cdef int64_t i, result = 0
    # note that calling np.anything within a nogil block doesn't work:
    #cdef int64_t[::1] test = np.zeros(10)
    #test = np.resize(test, (20,))
    for i in range(x.shape[0]):
        result += x[i]
    return result


def xcorr(int64_t[::1] x,
          int64_t[::1] y,
          int64_t[::1] trange):
    """Calculate cross-correlation of timepoints in x with y, constrained to lower
    and upper bounds in trange. Assume timepoints in x and y are sorted. Return spike times
    of y relative to x."""
    cdef int64_t ntx, nty, loti, dtsi, xti, yti, maxxti, maxyti, t, dt
    cdef int64_t low = trange[0]
    cdef int64_t high = trange[1]
    cdef int64_t DTSALLOCSIZE = 1000000
    ntx = x.shape[0]
    nty = y.shape[0]
    maxxti = ntx - 1
    maxyti = nty - 1
    cdef int64_t[::1] dts = np.zeros(DTSALLOCSIZE, dtype=np.int64)
    cdef int64_t maxdtsi = dts.shape[0] - 1

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
        dt = y[yti] - t # dt is y relative to x
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
    return np.asarray(dts[:dtsi]) # trim it down, convert memory view slice to array


def sct(int8_t[:, ::1] c,
        int64_t[::1] t,
        int64_t[:, ::1] tranges,
        int8_t highval):
    """Calculate all pairwise spike correlations of codes in 2D array c for every trange
    in tranges. Rows in c are neurons, columns are time bins. t are the bin times"""
    cdef int64_t nn = c.shape[0] # number of neurons
    cdef int64_t nt = c.shape[1] # number of time bins
    cdef int64_t ntranges = tranges.shape[0]
    cdef int64_t i, j, trangei, sti
    cdef int64_t[:, ::1] tis = np.searchsorted(t, tranges) # ntranges x 2 array
    np_tis = np.asarray(tis)
    # calc number of slice time bins for each trange:
    np_nst = np_tis[:, 1] - np_tis[:, 0]
    cdef int64_t maxslice = np_nst.max()
    cdef int64_t[::1] nst = np_nst # has length ntranges
    cdef int8_t[:, :, ::1] cslices = np.empty((ntranges, nn, maxslice), dtype=np.int8)
    cdef float64_t[:, ::1] means = np.zeros((ntranges, nn))
    cdef float64_t[:, ::1] stds = np.zeros((ntranges, nn))
    cdef int64_t[:, ::1] nhigh = np.zeros((ntranges, nn), dtype=np.int64)
    cdef int64_t lo, hi

    # pre-calc some arrays for use in main trange loop
    for trangei in prange(ntranges, nogil=True, schedule='dynamic'):
        lo, hi = tis[trangei, 0], tis[trangei, 1]
        cslices[trangei, :, :nst[trangei]] = c[:, lo:hi]
        mean_int8_axis1(cslices[trangei], nst[trangei], means[trangei])
        std_int8_axis1(cslices[trangei], nst[trangei], means[trangei], stds[trangei])
        # count up number of high states for each neuron in each trange, used later
        # for weighted average of sc(t) across neurons:
        for i in range(nn):
            for sti in range(nst[trangei]):
                if cslices[trangei, i, sti] == highval:
                    nhigh[trangei, i] += 1
    '''
    print('cython:')
    print('nst:')
    print(np.asarray(nst))
    print('cslices:')
    print(np.asarray(cslices))
    print('means:')
    print(np.asarray(means))
    print('stds:')
    print(np.asarray(stds))
    print('nhigh:')
    print(np.asarray(nhigh))
    '''
    cdef int64_t npairs = nn * (nn - 1) / 2
    cdef float64_t[:, ::1] corrs = np.empty((ntranges, npairs))
    cdef int64_t[:, ::1] counts = np.empty((ntranges, npairs), dtype=np.int64)
    cdef int64_t pairi, sumprod # sum of products
    cdef float64_t numer, denom
    """pairi and sumprod can't be incremented in place in this prange, because of Cython's
    automatic inferral of thead-locality and reductions. From the docs: "If you use an
    inplace operator on a variable, it becomes a reduction, meaning that the values from the
    thread-local copies of the variable will be reduced with the operator and assigned to
    the original variable after the loop." After incrementing (and causing a reduction),
    Cython doesn't allow reading the variable later within the loop, raising this error:
    "Cannot read reduction variable in loop body". Assigning to it initially in a `with
    parallel()` block doesn't help - that only works for buffers (I think). There are two
    ways around this currently: 1) don't do in-place operations; 2) abstract the whole loop
    out into its own function, where the thread-local and reductions rules don't apply. For
    now, I've done the former. See:
    * http://docs.cython.org/src/userguide/parallelism.html
    * https://groups.google.com/forum/?fromgroups=#!topic/cython-users/Ady-DdWu6rE
    """
    for trangei in prange(ntranges, nogil=True, schedule='dynamic'):
        pairi = 0
        for i in range(nn):
            for j in range(i+1, nn):
                # accumulate element-wise product of c[trangei] for this neuron pair:
                sumprod = 0
                for sti in range(0, nst[trangei]):
                    #sumprod += cslices[trangei, i, sti] * cslices[trangei, j, sti]
                    sumprod = sumprod + cslices[trangei, i, sti] * cslices[trangei, j, sti]
                # (mean of product - product of means) / product of stds:
                numer = (<float64_t>sumprod / nst[trangei]
                         - means[trangei, i] * means[trangei, j])
                denom = stds[trangei, i] * stds[trangei, j]
                # all codes values for at least one neuron must've been identical,
                # leading to 0 std, call that 0 code correlation:
                if denom == 0.0:
                    corrs[trangei, pairi] = 0.0
                else:
                    corrs[trangei, pairi] = numer / denom
                # store sum of high code counts of this pair:
                counts[trangei, pairi] = nhigh[trangei, i] + nhigh[trangei, j]
                #pairi += 1 # inc for next loop
                pairi = pairi + 1 # inc for next loop

    return np.asarray(corrs.T), np.asarray(counts.T) # pairs in rows, tranges in columns

'''
cdef double mean_int8(int8_t[::1] x) nogil:
    """Return mean of 1D int8 array"""
    cdef int i
    cdef int n = x.shape[0]
    cdef double sum = 0.0
    for i in range(n):
        sum += x[i]
    return sum / n
'''
cdef void mean_int8_axis1(int8_t[:, ::1] x, int64_t n, float64_t[::1] means) nogil:
    """Store in `means` the mean of 2D int8 array x along axis 1.
    Consider only the first `n` values of x along axis 1.
    Assume `means` is initialized to zeros"""
    cdef int i, j, m
    m = x.shape[0]
    for i in range(m):
        #means[i] = 0.0 # shouldn't be any need to clear
        for j in range(n):
            means[i] += x[i, j]
        means[i] /= n

cdef void std_int8_axis1(int8_t[:, ::1] x, int64_t n, float64_t[::1] means,
                         float64_t[::1] stds) nogil:
    """Store in `stds` the standard deviation of 2D int8 array x along axis 1.
    Consider only the first `n` values of x along axis 1.
    `means` holds mean of each row in x.
    Assume `stds` is initialized to zeros"""
    cdef int i, j, m
    cdef double d
    m = x.shape[0]
    for i in range(m):
        #stds[i] = 0.0 # shouldn't be any need to clear
        for j in range(n):
            d = x[i, j] - means[i]
            stds[i] += d * d
        stds[i] /= n
        stds[i] = sqrt(stds[i])
