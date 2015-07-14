"""Some PSTH-related functions used by multiple scripts"""

from __future__ import division, print_function

import sys

import pylab as pl
from pylab import get_current_fig_manager as gcfm
import numpy as np

from core import argfwhm, dist

spykepath = '/home/mspacek/dev/spyke/' # where spyke (http://spyke.github.io) is installed
sys.path.append(spykepath)
from spyke import gac

EPS = np.spacing(1) # epsilon, smallest representable non-zero number


def get_nids_psths(rec, strange, kind='responsive', blank=False,
                   binw=0.02, tres=0.0001, gauss=True, medianx=2, minthresh=3):
    """Return responsive or active nids and corresponding PSTHs for strange in rec"""
    if strange == None:
        tranges = None
    else:
        tranges = [strange] # tranges kwarg expects a list
    if kind == 'responsive': # start with all nids, then weed out unresponsive ones
        nids = rec.get_nids(tranges=tranges, kind='all')
    else: # kind is 'active'
        nids = rec.get_nids(tranges=tranges, kind='active')
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=blank, strange=strange,
                                 plot=False, binw=binw, tres=tres, gauss=gauss, norm='ntrials')
    if kind == 'responsive':
        rnids, rpsths = [], [] # nids and PSTHS to return
        for nid, psth, ts in zip(nids, psths, spikets):
            # run PSTH peak detection on this PSTH:
            baseline = medianx * np.median(psth)
            thresh = baseline + minthresh # peak detection threshold
            print("n%d" % nid, end='')
            peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh)
            npeaks = len(peakis)
            if npeaks == 0:
                continue # this PSTH has no peaks, this nid is unresponsive
            rnids.append(nid)
            rpsths.append(psth)
        print()
        return np.asarray(rnids), np.asarray(rpsths)
    else: # kind is 'active'
        return np.asarray(nids), np.asarray(psths)

def get_seps(nids, nd):
    """Build flattened array of distances between all unique pairs in nids, given neuron
    dict nd"""
    nn = len(nids)
    lti = np.tril_indices(nn, -1) # lower triangle (below diagonal) indices, ie unique pairs
    seps = []
    for nidi0, nidi1 in np.asarray(lti).T:
        sep = dist(nd[nids[nidi0]].pos, nd[nids[nidi1]].pos)
        seps.append(sep)
    seps = np.hstack(seps)
    return seps

def plot_psth(psthparams, nid, fmt='k-', alpha=0.8, ms=6, mew=2, ymax=None, yticks=None,
              figsize=(24, 7)):
    """Plot nid's corresponding PSTH in psthparams"""
    t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid]
    pl.figure(figsize=figsize)
    pl.plot(t, psth, fmt)
    # plot thresh and baseline levels:
    pl.axhline(y=thresh, c='r', ls='--')
    pl.axhline(y=baseline, c='e', ls='--')
    # mark peaks and their edges:
    if len(lis) > 0:
        pl.plot(t[lis], psth[lis], 'g+', alpha=alpha, ms=ms*1.5, mew=mew) # left edges
    if len(ris) > 0:
        pl.plot(t[ris-1], psth[ris-1], 'mx', alpha=alpha, ms=ms, mew=mew) # right edges
    if len(peakis) > 0:
        pl.plot(t[peakis], psth[peakis], 'ko', alpha=alpha, ms=ms, mec='none') # peaks
    pl.xlim(xmax=t[-1])
    pl.ylim(ymin=0, ymax=ymax)
    if yticks:
        pl.yticks(yticks)
    titlestr = 'n%d, thresh=%g, baseline=%g' % (nid, thresh, baseline)
    gcfm().window.setWindowTitle(titlestr)
    pl.gcf().tight_layout(pad=0.3) # crop figure to contents

def get_psth_peaks_gac(ts, t, psth, thresh, sigma=0.02, alpha=1.0, minpoints=5,
                       lowp=16, highp=84, checkthresh=True):
    """Extract PSTH peaks from spike times ts collapsed across trials, by clustering them
    using GAC. Then, optionally check each peak against its amplitude in the PSTH (and its
    time stamps t), to ensure it passes thresh"""

    ts2d = np.float32(ts[:, None]) # convert to 2D (one row per spike), contig float32
    # get cluster IDs and positions corresponding to spikets:
    cids, cpos = gac.gac(ts2d, sigma=sigma, alpha=alpha, minpoints=minpoints, returncpos=True)
    ucids = np.unique(cids) # unique cluster IDs
    ucids = ucids[ucids >= 0] # exclude junk cluster -1
    #npeaks = len(ucids) # but not all of them will necessarily cross the PSTH threshold
    peakis, lis, ris = [], [], []
    for ucid, pos in zip(ucids, cpos): # clusters are numbered in order of decreasing size
        spikeis, = np.where(cids == ucid)
        cts = ts[spikeis] # this cluster's spike times
        # search all spikes for argmax, same as using lowp=0 and highp=100:
        #li, ri = t.searchsorted([cts[0], cts[-1]])
        # search within percentiles for argmax:
        lt, rt = np.percentile(cts, [lowp, highp])
        # indices of all local peaks within percentiles in psth:
        li, ri = t.searchsorted([lt, rt])
        if li == ri:
            # start and end indices are identical, cluster probably falls before first or
            # after last spike time:
            assert li == 0 or li == len(psth)
            continue # no peak to be found in psth for this cluster
        localpsth = psth[li:ri]
        #allpeakiis, = argrelextrema(localpsth, np.greater)
        #if len(allpeakiis) == 0:
        #    continue # no peaks found for this cluster
        # find peakii closest to pos:
        #peakii = allpeakiis[abs((t[li + allpeakiis] - pos)).argmin()]
        # find biggest peak:
        #peakii = allpeakiis[localpsth[allpeakiis].argmax()]
        peakii = localpsth.argmax() # find max point
        if peakii == 0 or peakii == len(localpsth)-1:
            continue # skip "peak" that's really just a start or end point of localpsth
        peaki = li + peakii
        if checkthresh and psth[peaki] < thresh:
            continue # skip peak that doesn't meet thresh
        if peaki in peakis:
            continue # this peak has already been detected by a preceding, larger, cluster
        peakis.append(peaki)
        lis.append(li)
        ris.append(ri)
        print('.', end='') # indicate a peak has been found
    return np.asarray(peakis), np.asarray(lis), np.asarray(ris)

def get_psth_peaks_simple(t, psth, nid, widthmaxpoints, minthresh=3, medianx=2,
                          fwfraction=0.5):
    """Extract peaks from PSTH, simpler, faster, more robust method. Find contiguous ranges of
    baseline-exceeding points. Within each range, find the biggest value. If that value
    exceeds thresh, designate that as a peak. Slice out that baseline-exceeding range of data,
    and run argfwhm on it, with outer kwarg - ie search from the outer edges in when looking
    for FWHM. If baseline is so high that it exceeds FWHM for that peak, discard the peak.
    t is passed only so that it can be conveniently returned for plotting"""
    baseline = medianx * np.median(psth)
    thresh = baseline + minthresh # peak detection threshold
    # indices of all baseline-exceeding points in psth:
    baselineis, = np.where(psth >= (baseline+EPS)) # add EPS because baseline is often 0
    
    # find only those local peaks above baseline that are separated from each
    # other by at least one point below baseline:
    # indices into baselineis of breaks in baselineis, marking borders of contiguous ranges:
    splitis = np.where(np.diff(baselineis) > 1)[0] + 1
    if len(splitis) == 0:
        splitbaselineis = [] # this will cause nothing to happen in the for loop below
    else:
        # list of arrays of indices, representing contiguous ranges of baseline-exceeding psth:
        splitbaselineis = np.array_split(baselineis, splitis)
    peakis, lis, ris = [], [], []
    print("n%d" % nid, end='')
    for baselineis in splitbaselineis: # for each contiguous baseline-exceeding range of points
        localpsth = psth[baselineis] # slice that range of points out of the psth
        peakii = localpsth.argmax()
        if localpsth[peakii] < thresh: # peak in this range of points doesn't exceed thresh
            continue # skip to next range
        try: # get left and right FWHM indices:
            lii, rii = argfwhm(localpsth, peakii, fraction=fwfraction, method='outer')
        except ValueError: # peaki has no FWHM
            continue # skip to next range
        if (rii - lii) > widthmaxpoints: # peak is too wide
            continue # skip to next range
        offset = baselineis[0]
        peakis.append(offset + peakii)
        lis.append(offset + lii)
        ris.append(offset + rii)
        print('.', end='') # printed dot indicates a found peak
    peakis = np.asarray(peakis)
    lis = np.asarray(lis)
    ris = np.asarray(ris)

    return t, psth, thresh, baseline, peakis, lis, ris
'''
def old_get_psth_peaks(t, psth, nid):
    """Extract peaks from PSTH"""
    baseline = MEDIANX*np.median(psth)
    thresh = baseline + MINTHRESH # peak detection threshold
    #thresh = max(MINTHRESH, BASELINEX*baseline) # peak detection threshold

    # find all local peaks above thresh:
    allpeakis, = argrelextrema(psth, np.greater_equal) # indices of all local maxima in psth
    threshis, = np.where(psth >= thresh) # indices of all thresh exceeding points in psth
    peakis = np.intersect1d(threshis, allpeakis) # indices of thresh exceeding maxima

    if len(peakis) == 0:
        print('x%d' % nid, end='')
        peakis, lis, ris = None, None, None
        return t, psth, thresh, baseline, peakis, lis, ris
    else:
        print("n%d" % nid, end='')

    # sort peakis in decreasing order of peak amplitude:
    peakis = peakis[psth[peakis].argsort()[::-1]]

    # collect left and right edges of all peaks:
    lis, ris, rmpeakiis = [], [], []
    for peakii, peaki in enumerate(peakis):
        try:
            li, ri = argfwhm(psth, peaki, fraction=FWFRACTION) # left and right indices
        except ValueError: # peaki has no FWHM
            rmpeakiis.append(peakii) # mark peaki for removal
            continue # skip to next peaki
        lis.append(li)
        ris.append(ri)
    lis = np.asarray(lis)
    ris = np.asarray(ris)
    peakis = np.delete(peakis, rmpeakiis) # delete peaks marked for removal

    # For each peak, check if any others overlap with it in width. If so, discard them.
    # This rejects smaller nearby spurious peaks:
    skippeakis, keeppeakis = [], []
    keeplis, keepris = [], []
    for peaki, li, ri in zip(peakis, lis, ris):
        if peaki in skippeakis:
            continue # was trumped by a larger peak with overlapping width, skip it
        # skip peaks which fall within width of this peak:
        #skippeakiis = (peakis >= li) & (peakis <= ri) & (peakis != peaki) # bool array
        # skip peaks whose right edge overlaps with this peak's left edge, and whose
        # left edge overlaps with this peak's right edge
        skippeakiis = ((li <= ris) & (lis <= ri)) & (peakis != peaki)
        skippeakis.extend(peakis[skippeakiis])
        keeppeakis.append(peaki)
        keeplis.append(li)
        keepris.append(ri)
        print('.', end='')

    #skippeakis = np.asarray(skippeakis)
    peakis = np.asarray(keeppeakis)
    lis = np.asarray(keeplis)
    ris = np.asarray(keepris)

    return t, psth, thresh, baseline, peakis, lis, ris
'''
