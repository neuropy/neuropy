"""Detect response events within PSTHs of specific inactive neurons, and plot them. Most of
this is copied from psth_precision.py. Run from within neuropy using `run -i
scripts/psth_precision_inactive.py`"""

from __future__ import division, print_function

from core import argfwhm, intround

rec = ptc22.tr1.r08
strange = 1550e6, np.inf # r08 synched, us, end is ~ 2300s
nids = [5, 23, 24] # 3 example inactive yet responsive nids in ptc22.tr1.r08

EPS = np.spacing(1) # epsilon, smallest representable non-zero number

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
# 2.5 Hz thresh is 1 spike in the same 20 ms wide bin every 20 trials, assuming 0 baseline:
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
FWFRACTION = 0.5 # full width fraction of max
FWHMMAX = 200 # maximum FWHM, ms
FWHMMAXPOINTS = intround(FWHMMAX / 1000 / TRES) # maximum FWHM, number of PSTH timepoints

# plotting params:
PLOTPSTH = True
FIGSIZE = 3.14, 2
YMAX = 6 # Hz
YTICKS = 0, 3, 6
MS = 5

# copied from psth_precision.py:
def plot_psth(psthparams, nid, fmt='k-', ms=MS):
    t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid]
    figure(figsize=FIGSIZE)
    pl.plot(t, psth, fmt)
    # plot thresh and baseline levels:
    axhline(y=thresh, c='r', ls='--')
    axhline(y=baseline, c='e', ls='--')
    # mark peaks and their edges:
    #if lis != None:
    #    pl.plot(t[lis], psth[lis], 'co', ms=ms, mec='none') # left edges
    #if ris != None:
    #    pl.plot(t[ris], psth[ris], 'bo', ms=ms, mec='none') # right edges
    if peakis != None:
        pl.plot(t[peakis], psth[peakis], 'ro', ms=ms, mec='none') # peaks
    xlim(xmax=t[-1])
    ylim(ymin=0, ymax=YMAX)
    yticks(YTICKS)
    gcfm().window.setWindowTitle('n%d, thresh=%g, baseline=%g' % (nid, thresh, baseline))
    gcf().tight_layout(pad=0.3) # crop figure to contents

# copied from psth_precision.py:
def get_psth_peaks(t, psth, nid):
    """Extract peaks from PSTH, simpler, faster, more robust method. Find contiguous ranges of
    baseline-exceeding points. Within each range, find the biggest value. If that value
    exceeds thresh, designate that as a peak. Slice out that baseline-exceeding range of data,
    and run argfwhm on it, with outer kwarg - ie search from the outer edges in when looking
    for FWHM. If baseline is so high that it exceeds FWHM for that peak, discard the peak."""
    baseline = MEDIANX*np.median(psth)
    thresh = baseline + MINTHRESH # peak detection threshold
    # indices of all baseline-exceeding points in psth:
    baselineis, = np.where(psth >= (baseline+EPS)) # add EPS because baseline is often 0
    
    # find only those local peaks above baseline that are separated from each
    # other by at least one point below baseline:
    # indices into baselineis of breaks in baselineis, marking borders of contiguous ranges:
    splitis = np.where(np.diff(baselineis) > 1)[0] + 1
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
            lii, rii = argfwhm(localpsth, peakii, fraction=FWFRACTION, method='outer')
        except ValueError: # peaki has no FWHM
            continue # skip to next range
        if (rii - lii) > FWHMMAXPOINTS: # peak is too wide
            continue # skip to next range
        offset = baselineis[0]
        peakis.append(offset + peakii)
        lis.append(offset + lii)
        ris.append(offset + rii)
        print('.', end='') # printed dot indicates a found peak
    peakis = np.asarray(peakis)
    lis = np.asarray(lis)
    ris = np.asarray(ris)

    if len(peakis) == 0:
        peakis, lis, ris = None, None, None

    return t, psth, thresh, baseline, peakis, lis, ris


t, psths = rec.psth(nids=nids, natexps=False, strange=strange, plot=False,
                    binw=BINW, tres=TRES, norm='ntrials')

psthparams = {} # params returned for each PSTH by get_psth_peaks
for nid, psth in zip(nids, psths):
    psthparams[nid] = get_psth_peaks(t, psth, nid)
    if PLOTPSTH:
        plot_psth(psthparams, nid, fmt='k-')
    #t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid] # unpack

pl.show()
