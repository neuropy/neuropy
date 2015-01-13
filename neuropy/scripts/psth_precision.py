"""Detect response events within PSTHs, count them, measure the FWHM of each one, and plot
their distributions as a function of cortical state within each of the 2 natural scene movies
in ptc22.tr1"""

from __future__ import division
from scipy.signal import argrelextrema
from core import argfwhm, get_ssnids

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
## TODO: maybe reduce to 3?
MINTHRESH = 3 # peak detection thresh, Hz
#BASELINEX = 5 # PSTH baseline multiplier, Hz
FWFRACTION = 0.5 # full width fraction of max

# plotting params:
FWHMMIN, FWHMMAX, FWHMSTEP, XTICKSTEP = 0, 100, 2.5, 25
figsize = (3.5, 3.5) # inches

def get_psth_peaks(psthi, psth, plot=True):
    """Extract peaks from PSTH"""
    baseline = np.median(psth)
    thresh = baseline + MINTHRESH # peak detection threshold
    #thresh = max(MINTHRESH, BASELINEX*baseline) # peak detection threshold

    # find all local peaks above thresh:
    allpeakis, = argrelextrema(psth, np.greater_equal) # indices of all local maxima in psth
    threshis, = np.where(psth >= thresh) # indices of all thresh exceeding points in psth
    peakis = np.intersect1d(threshis, allpeakis) # indices of thresh exceeding maxima
    '''
    # alternatively, find only those local peaks above thresh that are separated from each
    # other by at least one point below threshold:
    splitis = np.where(np.diff(threshis) > 1)[0] + 1 # indices of thresh exceeding ranges
    splitthreshis = np.array_split(threshis, splitis) # list of arrays of contiguous indices
    peakis = []
    for threshis in splitthreshis: # for each contiguous thresh-exceeding range of points
        localpeakis = np.intersect1d(threshis, allpeakis) # indices of thresh exceeding maxima
        if len(localpeakis) == 0:
            continue # skip this range of points
        # keep only biggest maximum in this contiguous range:
        peakii = np.argmax(psth[localpeakis])
        peaki = localpeakis[peakii]
        peakis.append(peaki)
    peakis = np.asarray(peakis)
    '''
    if len(peakis) == 0:
        raise ValueError("psthi %d has no peaks" % psthi)

    # sort peakis in decreasing order of peak amplitude:
    peakis = peakis[psth[peakis].argsort()[::-1]]

    # collect left and right edges of all peaks:
    lis, ris = [], []
    for peaki in peakis:
        ## TODO: maybe FWHM should be taken relative to baseline instead of 0.
        ## Or, is it more fair to measure FWHM relative to 0, not to some baseline that's
        ## different for each cell?
        li, ri = argfwhm(psth, peaki, fraction=FWFRACTION) # left and right indices into psth
        lis.append(li)
        ris.append(ri)
    lis = np.asarray(lis)
    ris = np.asarray(ris)

    # For each peak, check if any others overlap with it in width. If so, discard them.
    # This rejects smaller nearby spurious peaks:
    skippeakis, keeppeakis = [], []
    keeplis, keepris = [], []
    fwhms = []
    for peaki, li, ri in zip(peakis, lis, ris):
        if peaki in skippeakis:
            continue # was trumped by a larger peak with overlapping width, skip it
        # skip peaks which fall within width of this peak:
        #skippeakiis = (peakis >= li) & (peakis <= ri) & (peakis != peaki) # bool array
        # skip peaks whose right edge overlaps with this peak's left edge, and whose
        # left edge overlaps with this peak's right edge
        skippeakiis = ((li <= ris) & (lis <= ri)) & (peakis != peaki)
        skippeakis.extend(peakis[skippeakiis])
        fwhm = (ri - li)
        fwhms.append(fwhm)
        if plot:
            keeppeakis.append(peaki)
            keeplis.append(li)
            keepris.append(ri)

    fwhms = np.asarray(fwhms) * TRES * 1000 # sec

    if plot:
        figure(figsize=(23, 7))
        plot(psth, 'k-')
        plot(skippeakis, psth[skippeakis], 'eo', mec='none') # discarded peaks
        ms = 10
        plot(keeplis, psth[keeplis], 'co', ms=ms, mec='none') # left edges
        plot(keepris, psth[keepris], 'bo', ms=ms, mec='none') # rigth edges
        plot(keeppeakis, psth[keeppeakis], 'ro', ms=ms, mec='none') # kept peaks
        gcf().tight_layout(pad=0.3) # crop figure to contents
        gcfm().window.setWindowTitle('%d, thresh=%f, baseline=%f' % (psthi, thresh, baseline))

    return fwhms

# get active neuron ids for each section of r08:
ssnids, recsecnids = get_ssnids(ptc22tr1r08s, strangesr08s)
## TODO: use responsive nids instead of active nids used above?

# calculate PSTHs for both sections of r08:
midbinss, psthss = [], []
for rec, nids, strange in zip(ptc22tr1r08s, recsecnids, strangesr08s):
    midbins, psths = rec.traster(nids=nids, natexps=False, strange=strange, plot=False,
                                 psth=True, binw=BINW, tres=TRES, norm='ntrials')
    midbinss.append(midbins)
    psthss.append(psths)
    #figure()
    #plot(midbins, psths.T, '-')

fwhms = {} # fwhm values for each recording section
for recseci, psths in enumerate(psthss):
    psthsfwhms = [] # fwhm values across PSTHs from this recording section
    for psthi, psth in enumerate(psths):
        try:
            psthfwhms = get_psth_peaks(psthi, psth, plot=False)
        except ValueError:
            continue # this psth has no peaks, skip it
        psthsfwhms.append(psthfwhms)
    fwhms[recseci] = np.hstack(psthsfwhms)

ticks = np.arange(FWHMMIN, FWHMMAX, XTICKSTEP)
bins = np.arange(FWHMMIN, FWHMMAX+FWHMSTEP, FWHMSTEP)
figure(figsize=figsize)
## TODO: try plotting distributions on log scale
n1 = hist(fwhms[1], bins=bins, color='r')[0] # second section in r08, synched
n0 = hist(fwhms[0], bins=bins, color='b')[0] # first section in r08, desynched
n = np.hstack([n0, n1])
xlim(xmin=FWHMMIN, xmax=FWHMMAX)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('FWHM (ms)')
ylabel('PSTH peak count')
gcfm().window.setWindowTitle('PSTH precision ptc22.tr1.r08')
tight_layout(pad=0.3)

show()
