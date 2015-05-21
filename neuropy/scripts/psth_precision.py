"""Detect response events within PSTHs, count them, measure the FWHM of each one, and plot
their distributions as a function of cortical state within each of the natural scene movie
recordings listed below.

Also, measure trial raster plot reliability as a function of cortical state, by taking all
pairwise correlations between all trials for a given neuron, and then averaging them to get a
single number between 0 and 1 representing reliability of that trial raster plot. See Goard &
Dan, 2009

Run from within neuropy using `run -i scripts/psth_precision.py`"""

from __future__ import division, print_function
from scipy.signal import argrelextrema
from scipy.stats import ttest_ind, chisquare, mannwhitneyu
from numpy import log10

from core import argfwhm, get_ssnids, sparseness, intround, ceilsigfig

spykepath = '/home/mspacek/dev/spyke/' # where spyke (http://spyke.github.io) is installed
sys.path.append(spykepath)
from spyke import gac

# mapping of recording to list of desynched and synched trange, in that order:
rec2tranges = {ptc17.tr2b.r58: [(0, 700e6), # desynched trange, 66 Hz refresh rate
                                (800e6, np.inf)], # synched trange, 66 Hz refresh rate
               ptc18.tr1.r38:  [(0, 425e6), # desynched trange, ends ~ trial 76
                                (550e6, np.inf)], # synched trange, starts ~ trial 98
               ptc18.tr2c.r58: [(0, 750e6), # desynched trange
                                (1000e6, np.inf)], # synched trange
               ptc22.tr1.r08:  [(0, 1500e6), # desynched trange
                                (1550e6, np.inf)], # synched trange
               ptc22.tr1.r10:  [(1480e6, np.inf), # desynched trange
                                (0, 1400e6)], # synched trange
               ptc22.tr4b.r49: [(0, 1475e6), # desynched trange
                                (1500e6, np.inf)], # synched trange
              }
# compare and sort recordings by their absname:
reccmp = lambda reca, recb: cmp(reca.absname, recb.absname)
urecs = sorted(rec2tranges, cmp=reccmp) # unique recordings, no repetition
urecnames = ' '.join([rec.absname for rec in urecs])

EPS = np.spacing(1) # epsilon, smallest representable non-zero number

NIDSKIND = 'all' # 'active' or 'all'

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH by convolving collapsed spike train with Gaussian kernel?
TRASTERBINW, TRASTERTRES = 0.02, 0.005 # trial raster bins, sec
BLANK = False # consider blank periods between trials?
WEIGHT = False # weight trials by spike count for reliability measure?
# 2.5 Hz thresh is 1 spike in the same 20 ms wide bin every 20 trials, assuming 0 baseline:
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
FWFRACTION = 0.5 # full width fraction of max
FWHMMAX = 200 # maximum FWHM, ms
FWHMMAXPOINTS = intround(FWHMMAX / 1000 / TRES) # maximum FWHM, number of PSTH timepoints

# plotting params:
PLOTPSTH = False
FWHMMIN, FWHMSTEP, FWHMTICKSTEP = 0, 10, 50
TSMAX, TSSTEP = 5.5, 0.25
HEIGHTMIN, HEIGHTMAX, HEIGHTSTEP, HEIGHTTICKSTEP = 0, 100, 5, 25
NSPARSBINS = 15
NRELBINS = 15
LOGNULLREL = -5
NULLREL = 10**LOGNULLREL
figsize = (3, 3) # inches

# copied to psth_precision_inactive.py:
def plot_psth(psthparams, nid, fmt='k-', ms=10):
    t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid]
    figure(figsize=(24, 7))
    pl.plot(t, psth, fmt)
    # plot thresh and baseline levels:
    axhline(y=thresh, c='r', ls='--')
    axhline(y=baseline, c='e', ls='--')
    # mark peaks and their edges:
    if lis != None:
        pl.plot(t[lis], psth[lis], 'co', ms=ms, mec='none') # left edges
    if ris != None:
        pl.plot(t[ris], psth[ris], 'bo', ms=ms, mec='none') # right edges
    if peakis != None:
        pl.plot(t[peakis], psth[peakis], 'ro', ms=ms, mec='none') # peaks
    xlim(xmax=t[-1])
    ylim(ymin=0)
    gcfm().window.setWindowTitle('n%d, thresh=%g, baseline=%g' % (nid, thresh, baseline))
    gcf().tight_layout(pad=0.3) # crop figure to contents
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

# copied to psth_precision_inactive.py:
def get_psth_peaks_simple(t, psth, nid):
    """Extract peaks from PSTH, simpler, faster, more robust method. Find contiguous ranges of
    baseline-exceeding points. Within each range, find the biggest value. If that value
    exceeds thresh, designate that as a peak. Slice out that baseline-exceeding range of data,
    and run argfwhm on it, with outer kwarg - ie search from the outer edges in when looking
    for FWHM. If baseline is so high that it exceeds FWHM for that peak, discard the peak.
    t is passed only so that it can be conveniently returned for plotting"""
    baseline = MEDIANX * np.median(psth)
    thresh = baseline + MINTHRESH # peak detection threshold
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
'''
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
    for ucid, pos in zip(ucids, cpos):
        spikeis, = np.where(cids == ucid)
        cts = ts[spikeis] # this cluster's spike times
        # search all spikes for argmax, same as using lowp=0 and highp=100:
        #li, ri = t.searchsorted([cts[0], cts[-1]])
        lt, rt = np.percentile(cts, [lowp, highp])
        li, ri = t.searchsorted([lt, rt]) # search within percentiles for argmax
        # indices of all local peaks within percentiles in psth
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
        peakis.append(peaki)
        lis.append(li)
        ris.append(ri)
        print('.', end='') # indicate a peak has been found
    return np.asarray(peakis), np.asarray(lis), np.asarray(ris)


# build corresponding lists of recs and stranges, even entries are desynched, odd are synched:
recs, stranges = [], [] # recs has repetitions, not unique
for rec in urecs: # iterate over sorted unique recs
    tranges = rec2tranges[rec]
    for trange in tranges:
        recs.append(rec)
        stranges.append(trange)
nrecsec = len(recs)
assert nrecsec == len(stranges)
assert nrecsec % 2 == 0 # should always be an even number of recording sections
fmts = ['b-', 'r-'] * (nrecsec // 2) # alternating plotting formats for desynched and synched

# get active or all neuron ids for each section of each recording:
recsecnids = get_ssnids(recs, stranges, kind=NIDSKIND)[1]

# calculate PSTHs for all sections of all recordings, collect data from each PSTH with at
# least 1 detected peak:
psthss, spiketss = [], []
psthparamsrecsec = [] # params returned for each PSTH, for each recording section
fwhmsrecsec = [] # fwhm values, for each recording section
tsrecsec = [] # peak times, for each recording section
heightsrecsec = [] # peak heights, for each recording section
sparsrecsec = [] # sparseness values of cells with at least 1 peak, for each recording section
relsrecsec = [] # reliability values of cells with at least 1 peak, for each recording section
nreplacedbynullrel = 0
for rec, nids, strange, fmt in zip(recs, recsecnids, stranges, fmts):
    psthparams = {} # various parameters for each PSTH
    psthsfwhms = [] # fwhm values of all nids in this recording section
    psthsts = [] # times of all peaks of all nids in this recording section
    psthsheights = [] # peak heights of all nids in this recording section
    n2sparseness = {} # nid:sparseness mapping for this recording section
    n2rel = {} # nid:reliability mapping for this recording section
    # psths is a regular 2D array, spikets is a 2D ragged array (list of arrays):
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK, strange=strange,
                                 plot=False, binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
    psthss.append(psths)
    spiketss.append(spikets)
    #figure(); #pl.plot(t, psths.T, '-')
    # n2count is needed for calculating reliability:
    n2count = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                             binw=TRASTERBINW, tres=TRASTERTRES)[0]
    for nid, psth, ts in zip(nids, psths, spikets):
        # calculate reliability for all nids regardless of PSTH peaks:
        cs = n2count[nid] # 2D array of spike counts over time, one row per trial
        rhos, weights = core.pairwisecorr(cs, weight=WEIGHT, invalid='ignore')
        # set rho to 0 for trial pairs with undefined rho (one or both trials with 0 spikes):
        nanis = np.isnan(rhos)
        rhos[nanis] = 0.0
        # for log plotting convenience, replace any mean rhos < NULLREL with NULLREL
        n2rel[nid] = np.mean(rhos)
        if n2rel[nid] < NULLREL:
            n2rel[nid] = NULLREL
            nreplacedbynullrel += 1
        # run PSTH peak detection:
        baseline = MEDIANX * np.median(psth)
        thresh = baseline + MINTHRESH # peak detection threshold
        print("n%d" % nid, end='')
        peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh)
        psthparams[nid] = t, psth, thresh, baseline, peakis, lis, ris
        #psthparams[nid] = get_psth_peaks(t, psth, nid)
        #t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid] # unpack
        if PLOTPSTH:
            plot_psth(psthparams, nid, fmt)
        if len(peakis) == 0:
            continue # this PSTH has no peaks, skip all subsequent measures
        fwhms = (ris - lis) * TRES * 1000 # ms
        psthsfwhms.append(fwhms)
        psthsts.append(peakis * TRES) # sec
        psthsheights.append(psth[peakis] - baseline) # peak height above baseline
        #psthsheights.append(psth[peakis]) # peak height above 0
        n2sparseness[nid] = sparseness(psth)
    psthparamsrecsec.append(psthparams)
    fwhmsrecsec.append(np.hstack(psthsfwhms))
    tsrecsec.append(np.hstack(psthsts))
    heightsrecsec.append(np.hstack(psthsheights))
    sparsrecsec.append(n2sparseness)
    relsrecsec.append(n2rel)
    print('\n') # two newlines

# 0: desynched, 1: synched
fwhms = [np.hstack(fwhmsrecsec[0::2]), # even are desynched
         np.hstack(fwhmsrecsec[1::2])] # odd are synched

ts = [np.hstack(tsrecsec[0::2]),
      np.hstack(tsrecsec[1::2])]

heights = [np.hstack(heightsrecsec[0::2]),
           np.hstack(heightsrecsec[1::2])]

spars = [np.hstack([ sparsrecsec[i].values() for i in range(0, nrecsec, 2) ]),
         np.hstack([ sparsrecsec[i].values() for i in range(1, nrecsec, 2) ])]

rels = [np.hstack([ relsrecsec[i].values() for i in range(0, nrecsec, 2) ]),
        np.hstack([ relsrecsec[i].values() for i in range(1, nrecsec, 2) ])]

'''
# plot FWHM distributions:
ticks = np.arange(FWHMMIN, FWHMMAX, FWHMTICKSTEP)
bins = np.arange(FWHMMIN, FWHMMAX+FWHMSTEP, FWHMSTEP)
figure(figsize=figsize)
n1 = hist(fwhms[1], bins=bins, color='r')[0] # synched
n0 = hist(fwhms[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=FWHMMIN, xmax=FWHMMAX)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('peak FWHM (ms)')
ylabel('PSTH peak count')
#t, p = ttest_ind(fwhms[0], fwhms[1], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(fwhms[0], fwhms[1]) # 1-sided
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f ms' % fwhms[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f ms' % fwhms[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak FWHM %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
'''

# plot FWHM distributions in log space:
logmin, logmax = log10(10), log10(FWHMMAX)
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(fwhms[1], bins=bins, color='r')[0] # synched
n0 = hist(fwhms[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**logmin, xmax=10**logmax)
ylim(ymax=n.max()+10) # effectively normalizes the histogram
#xticks(ticks)
xscale('log')
xlabel('peak FWHM (ms)')
ylabel('PSTH peak count')
#t, p = ttest_ind(log10(fwhms[0]), log10(fwhms[1]), equal_var=False) # Welch's T-test
u, p = mannwhitneyu(log10(fwhms[0]), log10(fwhms[1])) # 1-sided
# display geometric means and p value:
text(0.98, 0.98, '$\mu$ = %.1f ms' % 10**(log10(fwhms[1]).mean()), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f ms' % 10**(log10(fwhms[0]).mean()), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak FWHM log %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# plot PSTH peak time distributions:
bins = np.arange(0, TSMAX+TSSTEP, TSSTEP)
figure(figsize=figsize)
n1 = hist(ts[1], bins=bins, color='r')[0] # synched
n0 = hist(ts[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=TSMAX)
ylim(ymax=n.max()+10) # effectively normalizes the histogram
#xticks(ticks)
xlabel('PSTH peak times (sec)')
ylabel('PSTH peak count')
#t, p = ttest_ind(fwhms[0], fwhms[1], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(ts[0], ts[1]) # 1-sided
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f ms' % ts[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f ms' % ts[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak times %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

'''
# plot peak height distribution:
ticks = np.arange(HEIGHTMIN, HEIGHTMAX, HEIGHTTICKSTEP)
bins = np.arange(HEIGHTMIN, HEIGHTMAX+HEIGHTSTEP, HEIGHTSTEP)
figure(figsize=figsize)
n1 = hist(heights[1], bins=bins, color='r')[0] # synched
n0 = hist(heights[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=HEIGHTMIN, xmax=HEIGHTMAX)
ylim(ymax=n.max()+10) # effectively normalizes the histogram
xticks(ticks)
xlabel('peak amplitude (Hz)')
ylabel('PSTH peak count')
#t, p = ttest_ind(heights[0], heights[1], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(heights[0], heights[1]) # 1-sided
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f Hz' % heights[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f Hz' % heights[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak amplitude %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
'''

# plot peak height distribution in log space:
logmin, logmax = log10(2), log10(200)
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(heights[1], bins=bins, color='r')[0] # synched
n0 = hist(heights[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**0, xmax=10**logmax)
ylim(ymax=n.max()+10) # effectively normalizes the histogram
#xticks(ticks)
xscale('log')
xlabel('peak amplitude (Hz)')
ylabel('PSTH peak count')
#t, p = ttest_ind(log10(heights[0]), log10(heights[1]), equal_var=False) # Welch's T-test
u, p = mannwhitneyu(log10(heights[0]), log10(heights[1])) # 1-sided
# display geometric means and p value:
text(0.98, 0.98, '$\mu$ = %.1f Hz' % 10**(log10(heights[1]).mean()), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f Hz' % 10**(log10(heights[0]).mean()), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak amplitude log %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)


def filterdict(d, key, fallback=0):
    try:
        return d[key]
    except KeyError:
        return fallback

# Scatter plot sparseness in neighbouring synched vs desynched periods.
# Missing values (nids active in one neighbouring period but not the other)
# are assigned a sparseness of 0 to indicate they're missing.
# Get supersets of nids, one set per state transition, with each nid having at least one peak
# during at least one state on either side of the transition:
scatnids = [ np.union1d(list(sparsrecsec[recseci]), list(sparsrecsec[recseci+1]))
             for recseci in range(0, nrecsec, 2) ]
ntrans = len(scatnids) # number of state transitions
assert ntrans == nrecsec / 2 # should be half as many transitions as recording sections
scatspars = [[], []]
for trani in range(ntrans):
    nids = scatnids[trani] # superset of nids for this state transition
    for nid in nids:
        scatspars[0].append(filterdict(sparsrecsec[2*trani], nid)) # desynched
        scatspars[1].append(filterdict(sparsrecsec[2*trani+1], nid)) # synched
scatspars = np.asarray(scatspars)
figure(figsize=figsize)
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(scatspars[1], scatspars[0], 'o', mec='k', mfc='None')
nbelowsparsyxline = (scatspars[1] > scatspars[0]).sum()
nabovesparsyxline = (scatspars[0] > scatspars[1]).sum()
fractionbelowsparsyxline = nbelowsparsyxline / (nbelowsparsyxline + nabovesparsyxline)
xlabel('synchronized sparseness')
ylabel('desynchronized sparseness')
xlim(-0.02, 1)
ylim(-0.02, 1)
sparsticks = (np.arange(0, 1+0.2, 0.2), ['0', '0.2', '0.4', '0.6', '0.8', '1'])
xticks(*sparsticks)
yticks(*sparsticks)
titlestr = 'sparseness scatter %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# plot sparseness distributions:
bins = np.linspace(0, 1, NSPARSBINS)
figure(figsize=figsize)
n1 = hist(spars[1], bins=bins, histtype='step', color='r')[0] # synched
n0 = hist(spars[0], bins=bins, histtype='step', color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=1)
ylim(ymax=35)
xticks(*sparsticks)
xlabel('sparseness')
ylabel('cell count')
#t, p = ttest_ind(spars[0], spars[1], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(spars[0], spars[1]) # 1-sided
# display means and p value:
text(0.03, 0.98, '$\mu$ = %.2f' % spars[1].mean(), # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.03, 0.90, '$\mu$ = %.2f' % spars[0].mean(), # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.03, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'sparseness %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# Scatter plot reliability in neighbouring synched vs desynched periods.
# Missing values (nids active in one neighbouring period but not the other)
# are assigned a reliability of NULLREL to indicate they're missing.
# Get supersets of nids, one set per state transition, with each nid having at least one peak
# during at least one state on either side of the transition:
scatnids = [ np.union1d(list(relsrecsec[recseci]), list(relsrecsec[recseci+1]))
             for recseci in range(0, nrecsec, 2) ]
ntrans = len(scatnids) # number of state transitions
assert ntrans == nrecsec / 2 # should be half as many transitions as recording sections
scatrels = [[], []]
for trani in range(ntrans):
    nids = scatnids[trani] # superset of nids for this state transition
    for nid in nids:
        scatrels[0].append(filterdict(relsrecsec[2*trani], nid, NULLREL)) # desynched
        scatrels[1].append(filterdict(relsrecsec[2*trani+1], nid, NULLREL)) # synched
scatrels = np.asarray(scatrels)
figure(figsize=figsize)
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(scatrels[1], scatrels[0], 'o', mec='k', mfc='None')
nbelowrelsyxline = (scatrels[1] > scatrels[0]).sum()
naboverelsyxline = (scatrels[0] > scatrels[1]).sum()
fractionbelowrelsyxline = nbelowrelsyxline / (nbelowrelsyxline + naboverelsyxline)
xlabel('synchronized reliability')
ylabel('desynchronized reliability')
logmin, logmax = LOGNULLREL, 0
xscale('log')
yscale('log')
xlim(10**(logmin-0.1), 10**logmax)
ylim(10**(logmin-0.1), 10**logmax)
# replace 10^0 label with 1 to save horizontal space:
ticks = 10**(np.arange(logmin, logmax+1.0, 1.0))
ticklabels = ['$10^{%d}$' % logtick for logtick in range(logmin, 0)]
ticklabels.append('1')
xticks(ticks, ticklabels)
yticks(ticks, ticklabels)
#xticks(10**(np.arange(logmin, logmax+1.0, 1.0)))
#yticks(10**(np.arange(logmin, logmax+1.0, 1.0)))
titlestr = 'reliability scatter %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# plot reliability distributions:
bins = np.logspace(LOGNULLREL, 0, NRELBINS)
figure(figsize=figsize)
nnsynchrels = rels[1][rels[1] > NULLREL] # non-null synched reliabilities
nndesynchrels = rels[0][rels[0] > NULLREL] # non-null desynched reliabilities
n0 = hist(nnsynchrels, bins=bins, histtype='step', color='r')[0] # synched
n1 = hist(nndesynchrels, bins=bins, histtype='step', color='b')[0] # desynched
n = np.hstack([n0, n1])
#xlim(xmin=0, xmax=xmax)
#ylim(ymax=n.max()) # effectively normalizes the histogram
xscale('log')
xlim(bins[0], bins[-1])
ylim(ymax=35)
xticks(ticks, ticklabels)
#xticks(np.arange(0, xmax, 0.2))
#yticks(np.arange(0, 20, 5))
xlabel('reliability')
ylabel('cell count')
#t, p = ttest_ind(rels[1], rels[0], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(log10(nndesynchrels), log10(nnsynchrels)) # 1-sided
# display geometric means and p value:
text(0.03, 0.98, '$\mu$ = %.1e' % 10**(log10(nnsynchrels).mean()), # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.03, 0.90, '$\mu$ = %.1e' % 10**(log10(nndesynchrels).mean()), # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.03, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'reliability %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# report recordings:
print('recordings: %s' % urecnames)

# report scatter fractions:
print('fractionbelowsparsyxline: %g' % fractionbelowsparsyxline)
print('fractionbelowrelsyxline: %f' % fractionbelowrelsyxline)

# report numbers of all and active PSTHs
for kind in ['all', 'active']:
    temprecsecnids = get_ssnids(recs, stranges, kind=kind)[1]
    npsths = sum(nids.shape[0] for nids in temprecsecnids) # number that met NIDSKIND criterion
    print('%s: npsths=%d' % (kind, npsths))

# report numbers of responsive PSTHs:
nresppsths = len(spars[0]) + len(spars[1])
print('responsive PSTHs (with at least 1 peak): %d' % nresppsths)
nresppsthsdesynched = len(spars[0])
print('responsive PSTHs in desynched periods: %d' % nresppsthsdesynched)
nresppsthssynched = len(spars[1])
print('responsive PSTHs in synched periods: %d' % nresppsthssynched)
# report chi-square results of numbers of responsive PSTHs in both states:
chi2, p = chisquare([nresppsthsdesynched, nresppsthssynched])
print('chi2=%.3g, p=%.3g' % (chi2, p))

# report chi-square results of numbers of peaks in both states:
ndesynchedpeaks, nsynchedpeaks = len(fwhms[0]), len(fwhms[1]) # peak counts
chi2, p = chisquare([ndesynchedpeaks, nsynchedpeaks])
print('peak counts:')
print('ndesynched=%d, nsynched=%d, chi2=%.3g, p=%.3g'
      % (ndesynchedpeaks, nsynchedpeaks, chi2, p))

nrelsrecsec = sum([len(d) for d in relsrecsec])
print('fraction replaced by NULLREL due to rel < NULLREL: %f'
      % (nreplacedbynullrel/nrelsrecsec))

pl.show()
