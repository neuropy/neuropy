"""Detect response events within PSTHs, count them, measure the FWHM of each one, and plot
their distributions as a function of cortical state within each of the 2 natural scene movies
in ptc22.tr1"""

from __future__ import division, print_function
from scipy.signal import argrelextrema
from scipy.stats import ttest_ind, chisquare
from numpy import log10

from core import argfwhm, get_ssnids, sparseness, intround

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

recs = ptc22tr1r08s + ptc22tr1r10s
stranges = strangesr08s + strangesr10s

EPS = np.spacing(1) # epsilon, smallest representable non-zero number

NIDSKIND = 'all' # 'active' or 'all'

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
# 2.5 Hz thresh is 1 spike in the same 20 ms wide bin every 20 trials, assuming 0 baseline:
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
FWFRACTION = 0.5 # full width fraction of max
FWHMMAX = 200 # maximum FWHM, ms
FWHMMAXPOINTS = intround(FWHMMAX / 1000 / TRES) # maximum FWHM, number of PSTH timepoints

# plotting params:
PLOTPSTH = False
FWHMMIN, FWHMSTEP, FWHMTICKSTEP = 0, 10, 50
HEIGHTMIN, HEIGHTMAX, HEIGHTSTEP, HEIGHTTICKSTEP = 0, 100, 5, 25
SPARSTEP = 0.1
figsize = (3, 3) # inches


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
        pl.plot(t[ris], psth[ris], 'bo', ms=ms, mec='none') # rigth edges
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
'''

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


## TODO: use only neurons qualitatively deemed responsive?
# get active or all neuron ids for each section of r08:
#ssnids, recsecnids = get_ssnids(ptc22tr1r08s, strangesr08s, kind=NIDSKIND)
# get active or all neuron ids for each section of both r08 and r10:
ssnids, recsecnids = get_ssnids(recs, stranges, kind=NIDSKIND)

# calculate PSTHs for both sections of both recordings:

ts, psthss = [], []
for rec, nids, strange in zip(recs, recsecnids, stranges):
    t, psths = rec.traster(nids=nids, natexps=False, strange=strange, plot=False,
                           psth=True, binw=BINW, tres=TRES, norm='ntrials')
    ts.append(t) # same time array for all PSTHs in this recording section
    psthss.append(psths)
    #figure()
    #pl.plot(t, psths.T, '-')

# collect data from each PSTH:
psthparamsrecsec = [] # params returned for each PSTH, for each recording section
fwhmsrecsec = [] # fwhm values, for each recording section
heightsrecsec = [] # peak heights, for each recording section
sparsrecsec = [] # sparseness values of cells with at least 1 peak, for each recording section
nidsrecsec = [] # nids of cells with at least 1 peak, for each recording section
psthsrecsec = [] # psths with at least 1 peak, for each recording section
for nids, t, psths, fmt in zip(recsecnids, ts, psthss, ['b-', 'r-', 'r-', 'b-']):
    psthparams = {} # params returned for each PSTH by get_psth_peaks
    psthsfwhms = [] # fwhm values from all PSTHs from this recording section
    psthsheights = [] # peak heights from all PSTHs from this recording section
    sparsenesses = [] # sparseness values for each PSTH in this recording section
    peaknids = [] # nids with peaks in this recording section
    peakpsths = [] # psths with peaks in this recording section
    for nid, psth in zip(nids, psths):
        psthparams[nid] = get_psth_peaks(t, psth, nid)
        if PLOTPSTH:
            plot_psth(psthparams, nid, fmt)
        t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid] # unpack
        if peakis == None:
            continue # this psth has no peaks, skip it
        fwhms = (ris - lis) * TRES * 1000 # sec
        psthsfwhms.append(fwhms)
        psthsheights.append(psth[peakis] - baseline) # peak height above baseline
        sparsenesses.append(sparseness(psth)) # for only those PSTHS with peaks
        peaknids.append(nid) # this nid had at least one peak
        peakpsths.append(psth)
    psthparamsrecsec.append(psthparams)
    fwhmsrecsec.append(np.hstack(psthsfwhms))
    heightsrecsec.append(np.hstack(psthsheights))
    sparsrecsec.append(np.hstack(sparsenesses))
    nidsrecsec.append(np.hstack(peaknids))
    psthsrecsec.append(np.hstack(peakpsths))
    print('\n') # two newlines

fwhms = [np.hstack([fwhmsrecsec[0], fwhmsrecsec[3]]), # desynched
         np.hstack([fwhmsrecsec[1], fwhmsrecsec[2]])] # synched

heights = [np.hstack([heightsrecsec[0], heightsrecsec[3]]), # desynched
           np.hstack([heightsrecsec[1], heightsrecsec[2]])] # synched

spars = [np.hstack([sparsrecsec[0], sparsrecsec[3]]), # desynched
         np.hstack([sparsrecsec[1], sparsrecsec[2]])] # synched

 
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
xlabel('FWHM (ms)')
ylabel('PSTH peak count')
t, p = ttest_ind(fwhms[0], fwhms[1], equal_var=False) # Welch's T-test
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f ms' % fwhms[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f ms' % fwhms[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % p,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('peak FWHM ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot FWHM distributions in log space:
logmin, logmax = log10(10), log10(FWHMMAX)
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(fwhms[1], bins=bins, color='r')[0] # synched
n0 = hist(fwhms[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**logmin, xmax=10**logmax)
ylim(ymax=n.max()) # effectively normalizes the histogram
#xticks(ticks)
xscale('log')
xlabel('FWHM (ms)')
ylabel('PSTH peak count')
t, p = ttest_ind(log10(fwhms[0]), log10(fwhms[1]), equal_var=False) # Welch's T-test
# display geometric means and p value:
text(0.98, 0.98, '$\mu$ = %.1f ms' % 10**(log10(fwhms[1]).mean()), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f ms' % 10**(log10(fwhms[0]).mean()), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % p,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('peak FWHM log ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot peak height distribution:
ticks = np.arange(HEIGHTMIN, HEIGHTMAX, HEIGHTTICKSTEP)
bins = np.arange(HEIGHTMIN, HEIGHTMAX+HEIGHTSTEP, HEIGHTSTEP)
figure(figsize=figsize)
n1 = hist(heights[1], bins=bins, color='r')[0] # synched
n0 = hist(heights[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=HEIGHTMIN, xmax=HEIGHTMAX)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('peak amplitude (Hz)')
ylabel('PSTH peak count')
t, p = ttest_ind(heights[0], heights[1], equal_var=False) # Welch's T-test
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f Hz' % heights[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f Hz' % heights[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % p,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('peak amplitude ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot peak height distribution in log space:
logmin, logmax = log10(2), log10(200)
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(heights[1], bins=bins, color='r')[0] # synched
n0 = hist(heights[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**0, xmax=10**logmax)
ylim(ymax=n.max()) # effectively normalizes the histogram
#xticks(ticks)
xscale('log')
xlabel('peak amplitude (Hz)')
ylabel('PSTH peak count')
t, p = ttest_ind(log10(heights[0]), log10(heights[1]), equal_var=False) # Welch's T-test
# display geometric means and p value:
text(0.98, 0.98, '$\mu$ = %.1f Hz' % 10**(log10(heights[1]).mean()), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f Hz' % 10**(log10(heights[0]).mean()), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % p,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('peak amplitude log ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot sparseness distributions:
ticks = np.arange(0, 1, 0.25)
bins = np.arange(0, 1+SPARSTEP, SPARSTEP)
figure(figsize=figsize)
n1 = hist(spars[1], bins=bins, color='r')[0] # synched
n0 = hist(spars[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=1)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('sparseness')
ylabel('PSTH count')
t, p = ttest_ind(spars[0], spars[1], equal_var=False) # Welch's T-test
# display means and p value:
text(0.03, 0.98, '$\mu$ = %.2f' % spars[1].mean(), # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.03, 0.90, '$\mu$ = %.2f' % spars[0].mean(), # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.03, 0.82, 'p < %.1g' % p,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('sparseness ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# report chi-square results of numbers of peaks:
ndesynched, nsynched = len(fwhms[0]), len(fwhms[1]) # peak counts
chi2, p = chisquare([ndesynched, nsynched]) # compare number of peaks in both states
print('ndesynched=%d, nsynched=%d, chi2=%.3g, p=%.3g' % (ndesynched, nsynched, chi2, p))

show()
