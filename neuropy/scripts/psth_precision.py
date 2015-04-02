"""Detect response events within PSTHs, count them, measure the FWHM of each one, and plot
their distributions as a function of cortical state within each of the 2 natural scene movies
in ptc22.tr1.

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
TRASTERBINW, TRASTERTRES = 0.02, 0.005 # trial raster bins, sec
BLANK = False # consider blank periods between trials for reliability measure?
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
'''
# copied to psth_precision_inactive.py:
def get_psth_peaks(t, psth, nid):
    """Extract peaks from PSTH, simpler, faster, more robust method. Find contiguous ranges of
    baseline-exceeding points. Within each range, find the biggest value. If that value
    exceeds thresh, designate that as a peak. Slice out that baseline-exceeding range of data,
    and run argfwhm on it, with outer kwarg - ie search from the outer edges in when looking
    for FWHM. If baseline is so high that it exceeds FWHM for that peak, discard the peak.
    t is passed only so that it can be conveniently returned for plotting"""
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


# get active or all neuron ids for each section of both r08 and r10:
ssnids, recsecnids = get_ssnids(recs, stranges, kind=NIDSKIND)

# calculate PSTHs for both sections of both recordings, collect data from each PSTH with at
# least 1 detected peak:
psthss = []
psthparamsrecsec = [] # params returned for each PSTH, for each recording section
fwhmsrecsec = [] # fwhm values, for each recording section
tsrecsec = [] # peak times, for each recording section
heightsrecsec = [] # peak heights, for each recording section
sparsrecsec = [] # sparseness values of cells with at least 1 peak, for each recording section
relsrecsec = [] # reliability values of cells with at least 1 peak, for each recording section
fmts = ['b-', 'r-', 'r-', 'b-']
nreplacedbynullrel = 0
for rec, nids, strange, fmt in zip(recs, recsecnids, stranges, fmts):
    psthparams = {} # params returned by get_psth_peaks of all nids in this recording section
    psthsfwhms = [] # fwhm values of all nids in this recording section
    psthsts = [] # times of all peaks of all nids in this recording section
    psthsheights = [] # peak heights of all nids in this recording section
    n2sparseness = {} # nid:sparseness mapping for this recording section
    n2rel = {} # nid:reliability mapping for this recording section
    t, psths = rec.psth(nids=nids, natexps=False, strange=strange, plot=False,
                        binw=BINW, tres=TRES, norm='ntrials')
    psthss.append(psths)
    #figure(); #pl.plot(t, psths.T, '-')
    n2count = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                             binw=TRASTERBINW, tres=TRASTERTRES)[0]
    for nid, psth in zip(nids, psths):
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
        psthparams[nid] = get_psth_peaks(t, psth, nid)
        if PLOTPSTH:
            plot_psth(psthparams, nid, fmt)
        t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid] # unpack
        if peakis == None:
            continue # this PSTH has no peaks, skip all subsequent measures
        fwhms = (ris - lis) * TRES * 1000 # ms
        psthsfwhms.append(fwhms)
        psthsts.append(peakis * TRES) # sec
        psthsheights.append(psth[peakis] - baseline) # peak height above baseline
        n2sparseness[nid] = sparseness(psth)
    psthparamsrecsec.append(psthparams)
    fwhmsrecsec.append(np.hstack(psthsfwhms))
    tsrecsec.append(np.hstack(psthsts))
    heightsrecsec.append(np.hstack(psthsheights))
    sparsrecsec.append(n2sparseness)
    relsrecsec.append(n2rel)
    print('\n') # two newlines

fwhms = [np.hstack([fwhmsrecsec[0], fwhmsrecsec[3]]), # desynched
         np.hstack([fwhmsrecsec[1], fwhmsrecsec[2]])] # synched

ts = [np.hstack([tsrecsec[0], tsrecsec[3]]), # desynched
      np.hstack([tsrecsec[1], tsrecsec[2]])] # synched

heights = [np.hstack([heightsrecsec[0], heightsrecsec[3]]), # desynched
           np.hstack([heightsrecsec[1], heightsrecsec[2]])] # synched

spars = [np.hstack([sparsrecsec[0].values(), sparsrecsec[3].values()]), # desynched
         np.hstack([sparsrecsec[1].values(), sparsrecsec[2].values()])] # synched

rels = [np.hstack([relsrecsec[0].values(), relsrecsec[3].values()]), # desynched
        np.hstack([relsrecsec[1].values(), relsrecsec[2].values()])] # synched
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
gcfm().window.setWindowTitle('peak FWHM ptc22.tr1.r08 ptc22.tr1.r10')
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
ylim(ymax=n.max()) # effectively normalizes the histogram
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
gcfm().window.setWindowTitle('peak FWHM log ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot PSTH peak time distributions:
bins = np.arange(0, TSMAX+TSSTEP, TSSTEP)
figure(figsize=figsize)
n1 = hist(ts[1], bins=bins, color='r')[0] # synched
n0 = hist(ts[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=TSMAX)
ylim(ymax=n.max()) # effectively normalizes the histogram
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
ylim(ymax=n.max()) # effectively normalizes the histogram
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
gcfm().window.setWindowTitle('peak amplitude ptc22.tr1.r08 ptc22.tr1.r10')
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
ylim(ymax=n.max()) # effectively normalizes the histogram
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
gcfm().window.setWindowTitle('peak amplitude log ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)


def filterdict(d, key, fallback=0):
    try:
        return d[key]
    except KeyError:
        return fallback

# nids with at least one peak during at least one state on either side of 1st transition:
nids0 = np.union1d(sorted(sparsrecsec[0]), sorted(sparsrecsec[1]))
# nids with at least one peak during at least one state on either side of 2nd transition:
nids1 = np.union1d(sorted(sparsrecsec[2]), sorted(sparsrecsec[3]))

# scatter plot sparseness in neighbouring synched vs desynched periods.
# Missing values (nids active in one neighbouring period but not the other)
# are assigned a sparseness of 0 to indicate they're missing:
# 1st transition: desynched to synched
desynchspars0 = [ filterdict(sparsrecsec[0], nid) for nid in nids0 ]
synchspars0 = [ filterdict(sparsrecsec[1], nid) for nid in nids0 ]
# 2nd transition: synched to desynched
synchspars1 = [ filterdict(sparsrecsec[2], nid) for nid in nids1 ]
desynchspars1 = [ filterdict(sparsrecsec[3], nid) for nid in nids1 ]
# combine all transitions into one plot:
synchspars = np.asarray(synchspars0 + synchspars1)
desynchspars = np.asarray(desynchspars0 + desynchspars1)
figure(figsize=figsize)
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(synchspars, desynchspars, 'o', mec='k', mfc='None')
nbelowsparsyxline = (synchspars > desynchspars).sum()
nabovesparsyxline = (desynchspars > synchspars).sum()
fractionbelowsparsyxline = nbelowsparsyxline / (nbelowsparsyxline + nabovesparsyxline)
print('fractionbelowsparsyxline: %g' % fractionbelowsparsyxline)
xlabel('synchronized sparseness')
ylabel('desynchronized sparseness')
xlim(-0.02, 1)
ylim(-0.02, 1)
sparsticks = (np.arange(0, 1+0.2, 0.2), ['0', '0.2', '0.4', '0.6', '0.8', '1'])
xticks(*sparsticks)
yticks(*sparsticks)
gcfm().window.setWindowTitle('sparseness scatter ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot sparseness distributions:
bins = np.linspace(0, 1, NSPARSBINS)
figure(figsize=figsize)
n1 = hist(spars[1], bins=bins, histtype='step', color='r')[0] # synched
n0 = hist(spars[0], bins=bins, histtype='step', color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=1)
ylim(ymax=20) # ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(*sparsticks)
xlabel('sparseness')
ylabel('cell count')
#t, p = ttest_ind(spars[0], spars[1], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(spars[0], spars[1]) # 1-sided
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.2f' % spars[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.2f' % spars[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('sparseness ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)


# nids active in at least one state on either side of 1st transition:
nids0 = np.union1d(sorted(relsrecsec[0]), sorted(relsrecsec[1]))
# nids active in at least one state on either side of 2nd transition:
nids1 = np.union1d(sorted(relsrecsec[2]), sorted(relsrecsec[3]))

# scatter plot reliability in neighbouring synched vs desynched periods.
# Missing values (nids active in one neighbouring period but not the other)
# are assigned a reliability of NULLREL to indicate they're missing:
# 1st transition: desynched to synched
desynchrels0 = [ filterdict(relsrecsec[0], nid, NULLREL) for nid in nids0 ]
synchrels0 = [ filterdict(relsrecsec[1], nid, NULLREL) for nid in nids0 ]
# 2nd transition: synched to desynched
synchrels1 = [ filterdict(relsrecsec[2], nid, NULLREL) for nid in nids1 ]
desynchrels1 = [ filterdict(relsrecsec[3], nid, NULLREL) for nid in nids1 ]
# combine all transitions into one plot:
synchrels = np.asarray(synchrels0 + synchrels1)
desynchrels = np.asarray(desynchrels0 + desynchrels1)
figure(figsize=figsize)
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(synchrels, desynchrels, 'o', mec='k', mfc='None')
nbelowrelsyxline = (synchrels > desynchrels).sum()
naboverelsyxline = (desynchrels > synchrels).sum()
fractionbelowrelsyxline = nbelowrelsyxline / (nbelowrelsyxline + naboverelsyxline)
print('fractionbelowrelsyxline: %f' % fractionbelowrelsyxline)
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
titlestr = 'reliability scatter ptc22.tr1.r08 ptc22.tr1.r10'
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
ylim(ymax=20)
xticks(ticks, ticklabels)
#xticks(np.arange(0, xmax, 0.2))
#yticks(np.arange(0, 20, 5))
xlabel('reliability')
ylabel('cell count')
#t, p = ttest_ind(rels[1], rels[0], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(log10(nndesynchrels), log10(nnsynchrels)) # 1-sided
# display geometric means and p value:
text(0.98, 0.98, '$\mu$ = %.1e' % 10**(log10(nnsynchrels).mean()), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1e' % 10**(log10(nndesynchrels).mean()), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'reliability ptc22.tr1.r08 ptc22.tr1.r10'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)


# report numbers of all and active PSTHs
for kind in ['all', 'active']:
    tempssnids, temprecsecnids = get_ssnids(recs, stranges, kind=kind)
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
