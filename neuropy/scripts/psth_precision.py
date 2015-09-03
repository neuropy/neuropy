"""Detect response events within PSTHs, count them, measure the width of each one, and plot
their distributions as a function of cortical state within each of the natural scene movie
recordings listed below.

Also, measure trial raster plot reliability as a function of cortical state, by taking all
pairwise correlations between all trials for a given neuron, and then averaging them to get a
single number between 0 and 1 representing reliability of that trial raster plot. See Goard &
Dan, 2009.

Also, measure sparseness of responsive PSTHs.

Run from within neuropy using `run -i scripts/psth_precision.py`"""

from __future__ import division, print_function
from scipy.signal import argrelextrema
from scipy.stats import ttest_ind, chisquare, mannwhitneyu
from numpy import log10

from core import get_ssnids, sparseness, intround, ceilsigfig, scatterbin

from psth_funcs import plot_psth, get_psth_peaks_gac

## TODO: also take average peak widths for each cell in each state, scatter plot them vs state

## TODO: plot difference in mean response reliability between synched and desynched, as a
## function of gaussian sigma

# sort recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted
urecnames = ' '.join([rec.absname for rec in urecs])

NIDSKIND = 'all' # 'active' or 'all'

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
TRASTERBINW, TRASTERTRES = 0.02, 0.001 # trial raster bins, sec
BLANK = False # consider blank periods between trials?
WEIGHT = False # weight trials by spike count for reliability measure?
# 2.5 Hz thresh is 1 spike in the same 20 ms wide bin every 20 trials, assuming 0 baseline:
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
WIDTHMAX = 200 # maximum peak width, ms
#WIDTHMAXPOINTS = intround(WIDTHMAX / 1000 / TRES) # maximum width, number of PSTH timepoints

# plotting params:
PLOTPSTH = False
#WIDTHMIN, WIDTHSTEP, WIDTHTICKSTEP = 0, 10, 50
TSMAX, TSSTEP = 5.5, 0.25
HEIGHTMIN, HEIGHTMAX, HEIGHTSTEP, HEIGHTTICKSTEP = 0, 100, 5, 25
NSPARSBINS = 15
NRELBINS = 15
LOGNULLREL = -3
NULLREL = 10**LOGNULLREL
NULLSPARS = 0
figsize = (3, 3) # inches


# build corresponding lists of recs and stranges, even entries are desynched, odd are synched:
recs, stranges = [], [] # recs has repetitions, not unique
for rec in urecs: # iterate over sorted unique recs
    tranges = REC2STATETRANGES[rec.absname]
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
widthsrecsec = [] # peak widths, for each recording section
tsrecsec = [] # peak times, for each recording section
heightsrecsec = [] # peak heights, for each recording section
depthsrecsec = [] # physical depths of units for each peak, for each recording section
relsrecsec = [] # reliability values of cells with at least 1 peak, for each recording section
sparsrecsec = [] # sparseness values of cells with at least 1 peak, for each recording section
neurondepthsrecsec = [] # physical depths of units with at least 1 peak, for each rec section
nreplacedbynullrel = 0
for rec, nids, strange, fmt in zip(recs, recsecnids, stranges, fmts):
    print(rec.absname)
    psthparams = {} # various parameters for each PSTH
    psthswidths = [] # peak width of all nids in this recording section
    psthsts = [] # times of all peaks of all nids in this recording section
    psthsheights = [] # peak heights of all nids in this recording section
    psthsdepths = [] # physical unit depths of each peak
    n2rel = {} # nid:reliability mapping for this recording section
    n2sparseness = {} # nid:sparseness mapping for this recording section
    n2depth = {} # nid:unit depth mapping for this recording section
    # psths is a regular 2D array, spikets is a 2D ragged array (list of arrays):
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK, strange=strange,
                                 plot=False, binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
    psthss.append(psths)
    spiketss.append(spikets)
    #figure(); #pl.plot(t, psths.T, '-')
    # n2count is needed for calculating reliability:
    n2count = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                             binw=TRASTERBINW, tres=TRASTERTRES, gauss=GAUSS)[0]
    for nid, psth, ts in zip(nids, psths, spikets):
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
        npeaks = len(peakis)
        if npeaks == 0:
            continue # this PSTH has no peaks, skip all subsequent measures
        # calculate peak precision:
        widths = (ris - lis) * TRES * 1000 # ms
        psthswidths.append(widths)
        psthsts.append(peakis * TRES) # sec
        psthsheights.append(psth[peakis] - baseline) # peak height above baseline
        #psthsheights.append(psth[peakis]) # peak height above 0
        depth = rec.alln[nid].pos[1] # y position on polytrode, microns from top
        psthsdepths.append(np.tile([depth], npeaks))
        # calculate reliability of responsive PSTHs:
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
        # calculate sparseness of responsive PSTHs:
        n2sparseness[nid] = sparseness(psth)
        n2depth[nid] = depth
    psthparamsrecsec.append(psthparams)
    widthsrecsec.append(np.hstack(psthswidths))
    tsrecsec.append(np.hstack(psthsts))
    heightsrecsec.append(np.hstack(psthsheights))
    depthsrecsec.append(np.hstack(psthsdepths))
    relsrecsec.append(n2rel)
    sparsrecsec.append(n2sparseness)
    neurondepthsrecsec.append(n2depth)
    print('\n') # two newlines

# 0: desynched, 1: synched
widths = [np.hstack(widthsrecsec[0::2]), # even are desynched
          np.hstack(widthsrecsec[1::2])] # odd are synched

ts = [np.hstack(tsrecsec[0::2]),
      np.hstack(tsrecsec[1::2])]

heights = [np.hstack(heightsrecsec[0::2]),
           np.hstack(heightsrecsec[1::2])]

depths = [np.hstack(depthsrecsec[0::2]),
          np.hstack(depthsrecsec[1::2])]

rels = [[], []]
for i in range(nrecsec):
    n2r = relsrecsec[i] # nid to reliability mapping for this recsec
    rels[i%2].append( [ n2r[nid] for nid in sorted(n2r) ] ) # nid order
rels[0] = np.hstack(rels[0])
rels[1] = np.hstack(rels[1])

spars = [[], []]
for i in range(nrecsec):
    n2s = sparsrecsec[i] # nid to sparseness value mapping for this recsec
    spars[i%2].append( [ n2s[nid] for nid in sorted(n2s) ] ) # nid order
spars[0] = np.hstack(spars[0])
spars[1] = np.hstack(spars[1])

ndepths = [[], []]
for i in range(nrecsec):
    n2d = neurondepthsrecsec[i] # nid to unit depth mapping for this recsec
    ndepths[i%2].append( [ n2d[nid] for nid in sorted(n2d) ] ) # nid order
ndepths[0] = np.hstack(ndepths[0])
ndepths[1] = np.hstack(ndepths[1])


'''
# plot peak width distributions:
ticks = np.arange(WIDTHMIN, WIDTHMAX, WIDTHTICKSTEP)
bins = np.arange(WIDTHMIN, WIDTHMAX+WIDTHSTEP, WIDTHSTEP)
figure(figsize=figsize)
n1 = hist(widths[1], bins=bins, color='r')[0] # synched
n0 = hist(widths[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=WIDTHMIN, xmax=WIDTHMAX)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('event width (ms)')
ylabel('event count')
#t, p = ttest_ind(widths[0], widths[1], equal_var=False) # Welch's T-test
u, p = mannwhitneyu(widths[0], widths[1]) # 1-sided
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f ms' % widths[1].mean(), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f ms' % widths[0].mean(), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak width %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
'''

# plot peak width distributions in log space:
logmin, logmax = 0.5, log10(WIDTHMAX)
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(widths[1], bins=bins, color='r')[0] # synched
n0 = hist(widths[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**logmin, xmax=10**logmax)
ylim(ymax=n.max()+10)
#xticks(ticks)
xscale('log')
xlabel('event width (ms)')
ylabel('event count')
#t, p = ttest_ind(log10(widths[0]), log10(widths[1]), equal_var=False) # Welch's T-test
u, p = mannwhitneyu(log10(widths[0]), log10(widths[1])) # 1-sided
# display geometric means and p value:
text(0.03, 0.98, '$\mu$ = %.1f ms' % 10**(log10(widths[1]).mean()), # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.03, 0.90, '$\mu$ = %.1f ms' % 10**(log10(widths[0]).mean()), # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.03, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'peak width log %s' % urecnames
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
xlabel('event times (sec)')
ylabel('event count')
#t, p = ttest_ind(widths[0], widths[1], equal_var=False) # Welch's T-test
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


# plot peak height distribution:
ticks = np.arange(0, 25, 5)
bins = np.arange(0, 25, 1)
figure(figsize=figsize)
n1 = hist(heights[1], bins=bins, color='r')[0] # synched
n0 = hist(heights[0], bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
#xlim(xmin=HEIGHTMIN, xmax=HEIGHTMAX)
ylim(ymax=n.max()+10) # effectively normalizes the histogram
xticks(ticks)
xlabel('event amplitude (Hz)')
ylabel('event count')
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
xlabel('event amplitude (Hz)')
ylabel('event count')
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


# scatter plot peak width vs unit depth, in both cortical states:
figure(figsize=figsize)
plot(depths[0], widths[0], 'b.', ms=2) # desynched
plot(depths[1], widths[1], 'r.', ms=2) # synched
# plot trends, note that the y axis (widths) will be plotted logarithmically:
edges = np.arange(0, 1400+200, 200)
middd, logmeandw, logstddw = scatterbin(depths[0], log10(widths[0]), edges, xaverage=None)
midsd, logmeansw, logstdsw = scatterbin(depths[1], log10(widths[1]), edges, xaverage=None)
# this is tricky, but correct. Gives equal sized log error bars above and below each point:
desyncherr = [10**logmeandw-10**(logmeandw-logstddw), 10**(logmeandw+logstddw)-10**logmeandw]
syncherr = [10**logmeansw-10**(logmeansw-logstdsw), 10**(logmeansw+logstdsw)-10**logmeansw]
errorbar(middd, 10**logmeandw, yerr=desyncherr, fmt='b.-', ms=10, lw=2, zorder=999)
errorbar(midsd, 10**logmeansw, yerr=syncherr, fmt='r.-', ms=10, lw=2, zorder=999)
xlim(0, 1400)
ylim(5, 200)
xticks(np.arange(0, 1200+300, 300))
yscale('log')
xlabel('unit depth ($\mu$m)')
ylabel('event width (ms)')
titlestr = 'peak depth log %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)


def filterdict(d, key, fallback=0):
    try:
        return d[key]
    except KeyError:
        return fallback


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
scatrels = np.asarray(scatrels).T # nrows x 2 cols
figure(figsize=figsize)
truerows = (scatrels != NULLREL).all(axis=1) # exclude NULLREL rows
falserows = (scatrels == NULLREL).any(axis=1)
scatrelstrue = scatrels[truerows]
scatrelsfalse = scatrels[falserows]
# report numbers, fractions and chi2 p values for reliability scatter plot.
nbelowrelsyxline = (scatrels[:, 1] > scatrels[:, 0]).sum()
naboverelsyxline = (scatrels[:, 0] > scatrels[:, 1]).sum()
fractionbelowrelsyxline = nbelowrelsyxline / (nbelowrelsyxline + naboverelsyxline)
chi2, p = chisquare([naboverelsyxline, nbelowrelsyxline])
pstring = '$p<%g$' % ceilsigfig(p)
print('nbelowrelsyxline=%d, naboverelsyxline=%d, fractionbelowrelsyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowrelsyxline, naboverelsyxline,
                             fractionbelowrelsyxline, chi2, p))
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(scatrelstrue[:, 1], scatrelstrue[:, 0], 'o', mec='k', mfc='None')
plot(scatrelsfalse[:, 1], scatrelsfalse[:, 0], 'o', mec='e', mfc='None')
xlabel('synchronized reliability')
ylabel('desynchronized reliability')
logmin, logmax = LOGNULLREL, 0
xscale('log')
yscale('log')
xlim(10**(logmin-0.05), 10**logmax)
ylim(10**(logmin-0.05), 10**logmax)
# replace 10^0 label with 1 to save horizontal space:
ticks = 10**(np.arange(logmin, logmax+1.0, 1.0))
ticklabels = ['$10^{%d}$' % logtick for logtick in range(logmin, 0)]
ticklabels.append('1')
xticks(ticks, ticklabels)
yticks(ticks, ticklabels)
#xticks(10**(np.arange(logmin, logmax+1.0, 1.0)))
#yticks(10**(np.arange(logmin, logmax+1.0, 1.0)))
text(0.03, 0.98, '%s' % pstring, horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
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
ylim(ymax=30)
xticks(ticks, ticklabels)
#xticks(np.arange(0, xmax, 0.2))
#yticks(np.arange(0, 20, 5))
xlabel('reliability')
ylabel('unit count')
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


# scatter plot reliability vs unit depth, in both cortical states:
figure(figsize=figsize)
plot(ndepths[0], rels[0], 'b.', ms=2) # desynched
plot(ndepths[1], rels[1], 'r.', ms=2) # synched
# plot trends, note that the y axis (widths) will be plotted logarithmically:
edges = np.arange(0, 1400+200, 200)
middd, logmeandr, logstddr = scatterbin(ndepths[0], log10(rels[0]), edges, xaverage=None)
midsd, logmeansr, logstdsr = scatterbin(ndepths[1], log10(rels[1]), edges, xaverage=None)
# this is tricky, but correct. Gives equal sized log error bars above and below each point:
desyncherr = [10**logmeandr-10**(logmeandr-logstddr), 10**(logmeandr+logstddr)-10**logmeandr]
syncherr = [10**logmeansr-10**(logmeansr-logstdsr), 10**(logmeansr+logstdsr)-10**logmeansr]
errorbar(middd, 10**logmeandr, yerr=desyncherr, fmt='b.-', ms=10, lw=2, zorder=999)
errorbar(midsd, 10**logmeansr, yerr=syncherr, fmt='r.-', ms=10, lw=2, zorder=999)
xlim(0, 1400)
ylim(1e-3, 1)
xticks(np.arange(0, 1200+300, 300))
yscale('log')
xlabel('unit depth ($\mu$m)')
ylabel('reliability')
titlestr = 'reliability depth log %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)


# Scatter plot sparseness in neighbouring synched vs desynched periods.
# Missing values (nids active in one neighbouring period but not the other)
# are assigned a sparseness of NULLSPARS to indicate they're missing.
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
        scatspars[0].append(filterdict(sparsrecsec[2*trani], nid, NULLSPARS)) # desynched
        scatspars[1].append(filterdict(sparsrecsec[2*trani+1], nid, NULLSPARS)) # synched
scatspars = np.asarray(scatspars).T # nrows x 2 cols
figure(figsize=figsize)
truerows = (scatspars != NULLSPARS).all(axis=1) # exclude NULLSPARS rows
falserows = (scatspars == NULLSPARS).any(axis=1)
scatsparstrue = scatspars[truerows]
scatsparsfalse = scatspars[falserows]
# report numbers, fractions and chi2 p values for sparseness scatter plot.
nbelowsparsyxline = (scatspars[:, 1] > scatspars[:, 0]).sum()
nabovesparsyxline = (scatspars[:, 0] > scatspars[:, 1]).sum()
fractionbelowsparsyxline = nbelowsparsyxline / (nbelowsparsyxline + nabovesparsyxline)
chi2, p = chisquare([nabovesparsyxline, nbelowsparsyxline])
pstring = '$p<%g$' % ceilsigfig(p)
print('nbelowsparsyxline=%d, nabovesparsyxline=%d, fractionbelowsparsyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowsparsyxline, nabovesparsyxline,
                             fractionbelowsparsyxline, chi2, p))
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(scatsparstrue[:, 1], scatsparstrue[:, 0], 'o', mec='k', mfc='None')
plot(scatsparsfalse[:, 1], scatsparsfalse[:, 0], 'o', mec='e', mfc='None')
xlabel('synchronized sparseness')
ylabel('desynchronized sparseness')
xlim(-0.02, 1)
ylim(-0.02, 1)
sparsticks = (np.arange(0, 1+0.2, 0.2), ['0', '0.2', '0.4', '0.6', '0.8', '1'])
xticks(*sparsticks)
yticks(*sparsticks)
text(0.03, 0.98, '%s' % pstring, horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
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
ylabel('unit count')
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


# scatter plot sparseness vs unit depth, in both cortical states:
figure(figsize=figsize)
plot(ndepths[0], spars[0], 'b.', ms=2) # desynched
plot(ndepths[1], spars[1], 'r.', ms=2) # synched
edges = np.arange(0, 1400+200, 200)
middd, meands, stdds = scatterbin(ndepths[0], spars[0], edges, xaverage=None)
midsd, meanss, stdss = scatterbin(ndepths[1], spars[1], edges, xaverage=None)
errorbar(middd, meands, yerr=stdds, fmt='b.-', ms=10, lw=2, zorder=999)
errorbar(midsd, meanss, yerr=stdss, fmt='r.-', ms=10, lw=2, zorder=999)
xlim(0, 1400)
ylim(0, 1)
xticks(np.arange(0, 1200+300, 300))
xlabel('unit depth ($\mu$m)')
ylabel('sparseness')
titlestr = 'sparseness depth %s' % urecnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)


# report recordings:
print('recordings: %s' % urecnames)

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
ndesynchedpeaks, nsynchedpeaks = len(widths[0]), len(widths[1]) # peak counts
chi2, p = chisquare([ndesynchedpeaks, nsynchedpeaks])
print('peak counts:')
print('ndesynched=%d, nsynched=%d, chi2=%.3g, p=%.3g'
      % (ndesynchedpeaks, nsynchedpeaks, chi2, p))

nrelsrecsec = sum([len(d) for d in relsrecsec])
print('fraction replaced by NULLREL due to rel < NULLREL: %f'
      % (nreplacedbynullrel/nrelsrecsec))

# report total recording duration:
print()
dt = 0
for rec in urecs:
    dt += rec.dt
dtmin = dt / 1e6 / 60
dthour = dtmin / 60
print('total recording duration: %.1f min, %.1f h' % (dtmin, dthour))

# report recording duration in each state, in min:
desynchedstranges = np.asarray(stranges[0::2]) # even are desynched
synchedstranges = np.asarray(stranges[1::2]) # odd are synched
tdesynched = np.diff(desynchedstranges, axis=1).sum() / 1e6 / 60 # in min
tsynched = np.diff(synchedstranges, axis=1).sum() / 1e6 / 60 # in min
print('tsynched=%.1f min, tdesynched=%.1f min' % (tsynched, tdesynched))

print()
print('per-track unit counts:')
trackunits = {}
for rec in urecs:
    try: units = trackunits[rec.tr.absname]
    except KeyError: units = []
    trackunits[rec.tr.absname] = np.unique(np.append(units, list(rec.alln)))
total = 0
for trackname in trackunits:
    n = len(trackunits[trackname])
    total += n
    print(trackname, n)
print('total units', total)

pl.show()
