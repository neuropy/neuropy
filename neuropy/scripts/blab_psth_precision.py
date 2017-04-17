"""For blab animals, detect response events within PSTHs, count them, measure the width of
each one, and plot their distributions as a function of cortical state within each of the
natural scene movie recordings listed below.

Also, measure trial raster plot reliability as a function of cortical state, by taking all
pairwise correlations between all trials for a given neuron, and then averaging them to get a
single number between 0 and 1 representing reliability of that trial raster plot. See Goard &
Dan, 2009.

Also, measure sparseness of responsive PSTHs.

Run from within neuropy using `run -i scripts/blab_psth_precision.py`"""

from scipy.signal import argrelextrema
from scipy.stats import ttest_ind, chisquare, mannwhitneyu, linregress
from scipy.optimize import leastsq

from numpy import log10

import matplotlib.pyplot as plt

from core import get_ssnids, sparseness, intround, ceilsigfig, scatterbin, g

from psth_funcs import plot_psth, get_psth_peaks_gac

## TODO: also take average peak widths for each cell in each state, scatter plot them vs state

## TODO: plot difference in mean response reliability between synched and desynched, as a
## function of gaussian sigma

# sort recordings by their absname:
recs = [ eval(recname) for recname in sorted(REC2STATE2TRANGES) ] # unique, no reps, sorted
recnames = ' '.join([rec.absname for rec in recs])
states = ['d', 's'] # desynched, synched

'''
# saccade times manually read off of global motion plots, all exceeded 60 deg/sec and dropped
# back down below 60 within ~0.1 sec:
saccades = {
'ptc17.tr2b.r58': [0.45, 1.125, 1.725, 2.1, 2.43, 2.94, 3.405],
'ptc18.tr1.r38': [0.18, 0.435, 1.02, 1.605, 3.195, 3.615],
'ptc18.tr2c.r58': [0.45, 1.125, 1.725, 2.1, 2.43, 2.94, 3.405],
'ptc22.tr1.r08': [0.465, 0.75, 1.095, 1.635, 1.92, 2.205, 2.64, 2.91, 3.21, 3.51, 3.72, 3.9,
                  4.23, 4.395],
'ptc22.tr1.r10': [0.36, 0.705, 1.02, 1.59, 2.16, 2.4, 2.79, 3.885, 4.38],
'ptc22.tr4b.r49': [0.9, 1.305, 1.65, 1.95, 3.33, 3.855, 4.185]
}
'''
BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
TRASTERBINW, TRASTERTRES = 0.02, 0.001 # trial raster bins, sec
BLANK = False # consider blank periods between trials?
WEIGHT = False # weight trials by spike count for reliability measure?
# 2.5 Hz thresh is 1 spike in the same 20 ms wide bin every 20 trials, assuming 0 baseline:
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
WIDTHMAX = 300 # maximum peak width, ms
#WIDTHMAXPOINTS = intround(WIDTHMAX / 1000 / TRES) # maximum width, number of PSTH timepoints

# plotting params:
PLOTPSTH = False
OVERPLOTPSTH = False
#WIDTHMIN, WIDTHSTEP, WIDTHTICKSTEP = 0, 10, 50
TSMAX, TSSTEP = 5.5, 0.25
HEIGHTMIN, HEIGHTMAX = 0, 100
NSPARSBINS = 15
NRELBINS = 15
LOGNULLREL = -3
NULLREL = 10**LOGNULLREL
NULLSPARS = 0
figsize = (3, 3) # inches
DEPTHSRANGE = np.array([0, 1400]) # um

POOLOVEROPTO = True

psthss = {'d': [], 's': []}
spiketss = {'d': [], 's': []}

# init measurements collected across all recordings and movies, indexed by state:
allrnids = {'d': [], 's': []}
allnrnids = {'d': [], 's': []}
allrpsths = {'d': [], 's': []}
allpeaktimes = {'d': [], 's': []}
allpeakwidths = {'d': [], 's': []}
allpeakheights = {'d': [], 's': []}
allpeaknspikes = {'d': [], 's': []}
allpsthsdepths = {'d': [], 's': []}
alln2rel = []
alln2spars = []
alln2depth = []

nreplacedbynullrel = 0
ntotpeaks = 0 # num total peaks detected
# iterate over recordings:
for rec in recs:
    print(rec.absname)
    nids = sorted(rec.n)
    e = rec.e0
    # each row is trial indices for its matching movie:
    trialiss = e.p['seqnums'] - 1 # convert seqnums from 1-based to 0-based
    movienames = e.p['movie']
    umovienames = np.unique(movienames)
    # handle opto trials:
    if len(movienames) != len(umovienames):
        # some movies were displayed multiple times, probably in combination with opto stim:
        optopari, = np.where(e.p['parnames'] == 'opto')
        optovals = e.p['pars'][optopari]
        uoptovals = np.unique(optovals)
        if len(uoptovals) > 1:
            print('found multiple unique opto values: %r' % optovals)
        if POOLOVEROPTO:
            # pool trials over opto values, thereby mixing opto and non-opto trials
            # in the same raster plot:
            print('pooling over opto values')
            oldtrialiss = trialiss.copy()
            trialiss = []
            for umoviename in umovienames:
                movieis, = np.where(movienames == umoviename)
                pooledtrialis = np.sort(np.hstack(oldtrialiss[movieis]))
                trialiss.append(pooledtrialis)
    # iterate over movies:
    for trialis, moviename in zip(trialiss, umovienames):
        print('  movie: %s' % moviename)
        mvttranges = e.ttranges[trialis] # trial time ranges for this movie
        rnids = {'d': [], 's': []} # responsive neuron IDs for this movie
        nrnids = {'d': [], 's': []} # nonresponsive neuron IDs for this movie
        rpsths = {'d': [], 's': []} # responsive PSTHs for this movie
        peaktimes = {'d': [], 's': []} # peak times of all nids for this movie
        peakwidths = {'d': [], 's': []} # peak widths of all nids for this movie
        peakheights = {'d': [], 's': []} # peak heights of all nids for this movie
        peaknspikes = {'d': [], 's': []} # spike counts of each peak, normalized by ntrials
        psthsdepths = {'d': [], 's': []} # physical unit depths of each peak
        n2rel = {'d': {}, 's': {}} # nid:reliability mapping for this movie
        n2spars = {'d': {}, 's': {}} # nid:sparseness mapping for this movie
        n2depth = {'d': {}, 's': {}} # nid:unit depth mapping for this movie
        # iterate over cortical state:
        for state in states:
            ttranges = [] # build up trial tranges for this movie and state
            statetranges = REC2STATE2TRANGES[rec.absname][state]
            for statetrange in statetranges: # can be multiple time ranges for a given state
                mvttri0 = mvttranges[:, 0].searchsorted(statetrange[0])
                mvttri1 = mvttranges[:, 1].searchsorted(statetrange[1])
                ttranges.append(mvttranges[mvttri0:mvttri1])
            ttranges = np.vstack(ttranges)
            ntrials = len(ttranges)
            print('    state=%s, ntrials=%d' % (state, ntrials))
            psthparams = {} # various parameters for each PSTH
            # psths is a regular 2D array, spikets is a 2D ragged array (list of arrays):
            t, psths, spikets = rec.psth(nids=nids, ttranges=ttranges, natexps=False,
                                         blank=BLANK, strange=None, plot=False,
                                         binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
            # useful for later inspection:
            psthss[state].append(psths)
            spiketss[state].append(spikets)
            if OVERPLOTPSTH:
                figure()
                plt.plot(t, psths.T, '-') # plot all PSTHs for this movie and state
                show()
            # n2count is needed for calculating reliability:
            n2count = rec.bintraster(nids=nids, ttranges=ttranges, natexps=False,
                                     blank=BLANK, strange=None,
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
                    show()
                npeaks = len(peakis)
                if npeaks == 0:
                    nrnids[state].append(nid) # save nonresponsive nids by state
                    continue # this PSTH has no peaks, skip all subsequent measures
                ntotpeaks += npeaks
                rnids[state].append(nid) # save responsive nids by state
                rpsths[state].append(psth) # save responsive PSTH by state
                peaktimes[state].append(peakis * TRES) # save peak times, s
                # calculate peak precision:
                widths = (ris - lis) * TRES * 1000 # ms
                peakwidths[state].append(widths)
                peakheights[state].append(psth[peakis] - baseline) # peak height above baseline
                #peakheights[state].append(psth[peakis]) # peak height above 0
                lspikeis, rspikeis = ts.searchsorted(lis*TRES), ts.searchsorted(ris*TRES)
                nspikes = rspikeis - lspikeis # nspikes of each detected peak
                assert (nspikes > 0).all()
                peaknspikes[state].append(nspikes / ntrials) # nspikes per peak per trial
                depth = rec.alln[nid].pos[1] # y position on polytrode, microns from top
                psthsdepths[state].append(np.tile([depth], npeaks))
                # calculate reliability of responsive PSTHs:
                cs = n2count[nid] # 2D array of spike counts over trial time, one row per trial
                rhos, weights = core.pairwisecorr(cs, weight=WEIGHT, invalid='ignore')
                # set rho to 0 for trial pairs with undefined rho (one or both trials with
                # 0 spikes):
                nanis = np.isnan(rhos)
                rhos[nanis] = 0.0
                # for log plotting convenience, replace any mean rhos < NULLREL with NULLREL
                n2rel[state][nid] = np.mean(rhos)
                if n2rel[state][nid] < NULLREL:
                    n2rel[state][nid] = NULLREL
                    nreplacedbynullrel += 1
                # calculate sparseness of responsive PSTHs:
                n2spars[state][nid] = sparseness(psth)
                n2depth[state][nid] = depth
            print()
            # save measurements from this recording and movie:
            allrnids[state].append(rnids[state])
            allnrnids[state].append(nrnids[state])
            allrpsths[state].append(np.hstack(rpsths[state]))
            allpeaktimes[state].append(np.hstack(peaktimes[state]))
            allpeakwidths[state].append(np.hstack(peakwidths[state]))
            allpeakheights[state].append(np.hstack(peakheights[state]))
            allpeaknspikes[state].append(np.hstack(peaknspikes[state]))
            allpsthsdepths[state].append(np.hstack(psthsdepths[state]))
        alln2rel.append(n2rel)
        alln2spars.append(n2spars)
        alln2depth.append(n2depth)

# concatenate measurements across all recordings and movies
for state in states:
    allrnids[state] = np.hstack(allrnids[state])
    allnrnids[state] = np.hstack(allnrnids[state])
    allrpsths[state] = np.hstack(allrpsths[state])
    allpeaktimes[state] = np.hstack(allpeaktimes[state])
    allpeakwidths[state] = np.hstack(allpeakwidths[state])
    allpeakheights[state] = np.hstack(allpeakheights[state])
    allpeaknspikes[state] = np.hstack(allpeaknspikes[state])
    allpsthsdepths[state] = np.hstack(allpsthsdepths[state])

assert ntotpeaks == len(allpeaktimes['d']) + len(allpeaktimes['s'])
print('found %d peaks in total' % ntotpeaks)



# plot peak width distributions in log space:
logmin, logmax = 0.5, log10(WIDTHMAX)
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(allpeakwidths['s'], bins=bins, histtype='step', color='r')[0] # synched
n0 = hist(allpeakwidths['d'], bins=bins, histtype='step', color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**logmin, xmax=10**logmax)
ymax = n.max() + 25
ylim(ymax=ymax)
#xticks(ticks)
xscale('log')
xlabel('event width (ms)')
ylabel('event count')
# Welch's T-test:
#t, p = ttest_ind(log10(allpeakwidths['d']), log10(allpeakwidths['s']), equal_var=False)
u, p = mannwhitneyu(log10(allpeakwidths['d']), log10(allpeakwidths['s'])) # 1-sided
smean = 10**(log10(allpeakwidths['s']).mean()) # geometric
dmean = 10**(log10(allpeakwidths['d']).mean())
# display geometric means and p value:
text(0.03, 0.98, '$\mu$ = %.1f ms' % smean, # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.03, 0.90, '$\mu$ = %.1f ms' % dmean, # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.03, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
# arrow doesn't display correctly on log axis, use annotate instead:
annotate('', xy=(smean, (6/7)*ymax), xycoords='data', # synched
             xytext=(smean, ymax), textcoords='data',
             arrowprops=dict(fc='r', ec='none', width=1.3, headwidth=7, frac=0.5))
annotate('', xy=(dmean, (6/7)*ymax), xycoords='data', # desynched
             xytext=(dmean, ymax), textcoords='data',
             arrowprops=dict(fc='b', ec='none', width=1.3, headwidth=7, frac=0.5))
titlestr = 'peak width log %s' % recnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
show()



# Scatter plot reliability of each unit's responses in synched vs desynched for each movie.
# Missing values (nids active in one state but not the other) are assigned a reliability of
# NULLREL to indicate they're missing. For each movie, get a superset of nids, with each nid
# having at least one peak during at least one state:
scatrels = [[], []] # desynched, synched
for n2rel in alln2rel: # iterate over movies
    # get superset of nids for this movie
    nids = np.union1d(list(n2rel['d']), list(n2rel['s']))
    for nid in nids:
        scatrels[0].append(n2rel['d'].get(nid, NULLREL)) # desynched
        scatrels[1].append(n2rel['s'].get(nid, NULLREL)) # synched
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
pstring = 'p < %g' % ceilsigfig(p)
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
titlestr = 'reliability scatter %s' % recnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
show()



# Scatter plot sparseness of each unit's responses in synched vs desynched for each movie.
# Missing values (nids active in one state but not the other) are assigned a sparseness of
# NULLSPARS to indicate they're missing. For each movie, get a superset of nids, with each nid
# having at least one peak during at least one staten:
scatspars = [[], []] # desynched, synched
for n2spars in alln2spars: # iterate over movies
    # get superset of nids for this movie
    nids = np.union1d(list(n2spars['d']), list(n2spars['s']))
    for nid in nids:
        scatspars[0].append(n2spars['d'].get(nid, NULLSPARS)) # desynched
        scatspars[1].append(n2spars['s'].get(nid, NULLSPARS)) # synched
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
pstring = 'p < %g' % ceilsigfig(p)
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
titlestr = 'sparseness scatter %s' % recnames
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
show()
