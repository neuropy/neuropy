"""For blab animals, detect response events within PSTHs, count them, measure the width of
each one, and plot their distributions as a function of cortical state within each of the
natural scene movie recordings listed below.

Also, measure trial raster plot reliability as a function of cortical state, by taking all
pairwise correlations between all trials for a given neuron, and then averaging them to get a
single number between 0 and 1 representing reliability of that trial raster plot. See Goard &
Dan, 2009.

Also, measure sparseness of responsive PSTHs.

Run from within neuropy using `run -i scripts/psth_precision.py`"""

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
states = ['d', 's'] # desynch, synch

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
WIDTHMAX = 200 # maximum peak width, ms
#WIDTHMAXPOINTS = intround(WIDTHMAX / 1000 / TRES) # maximum width, number of PSTH timepoints

# plotting params:
PLOTPSTH = False
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

ts = {'d': [], 's': []}
psthss = {'d': [], 's': []}
spiketss = {'d': [], 's': []}
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
            print('    state, ntrials: %s, %d' % (state, ntrials))
            # psths is a regular 2D array, spikets is a 2D ragged array (list of arrays):
            t, psths, spikets = rec.psth(nids=nids, ttranges=ttranges, natexps=False,
                                         blank=BLANK, strange=None, plot=False,
                                         binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
            # useful for later inspection:
            ts[state].append(t)
            psthss[state].append(psths)
            spiketss[state].append(spikets)
            if PLOTPSTH:
                figure()
                plt.plot(t, psths.T, '-')
                show()
            # n2count is needed for calculating reliability:
            n2count = rec.bintraster(nids=nids, ttranges=ttranges, natexps=False,
                                     blank=BLANK, strange=None,
                                     binw=TRASTERBINW, tres=TRASTERTRES, gauss=GAUSS)[0]
