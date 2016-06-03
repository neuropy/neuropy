"""Extract optic flow vector fields from natural scene movies, average them to calculate
global motion within each movie for each specified recording. Adapted from
opencv/samples/python2/opt_flow.py. Run from within neuropy using `run -i
scripts/movie_motion_contrast_luminance.py`"""

from __future__ import division, print_function

import os

import cv2
import numpy as np

from scipy.stats import kurtosis, kurtosistest

import core
from core import intround, ceilsigfig

# sort recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted

pyr_scale = 0.5
levels = 3
winsize = 15
iterations = 3
poly_n = 5
poly_sigma = 1.2
flags = 0

FIGSIZE = (6, 3)
PLOTMOVIESIGNALS = True

MOVIEFRAMERATE = 60 # Hz, used for determining tmoviefilm

# plot example synched and desynched PSTH on movie motion plot, and point it out in
# PSTH-motion correlation plot:
EXAMPLERECNAME = 'ptc17.tr2b.r58'
EXAMPLENID = 48
EXAMPLEPSTH = [None, None]

# manbar center relative to screen center, extracted from NVS header of older recordings
# that didn't have the s.origDeg static param in their textheader:
ORIGDEGS = {'ptc17.tr2b.r58': (-3.2824440002441406, -3.054497003555298),
            'ptc18.tr1.r38': (-1.5044540166854858, -3.7839291095733643),
            'ptc18.tr2c.r58': (1.09414803981781, -0.5926640033721924)}

# calculate optic flow vector field between neighbouring pairs of frames, average their
# magnitudes to get global motion:
mot = {}
motspars = {} # sparseness of each motion signal
tmovie = {} # movie timepoints, times of all frames excluding the first
tmoviefilm = {} # movie timepoints, as filmed, times of all frames excluding the first
con = {} # contrast of each frame, excluding the first one
dcon = {} # change in contrast between successive frames
lum = {} # luminance of each frame, excluding the first one
dlum = {} # change in luminance between successive frames
for rec in urecs:
    name = rec.absname
    print(name)
    e0 = rec.e0
    movie = e0.e
    # load movie data for this recording, flips frames vertically by default
    # for bottom left origin:
    movie.load()
    degpermoviepix = e0.s.widthDeg / movie.ncellswide
    dt = e0.d.sweepSec # frame duration in seconds
    # converting all the frames from a list of arrays to a 3D numpy array is slow and
    # unnecessary. Take just the required frames before converting to an array:
    allframes = movie.frames
    frameis = np.asarray(e0.d.framei) # movie frame indices used by this recording
    frames = []
    for framei in frameis:
        frames.append(allframes[framei]) # dereference
    frames = np.asarray(frames)
    nframeintervals = len(frameis) - 1
    # calculate spatial limits of movie frames that were actually displayed:
    screenwidth = e0.I['SCREENWIDTHCM'] * e0.I['DEGPERCM'] # deg
    screenheight = e0.I['SCREENHEIGHTCM'] * e0.I['DEGPERCM']
    halfscreenwidth, halfscreenheight = screenwidth/2, screenheight/2
    moviewidth, movieheight = e0.s.widthDeg, e0.s.heightDeg # deg
    halfmoviewidth, halfmovieheight = moviewidth/2, movieheight/2
    # manbar center relative to screen center:
    try:
        xorigDeg, yorigDeg = e0.s.xorigDeg, e0.s.yorigDeg
    except AttributeError:
        xorigDeg, yorigDeg = ORIGDEGS[name]
    # movie center position, wrt screen center, deg:
    xpos, ypos = e0.d.xposDeg+xorigDeg, e0.d.yposDeg+yorigDeg
    # x and y indices into frames, spanning range of movie pixels that were on-screen,
    # assumes (x,y) origin of each movie frame is at bottom left:
    mvicenter_wrt_leftscredge = halfscreenwidth + xpos
    leftscredge_wrt_leftmviedge = halfmoviewidth - mvicenter_wrt_leftscredge
    x0i = intround(leftscredge_wrt_leftmviedge / degpermoviepix)
    x1i = intround((leftscredge_wrt_leftmviedge + screenwidth) / degpermoviepix)
    mvicenter_wrt_bottomscredge = halfscreenheight + ypos
    bottomscredge_wrt_bottommviedge = halfmovieheight - mvicenter_wrt_bottomscredge
    y0i = intround(bottomscredge_wrt_bottommviedge / degpermoviepix)
    y1i = intround((bottomscredge_wrt_bottommviedge + screenheight) / degpermoviepix)
    frames = frames[:, y0i:y1i, x0i:x1i]
    print('xis: %d:%d, yis: %d:%d' % (y0i, y1i, x0i, x1i))
    print('movie shape:', frames.shape)

    # optic flow magnitudes, in deg/sec, one per frame interval:
    mot[name] = np.zeros(nframeintervals)
    con[name] = np.zeros(nframeintervals)
    dcon[name] = np.zeros(nframeintervals)
    lum[name] = np.zeros(nframeintervals)
    dlum[name] = np.zeros(nframeintervals)
    frame0 = frames[0] # init
    for i, frame1 in enumerate(frames[1:]):
        flow = cv2.calcOpticalFlowFarneback(frame0, frame1, pyr_scale, levels, winsize,
                                            iterations, poly_n, poly_sigma, flags)
        mag, ang = cv2.cartToPolar(flow[:, :, 0], flow[:, :, 1]) # mag is in pix/frame
        # average over entire vector flow field in space, convert from pix/frame to deg/sec:
        mot[name][i] = mag.mean() * degpermoviepix / dt
        con[name][i] = frame1.std()
        dcon[name][i] = frame1.std() - frame0.std()
        lum[name][i] = frame1.mean()
        dlum[name][i] = frame1.mean() - frame0.mean()
        frame0 = frame1 # update for next iteration

    motspars[name] = core.sparseness(mot[name])
    # this doesn't measure the actual frame times, but there isn't any reason for their
    # actual display time to differ, on average, over all trials and recordings, vs how long
    # dimstim was told to display them for:
    tmovie[name] = np.arange(1, len(frames)) * dt
    tmoviefilm[name] = np.arange(1, len(frames)) / MOVIEFRAMERATE

    if not PLOTMOVIESIGNALS:
        continue

    # plot motion:
    figure(figsize=FIGSIZE)
    #plot(frameis[1:], mot[name], 'k-', lw=1.5)
    #xlabel('frame index')
    plot(tmovie[name], mot[name], 'k-', lw=1.5)
    xlabel('time (s)')
    ylabel('motion amplitude (deg/s)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_motion_%s' % name)
    tight_layout(pad=0.3)

    # plot contrast:
    figure(figsize=FIGSIZE)
    plot(tmovie[name], con[name], 'k-', lw=1.5)
    xlabel('time (s)')
    ylabel('contrast (?)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_contrast_%s' % name)
    tight_layout(pad=0.3)

    # plot delta contrast:
    figure(figsize=FIGSIZE)
    plot(tmovie[name], dcon[name], 'k-', lw=1.5)
    xlabel('time (s)')
    ylabel(r'$\Delta$ contrast (?)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_dcontrast_%s' % name)
    tight_layout(pad=0.3)

    # plot lum:
    figure(figsize=FIGSIZE)
    plot(tmovie[name], lum[name], 'k-', lw=1.5)
    xlabel('time (s)')
    ylabel('luminance (?)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_luminance_%s' % name)
    tight_layout(pad=0.3)

    # plot dlum:
    figure(figsize=FIGSIZE)
    plot(tmovie[name], dlum[name], 'k-', lw=1.5)
    xlabel('time (s)')
    ylabel(r'$\Delta$ luminance (?)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_dluminance_%s' % name)
    tight_layout(pad=0.3)


# plot motion distribution of unique movies and compare to normal distribution:
MOTIONBINW = 4 # deg/s
figure(figsize=(3, 3))
# map (movie name, framei0) tuple to motion signal (2 recs share the same movie frames):
mvi2mot = {}
for recname in sorted(mot):
    rec = eval(recname)
    mviname = os.path.basename(rec.e0.s.fname)
    framei0, framei1 = rec.e0.d.framei[0], rec.e0.d.framei[-1]
    print('%s: %s, frameis %d:%d' % (recname, mviname, framei0, framei1))
    mvi2mot[(mviname, framei0)] = mot[recname]
allmotion = np.hstack(list(mvi2mot.values()))
allmotion = np.hstack([allmotion, -allmotion]) # make it symmetric around 0
motionbins = np.arange(-300, 300+MOTIONBINW, MOTIONBINW) # deg/s, symmetric around 0
midbins = motionbins[:-1] + MOTIONBINW / 2
motioncount = np.histogram(allmotion, bins=motionbins)[0]
k = kurtosis(allmotion)
# kurtosistest() seems to use the method of Anscombe & Glynn (1983),
# http://biomet.oxfordjournals.org/content/70/1/227
z, p = kurtosistest(allmotion)
pstring = 'p < %g' % ceilsigfig(p)
# normally distributed signal with same std as data, to check that its kurtosis is 0:
#nsamples = 10000000
#normal = scipy.random.normal(0, allmotion.std(), nsamples)
#normalcount = np.histogram(normal, bins=motionbins)[0]
normalcount = core.g(0, allmotion.std(), midbins) # generate normal distrib directly
# normalize to get same probability mass:
normalcount = normalcount / normalcount.sum() * motioncount.sum()
plot(midbins, normalcount, marker=None, ls='-', c='0.7', lw=2)
plot(midbins, motioncount, marker=None, ls='-', c='k', lw=2)
text(0.98, 0.98, 'k = %.1f' % k, # kurtosis
     horizontalalignment='right', verticalalignment='top',
     transform=gca().transAxes, color='k')
text(0.98, 0.90, '%s' % pstring, # p-value of null (normal) hypothesis of kurtosis test
     horizontalalignment='right', verticalalignment='top',
     transform=gca().transAxes, color='k')
#k = kurtosis(normal)
#z, p = kurtosistest(normal)
#text(0.98, 0.82, 'k=%.1f' % k, # kurtosis
#     horizontalalignment='right', verticalalignment='top',
#     transform=gca().transAxes, color='e')
#text(0.98, 0.74, 'p=%.1g' % p, # p-value of null (normal) hypothesis of kurtosis test
#     horizontalalignment='right', verticalalignment='top',
#     transform=gca().transAxes, color='e')
xlabel('motion amplitude (deg/s)')
ylabel('frame count')
xlim(0, 300)
ylim(ymax=motioncount.max())
xticks([0, 100, 200, 300])
yticks([0, motioncount.max()])
gcfm().window.setWindowTitle('movie_global_motion_distrib')
tight_layout(pad=0.3)


"""Calculate PSTHs as in psth_precision.py, then correlate each one with its respective movie
motion, contrast and luminance signal."""

from scipy.stats import chisquare, mannwhitneyu, linregress, spearmanr

from psth_funcs import get_psth_peaks_gac

NIDSKIND = 'all' # 'active' or 'all'
BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
BLANK = False # consider blank periods between trials?
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
# correlation delay to use between stimulus and response for most plots:
CORRDELAY = 30 # ms
CORRDELAYS = range(-105, 165+15, 15) # for calculating correlation delay dependence
CORRDELAYS.remove(CORRDELAY) # remove from somewhere in the middle
# add to end of list so all subsequent variables are calculated from it
CORRDELAYS.append(CORRDELAY)
RHOMIN, RHOMAX = -1, 1
RHOBINS = np.arange(RHOMIN, RHOMAX+0.1, 0.1) # left edges + rightmost edge

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
states = ['desynch', 'synch'] * (nrecsec // 2) # alternating states

# get active or all neuron ids for each section of each recording:
recsecnids = core.get_ssnids(recs, stranges, kind=NIDSKIND)[1]

# calculate responsive PSTHs:
rts = {} # index into by [rec.absname][statei]
rpsths = {} # index into by [rec.absname][statei][nid]
for rec, nids, strange, state in zip(recs, recsecnids, stranges, states):
    print(rec.absname, state)
    statei = {'desynch': 0, 'synch': 1}[state]
    if rec.absname not in rts:
        rts[rec.absname] = [None, None] # init, index into using statei
        rpsths[rec.absname] = [{}, {}] # init, index into using statei
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK, strange=strange,
                                 plot=False, binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
    rts[rec.absname][statei] = t
    for nid, psth, ts in zip(nids, psths, spikets):
        # run PSTH peak detection:
        baseline = MEDIANX * np.median(psth)
        thresh = baseline + MINTHRESH # peak detection threshold
        peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh)
        if len(peakis) == 0: # not a responsive PSTH
            continue
        rpsths[rec.absname][statei][nid] = psth
    print('\n') # two newlines

# correlate each PSTH in each recording section with movie signals. Also, collect
# movie motion-PSTH sparseness pairs:
motrhostats, conrhostats, lumrhostats = [], [], []
for corrdelay in CORRDELAYS:
    motrhos = {'desynch': [], 'synch': []}
    conrhos = {'desynch': [], 'synch': []}
    lumrhos = {'desynch': [], 'synch': []}
    motscatspars = {'desynch': [], 'synch': []}
    NULLRHO = -1
    recsecscatmotrhos = {} # index into with [rec.absname][nid][statei], first 2 are dict keys
    recsecscatconrhos = {}
    recsecscatlumrhos = {}
    for rec, state in zip(recs, states):
        #print(rec.absname, state)
        statei = {'desynch': 0, 'synch': 1}[state]
        if rec.absname not in recsecscatmotrhos:
            recsecscatmotrhos[rec.absname] = {} # init
            recsecscatconrhos[rec.absname] = {}
            recsecscatlumrhos[rec.absname] = {}
        t = rts[rec.absname][statei]
        m, c, l = mot[rec.absname], con[rec.absname], lum[rec.absname]
        tm = tmovie[rec.absname]
        ms = motspars[rec.absname]
        psthis0d = t.searchsorted(tm) # PSTH indices closest to movie frame times, no delay
        dt = rec.e0.d.sweepSec # frame duration in seconds
        corrdelaysec = corrdelay / 1000 # convert from ms to s
        ntdelay = intround(corrdelaysec / dt) # ntimepoints to delay movie-PSTH correlation by
        if ntdelay >= 0:
            psthis = psthis0d[ntdelay:]
            movieis = np.arange(0, len(tm)-ntdelay)
        else: # ntdelay is -ve
            psthis = psthis0d[:ntdelay]
            movieis = np.arange(abs(ntdelay), len(tm))
        for nid in sorted(rpsths[rec.absname][statei]):
            psth = rpsths[rec.absname][statei][nid]
            motrho = core.corrcoef(m[movieis], psth[psthis])
            conrho = core.corrcoef(c[movieis], psth[psthis])
            lumrho = core.corrcoef(l[movieis], psth[psthis])
            motrhos[state].append(motrho)
            conrhos[state].append(conrho)
            lumrhos[state].append(lumrho)
            motscatspars[state].append([ms, core.sparseness(psth)])
            if nid not in recsecscatmotrhos[rec.absname]:
                recsecscatmotrhos[rec.absname][nid] = [NULLRHO, NULLRHO] # init list for this nid
                recsecscatconrhos[rec.absname][nid] = [NULLRHO, NULLRHO]
                recsecscatlumrhos[rec.absname][nid] = [NULLRHO, NULLRHO]
            recsecscatmotrhos[rec.absname][nid][statei] = motrho
            recsecscatconrhos[rec.absname][nid][statei] = conrho
            recsecscatlumrhos[rec.absname][nid][statei] = lumrho

            if rec.absname == EXAMPLERECNAME and nid == EXAMPLENID:
                EXAMPLEPSTH[statei] = psth[psthis0d]

    # convert to arrays:
    for state in states:
        motrhos[state] = np.asarray(motrhos[state])
        conrhos[state] = np.asarray(conrhos[state])
        lumrhos[state] = np.asarray(lumrhos[state])
        motscatspars[state] = np.asarray(motscatspars[state])
    # save mean and std of all rhos, as a function of delay:
    motrhostats.append([corrdelay, motrhos['desynch'].mean(), motrhos['desynch'].std(),
                                   motrhos['synch'].mean(), motrhos['synch'].std()])
    conrhostats.append([corrdelay, conrhos['desynch'].mean(), conrhos['desynch'].std(),
                                   conrhos['synch'].mean(), conrhos['synch'].std()])
    lumrhostats.append([corrdelay, lumrhos['desynch'].mean(), lumrhos['desynch'].std(),
                                   lumrhos['synch'].mean(), lumrhos['synch'].std()])

motrhostats = np.asarray(motrhostats)
conrhostats = np.asarray(conrhostats)
lumrhostats = np.asarray(lumrhostats)
sortis = motrhostats[:, 0].argsort() # sort by corrdelay for plotting
motrhostats = motrhostats[sortis]
conrhostats = conrhostats[sortis]
lumrhostats = lumrhostats[sortis]

# plot motion and PSTHs of example, on the same plot:
figure(figsize=FIGSIZE)
#plot(frameis[1:], mot[name], 'k-', lw=1.5)
#xlabel('frame index')
plot(tmovie[EXAMPLERECNAME], mot[EXAMPLERECNAME], 'k-', lw=1.5)
xlabel('time (s)')
ylabel('motion amplitude (deg/s)')
b = gca().twinx()
EXAMPLEPSTH[0] /= EXAMPLEPSTH[0].max() # normalize to arbitrary units
EXAMPLEPSTH[1] /= EXAMPLEPSTH[1].max()
b.plot(tmovie[EXAMPLERECNAME], EXAMPLEPSTH[0], 'b-', lw=1.5, alpha=0.5) # desynched
b.plot(tmovie[EXAMPLERECNAME], EXAMPLEPSTH[1], 'r-', lw=1.5, alpha=0.5) # synched
b.set_yticks([])
b.set_ylabel('firing rate (AU)')
mviname = os.path.basename(eval(EXAMPLERECNAME).e0.s.fname)
text(0.99, 0.98, '%s' % mviname, # movie file name
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('movie_global_motion_%s' % EXAMPLERECNAME)
tight_layout(pad=0.3)

motscatrhos, conscatrhos, lumscatrhos = [], [], []
for rec in urecs:
    motscatrhos.append(recsecscatmotrhos[rec.absname].values())
    conscatrhos.append(recsecscatconrhos[rec.absname].values())
    lumscatrhos.append(recsecscatlumrhos[rec.absname].values())

motscatrhos = np.vstack(motscatrhos)
conscatrhos = np.vstack(conscatrhos)
lumscatrhos = np.vstack(lumscatrhos)

# scatter plot motion-PSTH correlation in desynched vs synched state:
figure(figsize=(3, 3))
truerows = (motscatrhos != -1).all(axis=1) # exclude rows with -1 by collapsing across columns
falserows = (motscatrhos == -1).any(axis=1)
motscatrhostrue = motscatrhos[truerows]
motscatrhosfalse = motscatrhos[falserows]
# report numbers, fractions and chi2 p values for PSTH-motion scatter plot.
# Exclude units that were manually assigned a value of -1 in a state due to being
# nonresponsive in that state:
nbelowmotyxline = (motscatrhostrue[:, 1] > motscatrhostrue[:, 0]).sum()
nabovemotyxline = (motscatrhostrue[:, 0] > motscatrhostrue[:, 1]).sum()
fractionbelowmotyxline = nbelowmotyxline / (nbelowmotyxline + nabovemotyxline)
chi2, p = chisquare([nabovemotyxline, nbelowmotyxline])
pstring = 'p < %g' % ceilsigfig(p)
print('nbelowmotyxline=%d, nabovemotyxline=%d, fractionbelowmotyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowmotyxline, nabovemotyxline, fractionbelowmotyxline,
                             chi2, p))
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(motscatrhostrue[:, 1], motscatrhostrue[:, 0], 'o', mec='k', mfc='None')
plot(motscatrhosfalse[:, 1], motscatrhosfalse[:, 0], 'o', mec='e', mfc='None')
examplerhos = recsecscatmotrhos[EXAMPLERECNAME][EXAMPLENID]
#plot(examplerhos[1], examplerhos[0], 'o', mec='r', mfc='None')
# instead of colouring it, draw an arrow to highlight example point for CORRDELAY=0.030:
arrow(0.6, -0.35, -0.15, 0.15, head_width=0.08, head_length=0.12, length_includes_head=True,
      color='k')
text(0.02, 0.98, 'delay = %d ms' % CORRDELAY,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
text(0.02, 0.90, '%s' % pstring,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
xlabel('synchronized', color='r')
ylabel('desynchronized', color='b')
xlim(-1.035, 1)
ylim(-1.035, 1)
gcfm().window.setWindowTitle('movie_global_motion_correlation_scatter_%dms' % CORRDELAY)
tight_layout(pad=0.3)

# plot motion rho histograms:
figure(figsize=(3, 3))
dmean = motrhos['desynch'].mean()
smean = motrhos['synch'].mean()
u, p = mannwhitneyu(motrhos['desynch'], motrhos['synch']) # 1-sided
pstring = 'p < %g' % ceilsigfig(p)
nd = hist(motrhos['desynch'], bins=RHOBINS, histtype='step', color='b')[0]
ns = hist(motrhos['synch'], bins=RHOBINS, histtype='step', color='r')[0]
nmax = max(np.hstack([nd, ns]))
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
# draw arrows at means:
ah = nmax / 8 # arrow height
arrow(dmean, nmax, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='b')
arrow(smean, nmax, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='r')
xlim(xmin=-0.5, xmax=RHOMAX)
ylim(ymax=nmax*1.01)
# remove unnecessary decimal places:
rhoticks = ([-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1],
            ['-0.5', '-0.25', '0', '0.25', '0.5', '0.75', '1'])
xticks(*rhoticks)
yticks([0, nmax]) # turn off y ticks to save space
xlabel('PSTH-motion correlation')
ylabel('unit count')
text(0.98, 0.98, 'delay = %d ms' % CORRDELAY,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
text(0.98, 0.90, '$\mu$ = %.3f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.82, '$\mu$ = %.3f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.74, '%s' % pstring, color='k',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
gcfm().window.setWindowTitle('movie_global_motion_correlation_%dms' % CORRDELAY)
tight_layout(pad=0.3)

# plot rho vs motion stimulus-response delay:
figure(figsize=(3, 3))
plot(motrhostats[:, 0], motrhostats[:, 1], 'b.-') # desynch
plot(motrhostats[:, 0], motrhostats[:, 3], 'r.-') # synch
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
axhline(y=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw horizontal grey line at y=0
ymax = 0.117
ah = 0.021
arrow(30, ymax, 0, -ah, head_width=10, head_length=ah/2, length_includes_head=True,
      color='k')
xlim(-100, 130)
ylim(-0.05, ymax)
yticks(np.arange(-0.05, ymax, 0.05))
xlabel('delay (ms)')
ylabel('PSTH-motion correlation')
gcfm().window.setWindowTitle('movie_global_motion_rhovsdelay')
tight_layout(pad=0.3)

# scatter plot PSTH sparseness vs. movie motion sparseness. No significant correlation in
# either state:
figure(figsize=(3, 3))
plot(motscatspars['synch'][:, 0], motscatspars['synch'][:, 1], 'r.', ms=2)
plot(motscatspars['desynch'][:, 0], motscatspars['desynch'][:, 1], 'b.', ms=2)
xlim(0, 1)
ylim(0, 1)
r, p = spearmanr(motscatspars['synch'])
pstring = 'p < %g' % ceilsigfig(p)
text(0.02, 0.01, '$r_s$=%.2f, %s' % (r, pstring), color='r', alpha=0.7,
     transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
r, p = spearmanr(motscatspars['desynch'])
pstring = 'p < %g' % ceilsigfig(p)
text(0.02, 0.09, '$r_s$=%.2f, %s' % (r, pstring), color='b', alpha=0.7,
     transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
#xr = np.asarray(xlim()) # x range
#m, b, r, p, stderr = linregress(motscatspars['synch'])
#plot(xr, m*xr+b, 'r', ls='--', lw=2, alpha=0.7)
#text(0.02, 0.02, 'r=%.2f, p=%.1g' % (r, p), color='r',
     #transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
#m, b, r, p, stderr = linregress(motscatspars['desynch'])
#plot(xr, m*xr+b, 'b', ls='--', lw=2, alpha=0.7)
#text(0.02, 0.1, 'r=%.2f, p=%.1g' % (r, p), color='b',
     #transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
#allmotscatspars = np.vstack(motscatspars.values())
#m, b, r, p, stderr = linregress(allmotscatspars)
#plot(xr, m*xr+b, 'k', ls='--', lw=2, alpha=0.7)
#text(0.02, 0.18, 'r=%.2f, p=%.1g' % (r, p), color='k',
     #transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
xlabel('movie motion sparseness')
ylabel('PSTH sparseness')
gcfm().window.setWindowTitle('movie_PSTH_sparseness')
tight_layout(pad=0.3)

# scatter plot contrast-PSTH correlation in desynched vs synched state:
figure(figsize=(3, 3))
truerows = (conscatrhos != -1).all(axis=1) # exclude rows with -1 by collapsing across columns
falserows = (conscatrhos == -1).any(axis=1)
conscatrhostrue = conscatrhos[truerows]
conscatrhosfalse = conscatrhos[falserows]
# report numbers, fractions and chi2 p values for PSTH-contrast scatter plot.
# Exclude units that were manually assigned a value of -1 in a state due to being
# nonresponsive in that state:
nbelowconyxline = (conscatrhostrue[:, 1] > conscatrhostrue[:, 0]).sum()
naboveconyxline = (conscatrhostrue[:, 0] > conscatrhostrue[:, 1]).sum()
fractionbelowconyxline = nbelowconyxline / (nbelowconyxline + naboveconyxline)
chi2, p = chisquare([naboveconyxline, nbelowconyxline])
pstring = 'p < %g' % ceilsigfig(p)
print('nbelowconyxline=%d, naboveconyxline=%d, fractionbelowconyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowconyxline, naboveconyxline, fractionbelowconyxline,
                             chi2, p))
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(conscatrhostrue[:, 1], conscatrhostrue[:, 0], 'o', mec='k', mfc='None')
plot(conscatrhosfalse[:, 1], conscatrhosfalse[:, 0], 'o', mec='e', mfc='None')
text(0.02, 0.98, 'delay = %d ms' % CORRDELAY,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
text(0.02, 0.90, '%s' % pstring,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
xlabel('synchronized', color='r')
ylabel('desynchronized', color='b')
xlim(-1.035, 1)
ylim(-1.035, 1)
gcfm().window.setWindowTitle('movie_global_contrast_correlation_scatter_%dms' % CORRDELAY)
tight_layout(pad=0.3)

# plot contrast rho histograms:
figure(figsize=(3, 3))
dmean = conrhos['desynch'].mean()
smean = conrhos['synch'].mean()
u, p = mannwhitneyu(conrhos['desynch'], conrhos['synch']) # 1-sided
pstring = 'p < %g' % ceilsigfig(p)
nd = hist(conrhos['desynch'], bins=RHOBINS, histtype='step', color='b')[0]
ns = hist(conrhos['synch'], bins=RHOBINS, histtype='step', color='r')[0]
nmax = max(np.hstack([nd, ns]))
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
# draw arrows at means:
ah = nmax / 8 # arrow height
ymax = nmax * 1.15
arrow(dmean, nmax*1.05, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='b')
arrow(smean, nmax*1.05, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='r')
xlim(xmin=RHOMIN, xmax=RHOMAX)
ylim(ymax=ymax)
# remove unnecessary decimal places:
rhoticks = ([-1, -0.5, 0, 0.5, 1],
            ['-1', '-0.5', '0', '0.5', '1'])
xticks(*rhoticks)
yticks([0, nmax]) # turn off y ticks to save space
xlabel('PSTH-contrast correlation')
ylabel('unit count')
text(0.98, 0.98, 'delay = %d ms' % CORRDELAY,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
text(0.98, 0.90, '$\mu$ = %.3f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, '$\mu$ = %.3f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.74, '%s' % pstring, color='k',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
gcfm().window.setWindowTitle('movie_global_contrast_correlation_%dms' % CORRDELAY)
tight_layout(pad=0.3)

# plot rho vs contrast stimulus-response delay:
figure(figsize=(3, 3))
plot(conrhostats[:, 0], conrhostats[:, 1], 'b.-') # desynch
plot(conrhostats[:, 0], conrhostats[:, 3], 'r.-') # synch
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
axhline(y=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw horizontal grey line at y=0
ymin, ymax = -0.115, 0.04
ah = 0.02
arrow(45, ymin, 0, ah, head_width=10, head_length=ah/2, length_includes_head=True,
      color='k')
xlim(-70, 160)
ylim(ymin, ymax)
yticks(np.arange(-0.1, ymax+0.05, 0.05))
xlabel('delay (ms)')
ylabel('PSTH-contrast correlation')
gcfm().window.setWindowTitle('movie_global_contrast_rhovsdelay')
tight_layout(pad=0.3)

# scatter plot luminance-PSTH correlation in desynched vs synched state:
figure(figsize=(3, 3))
truerows = (lumscatrhos != -1).all(axis=1) # exclude rows with -1 by collapsing across columns
falserows = (lumscatrhos == -1).any(axis=1)
lumscatrhostrue = lumscatrhos[truerows]
lumscatrhosfalse = lumscatrhos[falserows]
# report numbers, fractions and chi2 p values for PSTH-luminance scatter plot.
# Exclude units that were manually assigned a value of -1 in a state due to being
# nonresponsive in that state:
nbelowlumyxline = (lumscatrhostrue[:, 1] > lumscatrhostrue[:, 0]).sum()
nabovelumyxline = (lumscatrhostrue[:, 0] > lumscatrhostrue[:, 1]).sum()
fractionbelowlumyxline = nbelowlumyxline / (nbelowlumyxline + nabovelumyxline)
chi2, p = chisquare([nabovelumyxline, nbelowlumyxline])
pstring = 'p < %g' % ceilsigfig(p)
print('nbelowlumyxline=%d, nabovelumyxline=%d, fractionbelowlumyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowlumyxline, nabovelumyxline, fractionbelowlumyxline,
                             chi2, p))
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(lumscatrhostrue[:, 1], lumscatrhostrue[:, 0], 'o', mec='k', mfc='None')
plot(lumscatrhosfalse[:, 1], lumscatrhosfalse[:, 0], 'o', mec='e', mfc='None')
text(0.02, 0.98, 'delay = %d ms' % CORRDELAY,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
text(0.02, 0.90, '%s' % pstring,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
xlabel('synchronized', color='r')
ylabel('desynchronized', color='b')
xlim(-1.035, 1)
ylim(-1.035, 1)
gcfm().window.setWindowTitle('movie_global_luminance_correlation_scatter_%dms' % CORRDELAY)
tight_layout(pad=0.3)

# plot luminance rho histograms:
figure(figsize=(3, 3))
dmean = lumrhos['desynch'].mean()
smean = lumrhos['synch'].mean()
u, p = mannwhitneyu(lumrhos['desynch'], lumrhos['synch']) # 1-sided
pstring = 'p < %g' % ceilsigfig(p)
nd = hist(lumrhos['desynch'], bins=RHOBINS, histtype='step', color='b')[0]
ns = hist(lumrhos['synch'], bins=RHOBINS, histtype='step', color='r')[0]
nmax = max(np.hstack([nd, ns]))
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
# draw arrows at means:
ah = nmax / 8 # arrow height
ymax = nmax * 1.15
arrow(dmean, nmax*1.05, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='b')
arrow(smean, nmax*1.05, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='r')
xlim(xmin=RHOMIN, xmax=RHOMAX)
ylim(ymax=ymax)
# remove unnecessary decimal places:
rhoticks = ([-1, -0.5, 0, 0.5, 1],
            ['-1', '-0.5', '0', '0.5', '1'])
xticks(*rhoticks)
yticks([0, nmax]) # turn off y ticks to save space
xlabel('PSTH-luminance correlation')
ylabel('unit count')
text(0.98, 0.98, 'delay = %d ms' % CORRDELAY,
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
text(0.98, 0.90, '$\mu$ = %.3f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, '$\mu$ = %.3f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.74, '%s' % pstring, color='k',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
gcfm().window.setWindowTitle('movie_global_luminance_correlation_%dms' % CORRDELAY)
tight_layout(pad=0.3)

# plot rho vs luminance stimulus-response delay:
figure(figsize=(3, 3))
plot(lumrhostats[:, 0], lumrhostats[:, 1], 'b.-') # desynch
plot(lumrhostats[:, 0], lumrhostats[:, 3], 'r.-') # synch
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
axhline(y=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw horizontal grey line at y=0
ymax = 0.02
ah = 0.01
arrow(45, ymax, 0, -ah, head_width=10, head_length=ah/2, length_includes_head=True,
      color='k')
xlim(-70, 160)
ylim(-0.06, ymax)
yticks(np.arange(-0.06, ymax+0.02, 0.02))
xlabel('delay (ms)')
ylabel('PSTH-luminance correlation')
gcfm().window.setWindowTitle('movie_global_luminance_rhovsdelay')
tight_layout(pad=0.3)

pl.show()


# save motion traces to .mat files for use in MATLAB in blab:
from scipy.io import savemat
MVI_1400_200_500 = np.asarray([tmoviefilm['ptc17.tr2b.r58'], mot['ptc17.tr2b.r58']])
MVI_1400_3300_3600 = np.asarray([tmoviefilm['ptc22.tr1.r10'], mot['ptc22.tr1.r10']])
MVI_1403_0_300 = np.asarray([tmoviefilm['ptc22.tr1.r08'], mot['ptc22.tr1.r08']])
# use blab names as keys:
motdict = {'MAS_1400_200_500': MVI_1400_200_500,
           'MAS_1400_B': MVI_1400_3300_3600,
           'MAS_1403_0_300': MVI_1403_0_300}
savemat('scripts/MAS_movie_motion.mat', motdict)
