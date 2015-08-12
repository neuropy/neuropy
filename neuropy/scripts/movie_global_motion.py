"""Extract optic flow vector fields from natural scene movies, average them to calculate
global motion within each movie for each specified recording. Adapted from
opencv/samples/python2/opt_flow.py. Run from within neuropy using `run -i
scripts/movie_global_motion.py`"""

from __future__ import division, print_function

import cv2
import numpy as np

from scipy.stats import kurtosis, kurtosistest

# mapping of recording to list of desynched and synched trange, in that order, copied from
# psth_precision.py:
rec2tranges = {ptc17.tr2b.r58: [(0, 700e6), # desynched trange, 66 Hz refresh rate
                                (800e6, 1117e6)], # synched trange, 66 Hz refresh rate
               ptc18.tr1.r38:  [(0, 425e6), # desynched trange, ends ~ trial 76
                                (550e6, 2243e6)], # synched trange, starts ~ trial 98
               ptc18.tr2c.r58: [(0, 750e6), # desynched trange
                                (1000e6, 2248e6)], # synched trange
               ptc22.tr1.r08:  [(0, 1500e6), # desynched trange
                                (1550e6, 2329e6)], # synched trange
               ptc22.tr1.r10:  [(1480e6, 2331e6), # desynched trange
                                (0, 1400e6)], # synched trange
               ptc22.tr4b.r49: [(0, 1475e6), # desynched trange
                                (1500e6, 2331e6)], # synched trange
              }
# compare and sort recordings by their absname:
reccmp = lambda reca, recb: cmp(reca.absname, recb.absname)
urecs = sorted(rec2tranges, cmp=reccmp) # unique recordings, no repetition, sorted

pyr_scale = 0.5
levels = 3
winsize = 15
iterations = 3
poly_n = 5
poly_sigma = 1.2
flags = 0

FIGSIZE = (6, 3)

# calculate optic flow vector field between neighbouring pairs of frames, average their
# magnitudes to get global motion:
motion, tmotion, motionspars = {}, {}, {}
for rec in urecs:
    name = rec.absname
    print(name)
    motion[name] = [] # optic flow magnitudes, in deg/sec, one per frame interval
    e0 = rec.e0
    movie = e0.e
    movie.load() # load movie data for this recording, flips frames vertically by default
    degpermoviepix = e0.s.widthDeg / movie.ncellswide
    dt = e0.d.sweepSec # frame duration in seconds
    frames = np.asarray(e0.e.frames)
    frameis = np.asarray(e0.d.framei) # movie frame indices used by this recording
    frames = frames[frameis] # dereference
    frame0 = frames[0] # init
    for frame1 in frames[1:]:
        ## TODO: if interested in flow direction, double-check vertical order of movie frames
        ## vs. what's expected by cv2.calcOpticalFlowFarneback. Should frames be flipped?
        flow = cv2.calcOpticalFlowFarneback(frame0, frame1, pyr_scale, levels, winsize,
                                            iterations, poly_n, poly_sigma, flags)
        mag, ang = cv2.cartToPolar(flow[:, :, 0], flow[:, :, 1]) # mag is in pix/frame
        # average over entire vector flow field in space, convert from pix/frame to deg/sec:
        motion[name].append(mag.mean() * degpermoviepix / dt)
        frame0 = frame1 # update for next iteration
    motion[name] = np.asarray(motion[name])
    # this doesn't measure the actual frame times, but there isn't any reason for their actual
    # display time to differ, on average, over all trials and recordings, vs how long dimstim
    # was told to display them for:
    tmotion[name] = np.arange(1, len(frames)) * dt
    motionspars[name] = core.sparseness(motion[name])
    figure(figsize=FIGSIZE)
    #plot(frameis[1:], motion[name], 'k-', lw=1.5)
    #xlabel('frame index')
    plot(tmotion[name], motion[name], 'k-', lw=1.5)
    xlabel('t (s)')
    ylabel('motion amplitude (deg/s)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_motion_%s' % name)
    tight_layout(pad=0.3)


# plot motion distribution and compare to normal distribution:
MOTIONBINW = 4 # deg/s
figure(figsize=(3, 3))
allmotion = np.hstack(list(motion.values()))
allmotion = np.hstack([allmotion, -allmotion]) # make it symmetric around 0
motionbins = np.arange(-300, 300+MOTIONBINW, MOTIONBINW) # deg/s, symmetric around 0
midbins = motionbins[:-1] + MOTIONBINW / 2
motioncount = np.histogram(allmotion, bins=motionbins)[0]
k = kurtosis(allmotion)
# kurtosistest() seems to use the method of Anscombe & Glynn (1983),
# http://biomet.oxfordjournals.org/content/70/1/227
z, p = kurtosistest(allmotion)
# normally distributed signal with same std as data, to check that its kurtosis is 0:
#nsamples = 10000000
#normal = scipy.random.normal(0, allmotion.std(), nsamples)
#normalcount = np.histogram(normal, bins=motionbins)[0]
normalcount = core.g(0, allmotion.std(), midbins) # generate normal distrib directly
# normalize to get same probability mass:
normalcount = normalcount / normalcount.sum() * motioncount.sum()
plot(midbins, normalcount, marker=None, ls='-', c='0.7', lw=2)
plot(midbins, motioncount, marker=None, ls='-', c='k', lw=2)
text(0.98, 0.98, 'k=%.1f' % k, # kurtosis
     horizontalalignment='right', verticalalignment='top',
     transform=gca().transAxes, color='k')
text(0.98, 0.90, 'p=%.1g' % p, # p-value of null (normal) hypothesis of kurtosis test
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
ylim(ymax=count.max())
xticks([0, 100, 200, 300])
yticks([0, count.max()])
gcfm().window.setWindowTitle('movie_global_motion_distrib')
tight_layout(pad=0.3)


"""Calculate PSTHs as in psth_precision.py, then correlate each one with its respective movie
motion signal."""

from scipy.stats import chisquare, mannwhitneyu, linregress

import core
from core import intround, ceilsigfig
from psth_funcs import get_psth_peaks_gac

NIDSKIND = 'all' # 'active' or 'all'
BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
BLANK = False # consider blank periods between trials?
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
CORRDELAY = 0.030 # correlation delay to use between stimulus and response, sec
CORRDELAYMS = intround(CORRDELAY*1000)
RHOMIN, RHOMAX = -0.5, 1

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
states = ['desynch', 'synch'] * (nrecsec // 2) # alternating states

# get active or all neuron ids for each section of each recording:
recsecnids = core.get_ssnids(recs, stranges, kind=NIDSKIND)[1]

# correlate each PSTH in each recording section with the movie motion signal. Also, collect
# movie motion-PSTH sparseness pairs:
rhos = {'desynch': [], 'synch': []}
spars = {'desynch': [], 'synch': []}
NULLRHO = -1
recsecscatrhos = {} # index into with [rec.absname][nid][statei], first 2 are dict keys
for rec, nids, strange, state in zip(recs, recsecnids, stranges, states):
    print(rec.absname, state)
    if rec.absname not in recsecscatrhos:
        recsecscatrhos[rec.absname] = {} # init
    if state == 'desynch':
        statei = 0
    else: # state == 'synch'
        statei = 1
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK, strange=strange,
                                 plot=False, binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
    m = motion[rec.absname]
    tm = tmotion[rec.absname]
    ms = motionspars[rec.absname]
    psthis = t.searchsorted(tm) # PSTH indices closest to movie frame times
    dt = rec.e0.d.sweepSec # frame duration in seconds
    ntdelay = intround(CORRDELAY / dt) # num timepoints to delay movie-PSTH correlation by
    if ntdelay >= 0:
        psthis = psthis[ntdelay:]
        motionis = np.arange(0, len(tm)-ntdelay)
    else: # ntdelay is -ve
        psthis = psthis[:ntdelay]
        motionis = np.arange(abs(ntdelay), len(tm))
    for nid, psth, ts in zip(nids, psths, spikets):
        # run PSTH peak detection:
        baseline = MEDIANX * np.median(psth)
        thresh = baseline + MINTHRESH # peak detection threshold
        peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh)
        if len(peakis) == 0: # not a responsive PSTH
            continue
        rho = core.corrcoef(m[motionis], psth[psthis])
        rhos[state].append(rho)
        spars[state].append([ms, core.sparseness(psth)])
        if nid not in recsecscatrhos[rec.absname]:
            recsecscatrhos[rec.absname][nid] = [NULLRHO, NULLRHO] # init list for this nid
        recsecscatrhos[rec.absname][nid][statei] = rho

    print('\n') # two newlines

scatrhos = []
for rec in urecs:
    scatrhos.append(recsecscatrhos[rec.absname].values())
    
scatrhos = np.vstack(scatrhos)

for state in ['desynch', 'synch']:
    rhos[state] = np.asarray(rhos[state])
    spars[state] = np.asarray(spars[state])

# scatter plot motion-PSTH correlation in desynched vs synched state:
figure(figsize=(3, 3))
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(scatrhos[:, 1], scatrhos[:, 0], 'o', mec='k', mfc='None')
text(0.02, 0.98, 'delay = %d ms' % CORRDELAYMS,
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
xlabel('synchronized', color='r')
ylabel('desynchronized', color='b')
xlim(-1.035, 1)
ylim(-1.035, 1)
gcfm().window.setWindowTitle('movie_global_motion_correlation_scatter_%dms' % CORRDELAYMS)
tight_layout(pad=0.3)
# report numbers, fractions and chi2 p values for PSTH-motion scatter plot:
nbelowmotionyxline = (scatrhos[:, 1] > scatrhos[:, 0]).sum()
nabovemotionyxline = (scatrhos[:, 0] > scatrhos[:, 1]).sum()
fractionbelowrelsyxline = nbelowmotionyxline / (nbelowmotionyxline + nabovemotionyxline)
chi2, p = chisquare([nabovemotionyxline, nbelowmotionyxline])
print('nbelowmotionyxline=%d, nabovemotionyxline=%d, fractionbelowrelsyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowmotionyxline, nabovemotionyxline, fractionbelowrelsyxline,
                             chi2, p))

# plot rho histograms:
figure(figsize=(3, 3))
dmean = rhos['desynch'].mean()
smean = rhos['synch'].mean()
u, p = mannwhitneyu(rhos['desynch'], rhos['synch']) # 1-sided
pstring = '$p<%g$' % ceilsigfig(p)
rhobins = np.arange(RHOMIN, RHOMAX+0.1, 0.1) # left edges + rightmost edge
nd = hist(rhos['desynch'], bins=rhobins, histtype='step', color='b')[0]
ns = hist(rhos['synch'], bins=rhobins, histtype='step', color='r')[0]
nmax = max(np.hstack([nd, ns]))
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
# draw arrows at means:
ah = nmax / 8 # arrow height
arrow(dmean, nmax, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='b')
arrow(smean, nmax, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='r')
xlim(xmin=RHOMIN, xmax=RHOMAX)
ylim(ymax=nmax*1.01)
# remove unnecessary decimal places:
rhoticks = ([-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1],
            ['-0.5', '-0.25', '0', '0.25', '0.5', '0.75', '1'])
xticks(*rhoticks)
yticks([0, nmax]) # turn off y ticks to save space
xlabel('PSTH-motion correlation')
ylabel('cell count')
text(0.98, 0.98, 'delay = %d ms' % CORRDELAYMS,
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
gcfm().window.setWindowTitle('movie_global_motion_correlation_%dms' % CORRDELAYMS)
tight_layout(pad=0.3)

# plot rho vs motion stimulus-response delay
#            delay  desynch  synch
rhovsdelay = [[-90, -0.019, -0.041],
              [-75, -0.009, -0.025],
              [-60,  0.000, -0.011],
              [-45,  0.010,  0.006],
              [-30,  0.020,  0.028],
              [-15,  0.031,  0.053],
              [  0,  0.039,  0.075],
              [ 15,  0.043,  0.090],
              [ 30,  0.044,  0.094],
              [ 45,  0.045,  0.085],
              [ 60,  0.048,  0.065],
              [ 75,  0.047,  0.036],
              [ 90,  0.041,  0.005],
              [105,  0.030, -0.022],
              [120,  0.019, -0.042]]
rhovsdelay = np.asarray(rhovsdelay)
figure(figsize=(3, 3))
plot(rhovsdelay[:, 0], rhovsdelay[:, 1], 'b.-')
plot(rhovsdelay[:, 0], rhovsdelay[:, 2], 'r.-')
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
ymax = 0.12
ah = 0.021
arrow(30, ymax, 0, -ah, head_width=10, head_length=ah/2, length_includes_head=True,
      color='k')
ylim(-0.05, ymax)
yticks(np.arange(-0.05, ymax, 0.05))
xlabel('delay (ms)')
ylabel('PSTH-motion correlation')
gcfm().window.setWindowTitle('movie_global_motion_rhovsdelay')
tight_layout(pad=0.3)

# scatter plot PSTH sparseness vs. movie motion sparseness. No significant correlation in
# either state:
figure(figsize=(3, 3))
plot(spars['synch'][:, 0], spars['synch'][:, 1], 'r.', ms=2)
plot(spars['desynch'][:, 0], spars['desynch'][:, 1], 'b.', ms=2)
xlim(0, 1)
ylim(0, 1)
xr = np.asarray(xlim()) # x range
m, b, r, p, stderr = linregress(spars['synch'])
plot(xr, m*xr+b, 'r', ls='--', lw=2, alpha=0.7)
text(0.02, 0.02, 'r=%.2f, p=%.1g' % (r, p), color='r',
     transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
m, b, r, p, stderr = linregress(spars['desynch'])
plot(xr, m*xr+b, 'b', ls='--', lw=2, alpha=0.7)
text(0.02, 0.1, 'r=%.2f, p=%.1g' % (r, p), color='b',
     transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
allspars = np.vstack(spars.values())
m, b, r, p, stderr = linregress(allspars)
plot(xr, m*xr+b, 'k', ls='--', lw=2, alpha=0.7)
text(0.02, 0.18, 'r=%.2f, p=%.1g' % (r, p), color='k',
     transform=gca().transAxes, horizontalalignment='left', verticalalignment='bottom')
xlabel('movie motion sparseness')
ylabel('PSTH sparseness')
gcfm().window.setWindowTitle('movie_PSTH_sparseness')
tight_layout(pad=0.3)

pl.show()
