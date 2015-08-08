"""Extract optic flow vector fields from natural scene movies, average them to calculate
global motion within each movie for each specified recording. Adapted from
opencv/samples/python2/opt_flow.py. Run from within neuropy using `run -i
scripts/movie_global_motion.py`"""

from __future__ import division, print_function

import cv2
import numpy as np

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

FIGSIZE = (8, 4)

# calculate optic flow vector field between neighbouring pairs of frames, average their
# magnitudes to get global motion:
motion, tmotion = {}, {}
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
    figure(figsize=FIGSIZE)
    plot(frameis[1:], motion[name], 'k-', lw=1.5)
    xlabel('frame index')
    #plot(tmotion[name], motion[name], 'k-', lw=1.5)
    #xlabel('t (s)')
    ylabel('motion amplitude (deg/s)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_motion_%s' % name)
    tight_layout(pad=0.3)


"""Calculate PSTHs as in psth_precision.py, then correlate each one with its respective movie
motion signal."""

## TODO: also scatter plot their correlation on desynched vs synched axes.

import core
from core import intround
from psth_funcs import get_psth_peaks_gac

NIDSKIND = 'all' # 'active' or 'all'
BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
BLANK = False # consider blank periods between trials?
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
CORRDELAY = 0.105 # correlation delay to use between stimulus and response, sec
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

# correlate each PSTH in each recording section with the movie motion signal:
rhos = {'desynch': [], 'synch': []}
for rec, nids, strange, state in zip(recs, recsecnids, stranges, states):
    print(rec.absname, state)
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK, strange=strange,
                                 plot=False, binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
    m = motion[rec.absname]
    tm = tmotion[rec.absname]
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
        rhos[state].append(core.corrcoef(m[motionis], psth[psthis]))
    print('\n') # two newlines

for state in ['desynch', 'synch']:
    rhos[state] = np.asarray(rhos[state])

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
xlabel('motion correlation')
ylabel('cell count')
text(0.98, 0.98, '$\mu$ = %.3f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.90, '$\mu$ = %.3f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.82, '%s' % pstring, color='k',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
gcfm().window.setWindowTitle('movie_global_motion_correlation_%dms' % intround(CORRDELAY*1000))
tight_layout(pad=0.3)

pl.show()
