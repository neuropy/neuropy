"""Load blab natmov spike times, LFP and runspeed from .mat files, plot rasters and specgrams"""

import numpy as np
import pylab as pl
from scipy.io import loadmat

from lfp import LFP
import core

# raster plot options:
rasfname = '/home/mspacek/dev/blab/natmov/results/PVCre_0113/s01/PVCre_0113_s01_e11_rasters.mat'
marker = '|'
s = 4 # marker size
alpha = 1
c = (0, 0, 0, alpha) # give the ticks some transparency
rasfigsize = 3, 5 # inches
axisbg = 'w'
xlim = 0, 5.5 # s
ylabel = False
title = False

# plot all rasters:
rasmat = loadmat(rasfname)
ename = os.path.splitext(os.path.basename(rasfname))[0]
spiketss = rasmat['spiketss'][0] # array of spike times (s), indexed by 0-based unit ID
for uidi, spikets in enumerate(spiketss): # iterate over units
    if spikets.size == 0:
        continue
    raster = spikets[0] # pull out of enclosing array to get array of arrays, one row per trial
    ntrials = len(raster)
    ts = np.vstack(raster)[:, 0] # row vector
    trialis = np.hstack([ np.tile(triali, len(raster[triali])) for triali in range(ntrials) ])

    f = pl.figure(figsize=rasfigsize)
    a = f.add_subplot(111, axisbg=axisbg)

    # plot 1-based trialis:
    a.scatter(ts, trialis+1, marker=marker, c=c, s=s, cmap=None)
    a.set_xlim(xlim)
    # -1 inverts the y axis, +1 ensures last trial is fully visible:
    a.set_ylim(ntrials+1, -1)
    # turn off annoying "+2.41e3" type offset on x axis:
    formatter = mpl.ticker.ScalarFormatter(useOffset=False)
    a.xaxis.set_major_formatter(formatter)
    a.set_xlabel("trial time (s)")
    if ylabel:
        a.set_ylabel("trial index") # sweep index order, not necessarily temporal order
    else:
        a.set_yticks([]) # turn off y ticks
    titlestr = ename + '_u%d' % (uidi+1) # 1-based unit ID
    gcfm().window.setWindowTitle(titlestr)
    if title:
        a.set_title(titlestr)
    f.tight_layout(pad=0.3) # crop figure to contents


# LFP plot options:
lfpfname = '/home/mspacek/dev/blab/natmov/results/PVCre_0113/s01/PVCre_0113_s01_e11_LFP.mat'
lfpchanis = [31] #range(32) # 0-based
f1 = 59 # Hz
relative2t0 = True
lfpfigsize = 5, 2 # inches
title = False
width, tres = 2, 0.5 # s

# plot LFP spectrograms:
lfp = LFP(None, None)
lfpmat = loadmat(lfpfname)
ename = os.path.splitext(os.path.basename(lfpfname))[0].rstrip('_LFP') + '_specgram'
lfp.data = lfpmat['lfp'] * 1000 # convert from mV to uV
tlfp = intround(lfpmat['tlfp'][0] * 1e6) # convert from s to nearest us
lfp.t0, lfp.t1 = tlfp[0], tlfp[-1]
lfp.tres = intround((np.diff(tlfp)).mean()) # us
lfp.sampfreq = intround(1e6 / lfp.tres) # Hz
for chani in lfpchanis:
    lfp.specgram(f1=f1, chanis=chani, width=width, tres=tres, relative2t0=relative2t0,
                 title=title, reclabel=False, figsize=lfpfigsize)
    titlestr = ename + '_c%d' % (chani+1) # 1-based chan ID
    gcfm().window.setWindowTitle(titlestr)


# plot runspeed as a color map:
width, tres = 10, 0.5 # s
minspeed = 1 # cm/s, otherwise considered at rest
runfname = '/home/mspacek/dev/blab/natmov/results/PVCre_0113/s01/PVCre_0113_s01_e11_runspeed.mat'
runmat = loadmat(runfname)
ename = os.path.splitext(os.path.basename(runfname))[0]
tspeed = runmat['tspeed'][0] # s
speed = runmat['speed'][0] # cm/s
t0 = lfp.t0 / 1e6 # sec, tspeed seems to start from 0, lfp starts a few ms later
t1 = lfp.t1 / 1e6
tranges = core.split_tranges([(t0, t1)], width, tres) # overlapping time bins
tiranges = tspeed.searchsorted(tranges)
nbins = len(tiranges)
tbinspeed = tranges[:, 0]
binspeed = np.zeros(nbins)
for bini, (t0i, t1i) in enumerate(tiranges):
    binspeed[bini] = speed[t0i:t1i].mean()

#figure()
#plot(tbinspeed, binspeed, '-')
binspeed.shape = -1, 1 # column vector
# turn it into a binary rest/run signal:
binspeed[binspeed < minspeed] = 0 # at rest
binspeed[binspeed >= minspeed] = 1 # running
#binspeed /= binspeed.max() # normalize to range [0, 1]
figure(figsize=(1, 10))
axis('off')
pl.imshow(binspeed, aspect=0.025, cmap='gray_r') # white=rest, black=run
titlestr = ename
gcfm().window.setWindowTitle(titlestr)
