"""Load blab natmov spike times, LFP and runspeed from .mat files, plot rasters and specgrams"""

import numpy as np
import pylab as pl
from scipy.io import loadmat

from lfp import LFP
import core
from core import toiter


## NOTE: this one has 2 interleaved movies: normal and contrast-inverted!
basepath = '/home/mspacek/blab/natstate/results/PVCre_0113/s01'

# raster plot options:
rasfname = 'PVCre_0113_s01_e11_rasters.mat'
rasfullfname = os.path.join(basepath, rasfname)
marker = '|'
s = 4 # marker size
alpha = 1
c = (0, 0, 0, alpha) # give the ticks some transparency
rasfigsize = 3, 5 # inches
axisbg = 'w'
xlim = 0, 5.5 # s
xticks = [0, 1, 2, 3, 4, 5]
ylabel = False
title = False

# plot all rasters:
rasd = loadmat(rasfullfname, squeeze_me=True) # dict
ename = os.path.splitext(rasfname)[0]
spiketss = rasd['spiketss'] # array of spike times (s), indexed by 0-based unit ID
for uid, raster in enumerate(spiketss): # iterate over units
    ntrials = len(raster)
    if ntrials == 0: # many of the entries are empty, presumably no spikes for most units?
        continue
    ts = np.hstack(raster) # flatten to 1D
    trialis = []
    for triali in range(ntrials):
        # due to squeeze_me=True, if there's only 1 spike time, it's a scalar, so use toiter()
        nspikes = len(toiter(raster[triali]))
        trialis.append(np.tile(triali, nspikes))
    trialis = np.hstack(trialis)

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
    a.set_xticks(xticks)
    a.set_xlabel("trial time (s)")
    if ylabel:
        a.set_ylabel("trial index") # sweep index order, not necessarily temporal order
    else:
        a.set_yticks([]) # turn off y ticks
    titlestr = ename + '_u%d' % (uid+1) # 1-based unit ID
    gcfm().window.setWindowTitle(titlestr)
    if title:
        a.set_title(titlestr)
    f.tight_layout(pad=0.3) # crop figure to contents
    show() # call within the neuron loop, to ensure that rasters are displayed in nid order

# LFP plot options:
lfpfname = 'PVCre_0113_s01_e11_LFP.mat'
lfpfullfname = os.path.join(basepath, lfpfname)
lfpchanis = [31] #range(32) # 0-based
f1 = 59 # Hz
relative2t0 = True
lfpfigsize = 5, 2 # inches
title = False
width, tres = 2, 0.5 # s

# plot LFP spectrograms:
lfp = LFP(None, None)
lfpd = loadmat(lfpfullfname, squeeze_me=True) # dict
ename = os.path.splitext(lfpfname)[0].rstrip('_LFP')
lfp.data = lfpd['lfp'] * 1000 # convert from mV to uV
tlfp = intround(lfpd['tlfp'] * 1e6) # convert from s to nearest us
lfp.t0, lfp.t1 = tlfp[0], tlfp[-1]
lfp.tres = intround((np.diff(tlfp)).mean()) # us
lfp.sampfreq = intround(1e6 / lfp.tres) # Hz
for chani in lfpchanis:
    lfp.specgram(f1=f1, chanis=chani, width=width, tres=tres, cm='jet', relative2t0=relative2t0,
                 title=title, reclabel=False, figsize=lfpfigsize)
    titlestr = ename + '_specgram_c%d' % (chani+1) # 1-based chan ID
    gcfm().window.setWindowTitle(titlestr)

# plot runspeed as a color map:
width, tres = 10, 0.5 # s
minspeed = 1 # cm/s, otherwise considered at rest
rsfname = 'PVCre_0113_s01_e11_runspeed.mat'
rsfullfname = os.path.join(basepath, rsfname)
rsd = loadmat(rsfullfname, squeeze_me=True)
ename = os.path.splitext(rsfname)[0]
tspeed = rsd['tspeed'] # s
speed = rsd['speed'] # cm/s
t0 = lfp.t0 / 1e6 # sec, tspeed seems to start from 0, lfp starts a few ms later
t1 = lfp.t1 / 1e6
tranges = core.split_tranges([(t0, t1)], width, tres) # overlapping time bins
tiranges = tspeed.searchsorted(tranges)
nbins = len(tiranges)
tbinspeed = tranges[:, 0]
binspeed = np.zeros(nbins)
for bini, (t0i, t1i) in enumerate(tiranges):
    binspeed[bini] = speed[t0i:t1i].mean()
show()

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
show()
