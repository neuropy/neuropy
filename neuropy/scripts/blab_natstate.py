"""Load blab natmov trial times and runspeed from .mat files, plot rasters and specgrams"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

from lfp import LFP
from core import recarray2dict


# specify blab mouse recording to analyze:
mname = 'Ntsr1-Cre_0174'
sid, eid = 2, 5
basepath = '/home/mspacek/data/blab/natstate'
mpath = os.path.join(basepath, mname)
ename = '%s_s%02d_e%02d' % (mname, sid, eid)
epath = os.path.join(mpath, 'tr%d' % sid, '%02d' % eid)
stimfname = ename + '_stim.mat'
rsfname = ename + '_runspeed.mat'
stimfullfname = os.path.join(epath, stimfname)
rsfullfname = os.path.join(epath, rsfname)

m = Animal(mpath)
m.load()
r = m.tr[sid].r[eid]
snids = sorted(r.n) # sorted neuron IDs
neurons = [ r.n[nid] for nid in snids ] # sorted list of neurons

# raster plot options:
marker = '|'
s = 4 # marker size
alpha = 1
c = (0, 0, 0, alpha) # give the ticks some transparency
rasfigsize = 3, 5 # inches
axisbg = 'w'
xlim = 0, 5.5 # s
ylabel = False
title = False

# load stimulus info:
stimd = loadmat(stimfullfname, squeeze_me=True) # dict
t0s, t1s = stimd['stimONTimes'], stimd['stimOFFTimes'] # sec
ntrials = len(t0s)
p = recarray2dict(stimd['p']) # stim parameters struct `p` converted to dict
# each row is trial indices for its matching movie:
trialiss = p['seqnums'] - 1 # convert seqnums from 1-based to 0-based
movienames = p['movie']

# plot rasters: iterate over movies, units, trials
for trialis, moviename in zip(trialiss, movienames):
    for neuron in neurons:
        spikes = neuron.spikes / 1e6 # spike times, in seconds
        ts, subtrialis = [], [] # x and y values for all spikes in this trial raster plot
        for triali in trialis:
            t0, t1 = t0s[triali], t1s[triali]
            s0i, s1i = spikes.searchsorted([t0, t1])
            # this unit's spike times relative to start of this trial, in seconds:
            t = spikes[s0i:s1i] - t0
            nspikes = len(t) # if nspikes == 0, append empty arrays to ts and subtrialis
            ts.append(t) # x values for this trial
            # generate 0-based y values for spikes in this trial:
            subtrialis.append(np.tile(triali, nspikes))

        # convert spike times and trial indices to flat arrays:
        ts, subtrialis = np.hstack(ts), np.hstack(subtrialis)

        f = plt.figure(figsize=rasfigsize)
        a = f.add_subplot(111, axisbg=axisbg)

        # plot 1-based trialis:
        a.scatter(ts, subtrialis+1, marker=marker, c=c, s=s, cmap=None)
        a.set_xlim(xlim)
        # -1 inverts the y axis, +1 ensures last trial is fully visible:
        a.set_ylim(ntrials+1, -1)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        a.set_xlabel("trial time (s)")
        if ylabel:
            a.set_ylabel("trial index") # trial index order, not necessarily temporal order
        else:
            a.set_yticks([]) # turn off y ticks
        titlestr = moviename + '_' + ename + '_n%d' % neuron.id
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents

        show() # call within the neuron loop, to ensure that rasters are displayed in nid order
