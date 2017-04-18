"""Load blab natmov trial times and runspeed from .mat files, plot rasters and specgrams"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

from lfp import LFP
from core import recarray2dict


# specify blab mouse recording to analyze:
mname, sid, eid = 'Ntsr1-Cre_0174', 2, 5 # awake, 3 interleaved movies
#mname, sid, eid = 'PVCre_0107', 1, 9 # anesth, 2 interleaved movies
#mname, sid, eid = 'PVCre_0113', 1, 11 # awake, 2 interleaved movies
POOLOVEROPTO = True
basepath = '/home/mspacek/data/blab/natstate'
mpath = os.path.join(basepath, mname)
ename = '%s_s%02d_e%02d' % (mname, sid, eid)
epath = os.path.join(mpath, 'tr%d' % sid, '%02d' % eid)
stimfname = ename + '_stim.mat'
rsfname = ename + '_runspeed.mat'
stimfullfname = os.path.join(epath, stimfname)
rsfullfname = os.path.join(epath, rsfname)

m = Animal(mpath)
exec('%s = m' % m.name) # add to local namespace via its short name
m.load()
r = m.tr[sid].r[eid]
snids = sorted(r.n) # sorted neuron IDs
neurons = [ r.n[nid] for nid in snids ] # sorted list of neurons

# raster plot options:
marker = '|'
s = 4 # marker size
alpha = 1
c = (0, 0, 0, alpha) # give the ticks some transparency
rasfigwidth = 3
rasfigheightoffset = 0.125 # inches
rasfigheightperntrials = (5 - 0.125) / 200 # inches per trial
fc = 'w' # axes face colour
xlim = 0, 6 # s, show full 1 s ITI
xticks = [0, 1, 2, 3, 4, 5]
# True: use xlim to designate duration to plot for each trial;
# False: use stimOFFTime to designate duration to plot for each trial:
forcedt = True
ylabel = False
title = False

# load stimulus info:
stimd = loadmat(stimfullfname, squeeze_me=True) # dict
t0s = stimd['stimONTimes'] # sec
if forcedt:
    t1s = t0s + xlim[1]
else:
    t1s = stimd['stimOFFTimes'] # sec
p = recarray2dict(stimd['p']) # stim parameters struct `p` converted to dict
# each row is trial indices for its matching movie:
trialiss = p['seqnums'] - 1 # convert seqnums from 1-based to 0-based
movienames = p['movie']
umovienames = np.unique(movienames)
# handle opto trials:
if len(movienames) != len(umovienames):
    # some movies were displayed multiple times, probably in combination with opto stim:
    optopari, = np.where(p['parnames'] == 'opto')
    optovals = p['pars'][optopari]
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

# plot rasters: iterate over movies, units, trials
nrasterplots = 0
for trialis, umoviename in zip(trialiss, umovienames):
    ntrials = len(trialis)
    rasfigheight = rasfigheightoffset + ntrials*rasfigheightperntrials
    rasfigsize = rasfigwidth, rasfigheight
    for neuron in neurons:
        spikes = neuron.spikes / 1e6 # spike times, in seconds
        ts, subtrialis = [], [] # x and y values for all spikes in this trial raster plot
        for subtriali, triali in enumerate(trialis):
            t0, t1 = t0s[triali], t1s[triali]
            s0i, s1i = spikes.searchsorted([t0, t1])
            # this unit's spike times relative to start of this trial, in seconds:
            t = spikes[s0i:s1i] - t0
            nspikes = len(t) # if nspikes == 0, append empty arrays to ts and subtrialis
            ts.append(t) # x values for this trial
            # generate 0-based y values for spikes in this trial:
            subtrialis.append(np.tile(subtriali, nspikes))

        # convert spike times and trial indices to flat arrays:
        ts, subtrialis = np.hstack(ts), np.hstack(subtrialis)

        f = plt.figure(figsize=rasfigsize)
        a = f.add_subplot(111, fc=fc)
        # plot 1-based trialis:
        a.scatter(ts, subtrialis+1, marker=marker, c=c, s=s, cmap=None)
        a.set_xlim(xlim)
        # -1 inverts the y axis, +1 ensures last trial is fully visible:
        a.set_ylim(ntrials+1, -1)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        a.set_xticks(xticks)
        a.set_xlabel("trial time (s)")
        if ylabel:
            a.set_ylabel("trial index") # trial index order, not necessarily temporal order
        else:
            a.set_yticks([]) # turn off y ticks
        titlestr = ename + '_' + umoviename + '_n%d' % neuron.id
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents

        show() # call within the neuron loop, to ensure that rasters are displayed in nid order
        nrasterplots += 1
print('created %d raster plots' % nrasterplots)

if os.path.exists(rsfullfname):
    # load runspeed info:
    rsd = loadmat(rsfullfname, squeeze_me=True) # dict
    tspeed, speed = rsd['tspeed'], rsd['speed'] # sec, cm/s, speed can contain nans

    # smooth runspeed using overlapping time bins:
    width, tres = 10, 0.5 # s
    minspeed = 1 # cm/s, otherwise considered at rest
    r.lfp.get_data() # make sure it's loaded so we can access t0 and t1
    t0 = r.lfp.t0 / 1e6 # sec, tspeed starts at t=0, lfp starts a few ms later
    t1 = r.lfp.t1 / 1e6
    tranges = core.split_tranges([(t0, t1)], width, tres) # overlapping time bins
    tiranges = tspeed.searchsorted(tranges)
    nbins = len(tiranges)
    tbinspeed = tranges[:, 0]
    binspeed = np.zeros(nbins)
    for bini, (t0i, t1i) in enumerate(tiranges):
        binspeed[bini] = np.nanmean(speed[t0i:t1i]) # handle nans with nanmean
    show()

    # plot runspeed as a color map:
    #figure()
    #plot(tbinspeed, binspeed, '-')
    binspeed.shape = -1, 1 # column vector
    if minspeed:
        # turn it into a binary rest/run signal:
        binspeed[binspeed < minspeed] = 0 # at rest
        binspeed[binspeed >= minspeed] = 1 # running
    f = figure(figsize=(1, 10))
    axis('off')
    plt.imshow(binspeed, aspect=0.025, cmap='gray_r') # white=rest, black=run
    #plt.imshow(binspeed, aspect=0.025, cmap='jet') # blue=rest, red=run
    titlestr = ename + '_runspeed'
    gcfm().window.setWindowTitle(titlestr)
    f.tight_layout(pad=0.3) # crop figure to contents
    show()

# plot LFP specgram:
chani = -2 # 0-based chan ID
r.lfp.specgram(f1=59, chanis=chani, width=2, tres=0.5, cm='jet', relative2t0=True, title=False,
               reclabel=False, figsize=None)
titlestr = ename + '_specgram_chan%d' % (chani+1) # 1-based chan ID
gcfm().window.setWindowTitle(titlestr)
show()
