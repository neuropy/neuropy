"""Plot natural scene movie trial rasters, and optionally runspeed and specgrams

Run from within neuropy using `run -i scripts/natstate_rasters.py`"""

import numpy as np
import matplotlib.pyplot as plt

## TODO: this script could iterate over multiple recordings, although that might
## generate excessive figures
## TODO: this script could probably be integrated into Recording.traster()
# specify recordings to analyze:
#catnsrecs = [ptc17.tr2b.r58, ptc18.tr1.r38, ptc18.tr2c.r58, ptc22.tr1.r08,
#             ptc22.tr1.r10, ptc22.tr4b.r49]
#mousensrecs = [nts174.tr2.r05, pvc107.tr1.r09, pvc113.tr1.r11]
#allrecs = catnsrecs + mousensrecs
rec = ptc22.tr1.r08
#rec = ptc22.tr1.r10
#rec = nts174.tr2.r05
#rec = pvc107.tr1.r09
#rec = pvc113.tr1.r11

#nids = None
nids = [17, 20, 32, 40, 74, 94]
#nids = [13, 23, 28, 34, 46, 52, 63, 66, 69, 79]

animal = rec.tr.animal

# raster plot options:
showstates = True # 'auto' or True, 'manual', False
sortstates = False # sort trials by cortical state
if sortstates: assert showstates
# for mouse data with opto trials, pool trials over opto values, thereby mixing opto
# and non-opto trials in the same raster plot, i.e., ignore opto state:
POOLOVEROPTO = True
marker = '|'
s = 4 # marker size
alpha = 1
c = (0, 0, 0, alpha) # give the ticks some transparency
slw = 5 # state line width
slalpha = 1 # state line alpha
slpos = -0.07 # x position of state line
rasfigwidth = 2.5
rasfigheightpadding = 0.675 # inches
rasfigheightperntrials = 0.01545 # inches
fc = 'w' # axes face colour

showITI = False
if animal.type == 'Mouse':
    if showITI:
        xlim = 0, 6 # s, 1 s ITI
        xticks = [0, 1, 2, 3, 4, 5]
    else:
        xlim = 0, 5 # s, without ITI
        xticks = [0, 1, 2, 3, 4]
elif animal.type == 'Cat':
    if showITI:
        xlim = 0, 5.5 # s, 1 s ITI
        xticks = [0, 1, 2, 3, 4, 5]
    else:
        xlim = 0, 4.5 # s, without ITI
        xticks = [0, 1, 2, 3, 4]
# True: use xlim to designate duration to plot for each trial;
# False: use stimOFFTime to designate duration to plot for each trial:
forcedt = True
ylabel = False
title = False

# start recording loop here:

# get stimulus info for this recording and its single experiment:
assert len(rec.e) == 1
e = rec.e0
t0s = e.ttranges[:, 0] # us
t1s = e.ttranges[:, 1] # us
if forcedt:
    t1s = t0s + xlim[1] * 1e6 # us
try:
    p = e.p # mouse stim params dict
except AttributeError:
    p = None
if p: # mouse recording with potentially interleaved movie trials and opto
    # trialiss is 2D array, each row is trial indices for its matching movie:
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
            trialiss = np.asarray(trialiss) # convert to 2D array
else: # cat recording with identical movie trials
    ntrials = len(e.ttranges)
    trialiss = np.arange(ntrials).reshape((1, -1)) # 0-based 2D array with 1 row
    # format movie name as for mouse:
    moviename = os.path.split(e.s['fname'])[-1]
    framerange = e.d['framei']
    moviename += '_%d-%d' % (framerange.start, framerange.stop)
    umovienames = [moviename]

if nids == None:
    nids = sorted(rec.n) # all neuron IDs, sorted
neurons = [ rec.n[nid] for nid in nids ] # list of neurons

if showstates:
    # show coloured state lines along left edge of raster plot
    si, tsi = rec.lfp.si(showstates=showstates, title=rec.absname+'.si()')
    stranges, states = rec.lfp.si_split(si, tsi) # state time ranges (sec) and values
    if sortstates:
        print('Sorting trials by state, excluding those that traverse states')
        # only keep trials that completely fall within a single state trange,
        # and sort trials by state instead of by time
        # first sort stranges by increasing state, use stable argsort to keep stranges
        # within a state sorted in time:
        statesortis = states.argsort(kind='mergesort')
        states = states[statesortis]
        stranges = stranges[statesortis]
        # filter and sort trialiss:
        strialiss = []
        for strange in stranges:
            strialis = np.where((strange[0] <= t0s/1e6) & (t1s/1e6 <= strange[1]))[0] # 0-based
            strialiss.append(strialis)
        strialis = np.concatenate(strialiss)
        oldtrialiss = trialiss # rename
        nconditions = len(oldtrialiss)
        # build up new trialiss as a list of lists,
        # can't use [[]]*nconditions because that makes copies!:
        trialiss = [ [] for cond in range(nconditions) ]
        for striali in strialis:
            # for each kept trial, find its row index in oldtrialiss, i.e. its movie/opto index,
            # there's probably a more vectorized way to do this:
            i = np.where(oldtrialiss == striali)[0]
            assert len(i) == 1
            i = int(i) # convert from array to int for indexing into list
            trialiss[i].append(striali)
    # find current state at start of each trial, and later use that to plot state colour
    # as a function of trial index:
    trial2state = {} # map trial index to state value
    for triali in np.concatenate(trialiss): # iterate over all trials, regardless of trial type
        t0 = t0s[triali] / 1e6 # sec
        statei, = np.where((stranges[:, 0] <= t0) & (t0 < stranges[:, 1]))
        assert len(statei) <= 1 # trial start should only match a single state trange
        if len(statei) == 1:
            state = states[statei[0]] # pull it out of the array, dereference
        else:
            state = None
        trial2state[triali] = state

# plot rasters: iterate over movies, units, trials
nrasterplots = 0
for trialtypei, (trialis, moviename) in enumerate(zip(trialiss, movienames)):
    ntrials = len(trialis) # number of trials for this movie
    rasfigheight = rasfigheightpadding + ntrials*rasfigheightperntrials
    rasfigsize = rasfigwidth, rasfigheight
    for neuron in neurons:
        spikes = neuron.spikes # spike times, us
        ts, subtrialis = [], [] # x and y values for all spikes in this trial raster plot
        for subtriali, triali in enumerate(trialis):
            t0, t1 = t0s[triali], t1s[triali]
            s0i, s1i = spikes.searchsorted([t0, t1])
            # this unit's spike times relative to start of this trial, in seconds:
            t = spikes[s0i:s1i] - t0
            nspikes = len(t) # if nspikes == 0, append empty arrays to ts and subtrialis
            ts.append(t) # x values for this trial, us
            # generate 0-based y values for spikes in this trial:
            subtrialis.append(np.tile(subtriali, nspikes))

        # convert spike times and trial indices to flat arrays:
        ts, subtrialis = np.hstack(ts), np.hstack(subtrialis)

        f = plt.figure(figsize=rasfigsize)
        a = f.add_subplot(111, fc=fc)

        if showstates:
            rasterstates = np.asarray([ trial2state[triali] for triali in trialis ])
            # don't try and plot state for trials with no state:
            validstateis = rasterstates != None
            statesubtrials = np.arange(ntrials)[validstateis] + 1 # 1-based
            # plot vertical lines just left of y axis:
            statelinepos = np.tile(slpos, len(statesubtrials))
        if showstates in [True, 'auto']:
            clrs = [ LFPPRBINCOLOURS[state] for state in rasterstates[validstateis] ]
            # the -0.5 and +0.5 center the lines around the raster ticks of single trials
            a.vlines(statelinepos, statesubtrials-0.5, statesubtrials+0.5,
                     colors=clrs, lw=slw, alpha=slalpha, clip_on=False)
        elif showstates == 'manual':
            raise NotImplementedError()
            '''
            dtrange, strange = np.asarray(REC2STATETRANGES[self.r.absname]) / 1e6
            dtrange = max(dtrange[0], t0), min(dtrange[1], t1) # clip desynch trange to t0, t1
            strange = max(strange[0], t0), min(strange[1], t1) # clip synch trange to t0, t1
            a.vlines(statelinepos, dtrange[0], dtrange[1], colors='b', lw=lw, alpha=alpha)
            a.vlines(statelinepos, strange[0], strange[1], colors='r', lw=lw, alpha=alpha)
            '''
        # plot 1-based trialis:
        a.scatter(ts/1e6, subtrialis+1, marker=marker, c=c, s=s, cmap=None) # in sec
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
        titlestr = rec.absname + '_' + moviename + '_n%d' % neuron.id
        if trialtypei == 1: # assume second row in trialiss (if any) is the opto trials
            titlestr += '_opto'
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents

        show() # call within the neuron loop, to ensure that rasters are displayed in nid order
        nrasterplots += 1

print('created %d raster plots' % nrasterplots)

## TODO: autoload runspeed info in Experiment.load_stim_mat
'''
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
'''

'''
# plot LFP specgram:
chani = -2 # 0-based chan ID
r.lfp.specgram(f1=59, chanis=chani, width=2, tres=0.5, cm='jet', relative2t0=True, title=False,
               reclabel=False, figsize=None)
titlestr = ename + '_specgram_chan%d' % (chani+1) # 1-based chan ID
gcfm().window.setWindowTitle(titlestr)
show()
'''
