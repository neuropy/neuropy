"""Defines the Recording class"""

from __future__ import division

import os
import StringIO
import random

import numpy as np
import scipy.stats

import pylab as pl
from pylab import get_current_fig_manager as gcfm
import matplotlib as mpl

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

import core
from core import (LFP, SpatialPopulationRaster, DensePopulationRaster, Codes, SpikeCorr,
                  binarray2int, nCrsamples, iterable, entropy_no_sing, lastcmd, intround,
                  tolist, rstrip, dictattr, warn, pmf, TAB)
from colour import CLUSTERCOLOURDICT
from experiment import Experiment
from sort import Sort
'''
# Good global setting for presentation plots:
pl.rcParams['axes.labelsize'] = 30
pl.rcParams['xtick.labelsize'] = 25
pl.rcParams['ytick.labelsize'] = 25
pl.rcParams['xtick.major.size'] = 7
pl.rcParams['ytick.major.size'] = 7
pl.rcParams['lines.markersize'] = 10
# use gca().set_position([0.15, 0.15, 0.8, 0.8]) or just the 'configure subplots' widget to
# make all the labels fit within the figure
'''

class BaseRecording(object):
    """A recording corresponds to a single SURF file, ie everything recorded between
    when the user hits record and when the user hits stop and closes the SURF file,
    including any pauses in between experiments within that recording. A recording
    can have multiple experiments, and multiple spike extractions, called sorts"""
    def __init__(self, path, track=None):
        self.level = 3 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # string buffer to print tree hierarchy to
        self.path = path
        self.tr = track
        if track != None:
            # update parent track's recording dict, in case self wasn't loaded by its parent
            track.r[self.id] = self
        self.e = dictattr() # store experiments in a dictionary with attrib access
        self.sorts = dictattr() # store sorts in a dictionary with attrib access

    def get_name(self):
        return os.path.split(self.path)[-1]

    name = property(get_name)

    def get_absname(self):
        """Absolute name including parent animal and track, kind of like absolute path, but
        more abbreviated, as one would enter it at the IPython prompt"""
        return '.'.join((self.tr.absname, 'r'+self.id))

    absname = property(get_absname)

    def get_id(self):
        # return the first word in the name, using - or whitespace as separators
        return self.name.split('-')[0].split(' ')[0]

    id = property(get_id)

    # shortcuts to various attribs and properties in default sort:
    n = property(lambda self: self.sort.n)
    qn = property(lambda self: self.sort.qn)
    alln = property(lambda self: self.sort.alln)
    nspikes = property(lambda self: self.sort.nspikes)
    nneurons = property(lambda self: self.sort.nneurons)
    nqneurons = property(lambda self: self.sort.nqneurons)
    nallneurons = property(lambda self: self.sort.nallneurons)
    datetime = property(lambda self: self.sort.datetime)
    pttype = property(lambda self: self.sort.pttype)
    chanpos = property(lambda self: self.sort.chanpos)

    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        if self.tr != None:
            self.tr.writetree(string)

    def load(self):
        treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)

        # Sorts from .ptcs files and .sort folders, and Experiments from .din files:
        allfdnames = os.listdir(self.path) # all file and dir names in self.path
        sortfdnames = []
        dinfnames = []
        lfpfnames = []
        for fdname in allfdnames:
            fullname = os.path.join(self.path, fdname)
            if os.path.isfile(fullname):
                if fdname.endswith('.ptcs'):
                    sortfdnames.append(fdname)
                elif fdname.endswith('.din'):
                    dinfnames.append(fdname)
                elif fdname.endswith('.lfp.zip'):
                    lfpfnames.append(fdname)
            elif os.path.isdir(fullname) and fdname.endswith('.sort'):
                sortfdnames.append(fdname)
        # sort filenames alphabetically, which should also be chronologically:
        sortfdnames.sort()
        dinfnames.sort()
        lfpfnames.sort()

        # load all Sorts, or just the most recent one:
        uns = get_ipython().user_ns
        if not uns['LOADALLSORTS']:
            sortfdnames = [sortfdnames[-1]] # just the most recent one
        for sortid, fdname in enumerate(sortfdnames):
            path = os.path.join(self.path, fdname)
            sort = Sort(path, id=sortid, recording=self)
            sort.load()
            self.sorts[sort.name] = sort # save it
            self.__setattr__('sort' + str(sort.id), sort) # add shortcut attrib
        # make last sort the default one
        self.sort = self.sorts[sortfdnames[-1]]

        # load all .din as Experiments:
        for expid, fname in enumerate(dinfnames): # expids follow order in dinfnames
            path = os.path.join(self.path, fname)
            experiment = Experiment(path, id=expid, recording=self)
            experiment.load()
            self.e[experiment.id] = experiment
            self.__setattr__('e' + str(experiment.id), experiment) # add shortcut attrib

        # load any LFP data from a .lfp.zip file:
        nlfpfiles = len(lfpfnames)
        if nlfpfiles == 0:
            pass
        elif nlfpfiles == 1:
            fullname = os.path.join(self.path, lfpfnames[0])
            self.lfp = LFP(self, fullname)
            #self.lfp.load() # for speed, don't load LFP data automatically
        else:
            raise RuntimeError("%d .lfp.zip files in %s, don't know which one to load"
                               % (nlfpfiles, self.path))

        if len(self.e) > 0:
            eids = list(self.e)
            e0, e1 = self.e[min(eids)], self.e[max(eids)]
            # start of the first experiment to end of the last one
            self.trange = e0.trange[0], e1.trange[1]
        else:
            # self.e is empty, no experiments in this recording, use first and last
            # spike across all neurons
            tranges = np.asarray([ n.trange for n in self.alln.values() ])
            self.trange = min(tranges[:, 0]), max(tranges[:, 1])

        # these are static, no need for properties:
        self.dt = self.trange[1] - self.trange[0] # duration (us)
        self.dtsec = self.dt / 1e6
        self.dtmin = self.dtsec / 60
        self.dthour = self.dtmin / 60

        self.calc_meanrates()

    def get_ordnids(self):
        """Return nids of active neurons in vertical spatial order, superficial to deep"""
        # numerical order, isn't necessarily the same as spatial order:
        nids = np.sort(self.n.keys())
        ypos = [ self.n[nid].pos[1] for nid in nids ]
        sortis = np.argsort(ypos) # superficial to deep
        return nids[sortis]

    def get_nids(self, tranges=None):
        """Find nids of neurons that are active in all tranges"""
        if tranges == None:
            return sorted(self.n) # return sorted nids of all active neurons
        # start with all neurons, even those with average rates below MINRATE over the
        # span of self. Remove them one by one if their average rates fall below MINRATE
        # in any trange in tranges
        #print('rid: %s' % self.id)
        tranges = np.asarray(tranges)
        assert tranges.ndim == 2 # 2D
        assert tranges.shape[1] == 2 # two columns
        uns = get_ipython().user_ns
        alln = self.alln
        nids = alln.keys()
        for trange in tranges:
            #print('trange: %r' % (trange,))
            dt = (trange[1] - trange[0]) / 1e6 # trange duration in sec
            assert dt >= 0
            nidi = 0
            while nidi < len(nids):
                #print('nids: %r' % nids)
                nid = nids[nidi]
                lo, hi = alln[nid].spikes.searchsorted(trange)
                nspikes = hi - lo # nspikes for nid in this trange
                meanrate = nspikes / dt # Hz
                if meanrate < uns['MINRATE']:
                    nids.remove(nid)
                    # new nid has slid into view, don't inc nidi
                else: # keep nid (for now), inc nidi
                    nidi += 1
        return np.sort(nids) # may as well sort them

    def mua(self, width=None, tres=None, neurons=None, smooth=False, plot=True):
        """Calculate and optionally plot multiunit activity as a function of time. neurons can
        be None, 'quiet', 'all', or a dict. `width' and `tres' of time bins are in seconds. If
        `smooth' is True, convolve with a smoothing window of width `width'"""
        if neurons == None:
            neurons = self.n # use active neurons
        elif neurons == 'quiet':
            neurons = self.qn # use quiet neurons
        elif neurons == 'all':
            neurons = self.alln # use all neurons
        nn = len(neurons)

        uns = get_ipython().user_ns
        if width == None:
            width = uns['MUAWIDTH'] # sec
        if tres == None:
            tres = uns['MUATRES'] # sec
        assert tres <= width
        width = intround(width * 1000000) # convert from sec to us
        tres = intround(tres * 1000000) # convert from sec to us

        nids = np.sort(neurons.keys())
        ys = np.array([ neurons[nid].pos[1] for nid in nids ]) # y positions of each neuron
        supis, midis, deepis = core.laminarity(ys)
        supnids = nids[supis]
        midnids = nids[midis]
        deepnids = nids[deepis]
        nsup = len(supnids)
        nmid = len(midnids)
        ndeep = len(deepnids)

        if nn == 0: allspikes = np.array([])
        else: allspikes = np.concatenate([ neurons[nid].spikes for nid in nids ])
        if nsup == 0: supspikes = np.array([])
        else: supspikes = np.concatenate([ neurons[nid].spikes for nid in supnids ])
        if nmid == 0: midspikes = np.array([])
        else: midspikes = np.concatenate([ neurons[nid].spikes for nid in midnids ])
        if ndeep == 0: deepspikes = np.array([])
        else: deepspikes = np.concatenate([ neurons[nid].spikes for nid in deepnids ])
        allspikes.sort() # sorted spikes from all neurons
        supspikes.sort() # sorted spikes from superficial neurons
        midspikes.sort() # sorted spikes from middle neurons
        deepspikes.sort() # sorted spikes from deep neurons

        if smooth:
            t0, t1 = self.trange
            # non-overlapping bin edges, including rightmost bin:
            edges = np.arange(t0, t1+tres, tres) # in us
            # get midpoint of each bin, convert from us to sec:
            t = (edges[:-1] + tres/2) / 1000000
        else:
            # potentially overlapping bin time ranges:
            tranges = core.split_tranges([self.trange], width, tres) # in us
            # get midpoint of each trange, convert from us to sec:
            t = tranges.mean(axis=1) / 1000000

        # in spikes/sec (Hz) per neuron:
        if smooth:
            allrates = self.calc_mua_smooth(allspikes, nn, edges, width)
            suprates = self.calc_mua_smooth(supspikes, nsup, edges, width)
            midrates = self.calc_mua_smooth(midspikes, nmid, edges, width)
            deeprates = self.calc_mua_smooth(deepspikes, ndeep, edges, width)
        else:
            allrates = self.calc_mua(allspikes, nn, tranges)
            suprates = self.calc_mua(supspikes, nsup, tranges)
            midrates = self.calc_mua(midspikes, nmid, tranges)
            deeprates = self.calc_mua(deepspikes, ndeep, tranges)
        rates = np.vstack([allrates, suprates, midrates, deeprates])
        n = nn, nsup, nmid, ndeep
        if plot:
            self.plot_mua(rates, t, n)
        return rates, t, n

    def calc_mua(self, spikes, nn, tranges):
        """Take sorted spike train from nn neurons and set of tranges and return mean
        spike rate as a function of time"""
        spikeis = spikes.searchsorted(tranges)
        counts = spikeis[:, 1] - spikeis[:, 0]
        widths = (tranges[:, 1] - tranges[:, 0]) / 1000000 # width of each trange, in sec
        if nn == 0:
            nn = 1 # prevent div by 0, 0 neurons result in 0 rates anyway
        return counts / widths / nn # in spikes/sec (Hz) per neuron

    def calc_mua_smooth(self, spikes, nn, edges, width):
        """Take sorted multiunit spike train, bin it, and return mean smoothed multiunit
        firing rate signal, in spikes/sec (Hz)"""
        spikehist = np.histogram(spikes, bins=edges)[0]
        tres = edges[1] - edges[0] # bin width (us)
        tressec = tres / 1000000
        nw = width / tres # window width, in number of bins
        window = np.hanning(nw)[nw/2:] # half a hanning window, causal (convolve flips it)
        if nn == 0:
            nn = 1 # prevent div by 0, 0 neurons result in 0 rates anyway
        # normalize by bin width and window area and nn to get spikes/sec (Hz) per neuron:
        return np.convolve(spikehist, window, mode='same') / tressec / window.sum() / nn
        # might be faster:
        #scipy.signal.fftconvolve(spikehist, window, mode='same') / tressec / window.sum() / nn

    def plot_mua(self, rates, t, n, figsize=(20, 6.5)):
        """Plot multiunit activity (all, sup, mid and deep firing rates) as a function of
        time"""
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.plot(t, rates[0], 'e', label='all (%d)' % n[0])
        a.plot(t, rates[1], 'r', label='superficial (%d)' % n[1])
        a.plot(t, rates[2], 'g', label='middle (%d)' % n[2])
        a.plot(t, rates[3], 'b', label='deep (%d)' % n[3])
        a.set_xlabel("time (sec)")
        a.set_ylabel("mean MUA (Hz)")
        # limit plot to duration of acquistion, in sec:
        t0, t1 = np.asarray(self.trange) / 1000000
        a.set_xlim(t0, t1)
        #a.autoscale(axis='x', enable=True, tight=True)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        uns = get_ipython().user_ns
        a.text(0.998, 0.99,
               '%s\n'
               'sup = %r um\n'
               'mid = %r um\n'
               'deep = %r um\n'
               % (self.name, uns['SUPRANGE'], uns['MIDRANGE'], uns['DEEPRANGE']),
               color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents

    def mua_si(self, smooth=False, chani=-1, ratio='L/(L+H)', figsize=(7.5, 6.5)):
        """Scatter plot multiunit activity vs LFP synchrony index"""
        mua, muat, n = self.mua(smooth=smooth, plot=False)
        si, sit = self.lfp.si(chani=chani, ratio=ratio, plot=False)
        # get common time resolution:
        if len(sit) > len(muat):
            siti = sit.searchsorted(muat)
            sitii = siti < len(sit) # prevent right side out of bounds indices into si
            muat = muat[sitii]
            mua = mua[:, sitii]
            siti = siti[sitii]
            sit = sit[siti]
            si = si[siti]
        else:
            muati = muat.searchsorted(sit)
            muatii = muati < len(muat) # prevent right side out of bounds indices into mua
            sit = sit[muatii]
            si = si[muatii]
            muati = muati[muatii]
            muat = muat[muati]
            mua = mua[:, muati]

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        ylim = mua.min(), mua.max()
        yrange = ylim[1] - ylim[0]
        extra = yrange * 0.03 # 3 %
        ylim = max(ylim[0]-extra, 0), ylim[1]+extra # don't go below 0
        sirange = np.array([0, 1])

        # plot separate regression lines for all, superficial, middle and deep cells:
        m0, b0, r0, p0, stderr0 = scipy.stats.linregress(si, mua[0])
        m1, b1, r1, p1, stderr1 = scipy.stats.linregress(si, mua[1])
        m2, b2, r2, p2, stderr2 = scipy.stats.linregress(si, mua[2])
        m3, b3, r3, p3, stderr3 = scipy.stats.linregress(si, mua[3])
        a.plot(sirange, m0*sirange+b0, 'e--')
        a.plot(sirange, m1*sirange+b1, 'r--')
        a.plot(sirange, m2*sirange+b2, 'g--')
        a.plot(sirange, m3*sirange+b3, 'b--')

        # scatter plot MUA vs SI:
        a.plot(si, mua[0], 'e.', label='all (%d), m=%.3f, r=%.3f' % (n[0], m0, r0))
        a.plot(si, mua[1], 'r.', label='sup (%d), m=%.3f, r=%.3f' % (n[1], m1, r1))
        a.plot(si, mua[2], 'g.', label='mid (%d), m=%.3f, r=%.3f' % (n[2], m2, r2))
        a.plot(si, mua[3], 'b.', label='deep (%d), m=%.3f, r=%.3f' % (n[3], m3, r3))

        a.set_xlim(0, 1)
        a.set_ylim(ylim)
        #a.autoscale(enable=True, axis='y', tight=True)
        a.set_xlabel("LFP synchrony index (%s)" % ratio)
        a.set_ylabel("mean MUA (Hz)")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        uns = get_ipython().user_ns
        a.text(0.998, 0.99,
               '%s\n'
               'sup = %r um\n'
               'mid = %r um\n'
               'deep = %r um\n'
               % (self.name, uns['SUPRANGE'], uns['MIDRANGE'], uns['DEEPRANGE']),
               color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents

    def calc_meanrates(self):
        """Calculate mean firing rates of all neurons in this recording"""
        RECNEURONPERIOD = get_ipython().user_ns['RECNEURONPERIOD']
        if RECNEURONPERIOD == 'recording':
            # calc n.meanrate using entire track duration:
            for n in self.alln.values():
                n.meanrate = n.nspikes / self.dtsec
        elif RECNEURONPERIOD == 'trange':
            # calc n.meanrate using duration between its first and last spike:
            for n in self.alln.values():
                if n.dtsec == 0:
                    n.meanrate = 0.0
                else:
                    n.meanrate = n.nspikes / n.dtsec
        else:
            raise ValueError("invalid value for RECNEURONPERIOD: %r" % RECNEURONPERIOD)

    def get_meanrates(self):
        """Return mean firing rates of all neurons in this recording"""
        return np.asarray([ n.meanrate for n in self.alln.values() ])

    meanrates = property(get_meanrates)

    def meanratepdf(self, bins=None, figsize=(7.5, 6.5)):
        """Plot histogram of mean firing rates"""
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        if bins == None:
            bins = np.arange(0, 1, 0.05)
        n, mr = np.histogram(self.meanrates, bins=bins, density=False)
        binwidth = mr[1] - mr[0] # take width of first bin
        a.bar(left=mr[:-1], height=n, width=binwidth, bottom=0, color='k', ec='k')
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        a.set_xlabel('mean firing rate (Hz)')
        a.set_ylabel('neuron count')
        f.tight_layout(pad=0.3) # crop figure to contents

    def collectcchs(self, nids, trange, bins, shiftcorrect=False, normalize=False):
        """Collect cross-correlation histograms for all pairs in nids. trange and bins are in
        us. If shiftcorrect, then shift correct the CCH by subtracting the CCH generated by
        shifting one spike train in each pair by some random amount. If normalize, weight each
        pair equally."""
        n = self.n
        nn = len(nids)
        hists = []
        if shiftcorrect:
            nshifts = int(shiftcorrect)
            # make sure no part of the shift corrector overlaps with trange:
            width = trange[1] - trange[0]
        for nii0 in range(nn):
            for nii1 in range(nii0+1, nn):
                spikes0 = n[nids[nii0]].spikes
                spikes1 = n[nids[nii1]].spikes
                dts = util.xcorr(spikes0, spikes1, trange) # spike time differences in us
                hist = np.histogram(dts, bins=bins)[0]
                if shiftcorrect:
                    # randomize amplitude and sign of shifts, shifting by anywhere from
                    # +/- width to 2*width
                    tshifts = intround(width + width*np.random.random(nshifts))
                    tshifts *= core.randsign(nshifts)
                    shifthist = []
                    for shifti in range(nshifts):
                        shiftdts = util.xcorr(spikes0, spikes1+tshifts[shifti], trange)
                        shifthist.append(np.histogram(shiftdts, bins=bins)[0])
                    shifthist = np.vstack(shifthist).mean(axis=0) # average shift predictor
                    hist -= shifthist # subtract average shift predictor
                # if we don't normalize, we treat our confidence in the CCH of a cell pair
                # proportionally to the number of spikes in that pair, which may be the
                # optimal thing to do. Otherwise, if we do normalize, we treat the CCH
                # of each pair equally, and therefore imply equal confidence in the
                # CCH of all pairs.
                if normalize:
                    hist = hist / hist.sum() # pmf: normalize so that sum of each hist is 1
                else:
                    hists.append(hist)
        hists = np.vstack(hists) # return array, with one row per CCH
        return hists

    def meancch(self, trange=(-100, 100), binw=2, shiftcorrect=False,
                shufflenids=False, subtract=False, figsize=(7.5, 6.5)):
        """Calculate mean cross-correlation histogram for all possible pairs of spike trains
        in self. trange and binw are in ms. If shiftcorrect, then shift correct the CCH by
        subtracting the CCH generated by shifting one spike train in each pair by some random
        amount. If subtract, take the mean CCH from shuffling nids and subtract it from the
        mean CCH from ordered nids, and plot a histogram of their difference"""
        ## TODO: add some kind of shuffle-correction, not over trials but by simply adding a
        ## random delay (greater than 2*max(trange)) to one of each pair of spike trains
        trange = np.asarray(trange) * 1000 # us
        assert trange[0] < trange[1]
        binw = binw * 1000 # us
        nids = self.get_ordnids() # in vertical spatial order
        """
        If shufflenids is nonzero, do that many runs, shuffling the nids on each run. nids by
        default are sorted by depth. shufflenids therefore allows examination of whether
        assymmetry in the resulting mean CCH is a result of some kind of causal regularity as
        a function of cell depth. To do this properly, when shuffling, should do this super
        averaging over a large number of runs, which would guarantee no assymmetry in the
        super averaged CCH:
        """
        nruns = 1
        if shufflenids:
            nruns = int(shufflenids)
        nn = len(nids)
        bins = np.arange(trange[0], trange[1]+binw, binw)
        hists = []
        for runi in range(nruns):
            if shufflenids:
                np.random.shuffle(nids) # in place
            hists.append(self.collectcchs(nids, trange, bins, shiftcorrect, normalize=False))
        #hists = np.vstack(hists).sum(axis=0)
        hists = np.vstack(hists).mean(axis=0)
        # disallow any negative bin counts from shift correction:
        if shiftcorrect:
            hists[hists < 0] = 0

        if subtract: # call collectcchs one more time, without shuffling
            nids = self.get_ordnids() # in vertical spatial order
            unshuffledhists = self.collectcchs(nids, trange, bins, shiftcorrect,
                                               normalize=False)
            unshuffledhists = unshuffledhists.mean(axis=0)
            hists = unshuffledhists - hists # unshuffled minus (potentially) shuffled

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.bar(bins[:-1]/1000, hists, width=binw/1000)
        a.set_xlim(trange/1000)
        if not subtract:
            a.set_ylim(ymin=0)
        a.set_xlabel('time (ms)')
        a.set_ylabel('mean bin count')
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        a.text(0.998, 0.99, '%s' % self.name, color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        f.tight_layout(pad=0.3)

    def pospdf(self, neurons=None, dim='y', nbins=10, a=None, stats=False, figsize=(7.5, 6.5)):
        """Plot PDF of cell positions ('x' or 'y') along the polytrode
        to get an idea of how cells are distributed in space"""
        if neurons == 'all':
            neurons = self.alln.values()
        elif neurons == 'quiet':
            neurons = self.qn.values()
        else:
            neurons = self.n.values()
        dimi = {'x':0, 'y':1}[dim]
        p = [ n.pos[dimi] for n in neurons ] # all position values
        nbins = max(nbins, 2*intround(np.sqrt(self.nneurons)))
        n, p = np.histogram(p, bins=nbins) # p includes rightmost bin edge
        binwidth = p[1] - p[0] # take width of first bin in p

        if stats:
            mean = np.mean(p)
            median = np.median(p)
            argmode = n.argmax()
            mode = p[argmode] + binwidth / 2 # middle of tallest bin
            stdev = np.std(p)

        if a == None:
            f = pl.figure(figsize=figsize)
            a = f.add_subplot(111)
        else: # add to existing axes
            a.hold(True)
            f = pl.gcf()

        # use CLUSTERCOLOURDICT for familiarity with len 10 1-based id to colour mapping,
        # remove any trailing alphas form id:
        color = CLUSTERCOLOURDICT[int(self.id.rstrip('bc'))]

        # exclude rightmost bin edge in p
        a.bar(left=p[:-1], height=n, width=binwidth, bottom=0, color=color, ec=color,
              yerr=None, xerr=None, capsize=3)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        a.set_xlabel('neuron %s position (um)' % dim)
        a.set_ylabel('neuron count')

        if stats:
            # add stuff to top right of plot:
            uns = get_ipython().user_ns
            a.text(0.99, 0.99, '%s\n'
                               'mean = %.3f\n'
                               'median = %.3f\n'
                               'mode = %.3f\n'
                               'stdev = %.3f\n'
                               'minrate = %.2f Hz\n'
                               'nneurons = %d\n'
                               'dt = %d min'
                               % (self.name, mean, median, mode, stdev,
                                  uns['MINRATE'], self.nneurons, intround(self.dtmin)),
                               transform = a.transAxes,
                               horizontalalignment='right',
                               verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        f.canvas.draw() # this is needed if a != None when passed as arg
        return a


class RecordingRevCorr(BaseRecording):
    """Mix-in class that defines reverse correlation related Recording methods"""
    def sta(self, **kwargs):
        assert len(self.e) == 1
        e = self.e[0]
        return e.sta(**kwargs)

    def stc(self, **kwargs):
        assert len(self.e) == 1
        e = self.e[0]
        return e.stc(**kwargs)


class RecordingRaster(BaseRecording):
    """Mix-in class that defines the raster related Recording methods"""
    def praster(self, t0=None, t1=None, neurons=None, norder=None, dense=True, units='sec'):
        """Create a population spike raster plot. neurons can be None, 'quiet', 'all', or a
        dict. norder can be a sequence of nids, designating what order to present them in
        the raster plot, from bottom to top. If set to True, the order is automatically
        determined by MDS of pairwise spike correlations. Default order is spatial"""
        if neurons == None:
            neurons = self.n # use active neurons
        elif neurons == 'quiet':
            neurons = self.qn # use quiet neurons
        elif neurons == 'all':
            neurons = self.alln # use all neurons
        if t0 == None:
            t0, t1 = self.trange # use full recording trange
        else:
            tx = {'us': 1, 'ms': 1000, 'sec': 1000000}[units]
            t0 *= tx # convert to us
            if t1 == None:
                t1 = t0 + 10000000
            else:
                t1 *= tx # convert to us
        trange = np.array([t0, t1])
        if norder == True:
            norder = self.cc().norder()
        if dense:
            return DensePopulationRaster(trange=trange, neurons=neurons, norder=norder,
                                         units=units, text=self.name)
        else:
            return SpatialPopulationRaster(trange=trange, neurons=neurons, norder=norder,
                                           units=units, text=self.name)

    def traster(self, nids=None, eid=0, t0=None, dt=None, blank=True, s=20,
                figsize=(7.5, None)):
        """Create a trial spike raster plot for given neurons, based on stimulus info
        in experiment eid. blank designates whether to include blank frames for trials in
        movie type stimuli"""
        if nids == None:
            nids = sorted(self.n.keys()) # use active neurons
        elif nids == 'quiet':
            nids = sorted(self.qn.keys()) # use quiet neurons
        elif nids == 'all':
            nids = sorted(self.alln.keys()) # use all neurons
        else:
            nids = tolist(nids) # use specified neurons
        e = self.e[0]
        if 'framei' in e.vs: # movie type of stimulus, where each frame is a sweep
            trialtype = 'dinrange' # one trial for every cycle of din values
        else:
            trialtype = 'dinval' # one trial per block of identical din values
        din = e.din
        times = din[:, 0] # sweep times
        sweepis = din[:, 1] # sweep indices
        uns = get_ipython().user_ns
        NULLDIN = uns['NULLDIN']
        # find unique sweep indices, excluding NULLDIN
        usweepis = np.unique(sweepis)
        usweepis = usweepis[usweepis != NULLDIN]
        # find tranges of all trials, either manually based on t0 & dt, or automatically
        # based on trialtype:
        if dt != None:
            # assume all trials of equal length dt, starting from t0
            dt *= 1000000 # convert from sec to us
            if t0 == None:
                t0i = np.where(sweepis != NULLDIN)[0][0] # first non-NULL sweepi in din
                t0 = times[t0i] # in us
            else:
                t0 *= 1000000 # convert from sec to us
            tlast = times[-1]
            t0s = np.arange(t0, tlast-dt, dt)
            t1s = np.arange(t0+dt, tlast, dt)
            assert (t1s - t0s == dt).all()
            tranges = np.column_stack((t0s, t1s))
        elif trialtype == 'dinrange':
            sw0, sw1 = usweepis[0], usweepis[-1] # first and last sweep index in each trial
            i0s, = np.where(sweepis == sw0) # screen refresh indices for sw0
            # indices into i0s of start of each trange, prepend i0s with a value (-2)
            # guaranteed to be non-consecutive with the first value in i0s:
            i0is, = np.where(np.diff(np.hstack(([-2], i0s))) != 1)
            i0s = i0s[i0is]
            t0s = times[i0s]
            if not blank:
                i1s, = np.where(sweepis == sw1) # screen refresh indices for sw1
                # indices into i1s of end of each trange, append i1s with a value (-2)
                # guaranteed to be non-consecutive with the last value in i1s:
                i1is, = np.where(np.diff(np.hstack((i1s, [-2]))) != 1)
                i1s = i1s[i1is]
            else: # include blank frames
                # alternate method: only use sw0 to designate start and end of each trial,
                # and therefore include any blank periods at the end of each trial as a
                # part of that trial:
                i1s = i0s[1:] # missing one more at end at this point
                di1s = np.diff(i1s)
                maxdi1 = max(di1s)
                # append one more index interval to end of i1s
                i1s = np.hstack((i1s, [i1s[-1]+maxdi1]))
            t1s = times[i1s]
            tranges = np.column_stack((t0s, t1s))
        elif trialtype == 'dinval':
            t0s, t1s = [], []
            for sweepi in usweepis: # ordered by sweepi
                i, = np.where(sweepis == sweepi) # screen refresh indices
                # indices into i of start of each trange, prepend i with a value (-2)
                # guaranteed to be non-consecutive with the first value in i:
                i0is, = np.where(np.diff(np.hstack(([-2], i))) != 1)
                # indices into i of end of each trange, append i with a value (-2)
                # guaranteed to be non-consecutive with the last value in i:
                i1is, = np.where(np.diff(np.hstack((i, [-2]))) != 1)
                i0s = i[i0is]
                i1s = i[i1is]
                t0s.append(times[i0s])
                t1s.append(times[i1s])
                # each sweepi's tranges could also be saved into a dict, indexed by sweepi,
                # as in the tuning curve code
            t0s = np.hstack(t0s)
            t1s = np.hstack(t1s)
            tranges = np.column_stack((t0s, t1s))
        dts = t1s - t0s
        maxdt = max(dts) # max trial duration
        # depth of nids from top of electrode
        ypos = np.array([ self.alln[nid].pos[1] for nid in nids ])
        supis, midis, deepis = core.laminarity(ypos)
        nn = len(nids)
        ntrials = len(tranges)
        figsize = figsize[0], 1 + ntrials / 36 # ~1/36th vertical inch per trial
        for nidi, nid in enumerate(nids):
            spikes = self.alln[nid].spikes
            trials = []
            trialis = []
            for triali, trange in enumerate(tranges):
                si0, si1 = spikes.searchsorted(trange)
                # slice out spikes that fall within trials, make them relative to start of
                # each trial, convert from us to sec
                ts = (spikes[si0:si1] - trange[0]) / 1e6
                nspikes = len(ts)
                if nspikes == 0:
                    continue
                trials.append(ts) # x values for this trial
                trialis.append(np.tile(triali+1, nspikes)) # 1-based y values for this trial
            trials = np.hstack(trials)
            trialis = np.hstack(trialis)
            if supis[nidi]: c = 'r'
            elif midis[nidi]: c = 'g'
            elif deepis[nidi]: c = 'b'
            else: c = 'y'

            f = pl.figure(figsize=figsize)
            a = f.add_subplot(111)
            a.scatter(trials, trialis, marker='|', c=c, s=s)
            a.set_xlim(0, maxdt / 1e6) # sec
            # -1 inverts the y axis, +1 ensures last trial is fully visible:
            a.set_ylim(ntrials+1, -1)
            # turn off annoying "+2.41e3" type offset on x axis:
            formatter = mpl.ticker.ScalarFormatter(useOffset=False)
            a.xaxis.set_major_formatter(formatter)
            a.set_xlabel("time (sec)")
            a.set_ylabel("trial")
            titlestr = lastcmd() + " nid%d nidi%d" % (nid, nidi)
            gcfm().window.setWindowTitle(titlestr)
            a.set_title(titlestr)
            f.tight_layout(pad=0.3) # crop figure to contents
            self.f = f

    def tune(self, nids=None, eid=0, var='ori', fixed=None, tdelay=None):
        """Plot tuning curves for given neurons, based on stimulus info in experiment eid"""
        if nids == None:
            nids = sorted(self.n.keys()) # use active neurons
        elif nids == 'quiet':
            nids = sorted(self.qn.keys()) # use quiet neurons
        elif nids == 'all':
            nids = sorted(self.alln.keys()) # use all neurons
        else:
            nids = tolist(nids) # use specified neurons
        for nid in nids:
            n = self.alln[nid]
            t = n.tune(eid=eid, tdelay=tdelay)
            t.plot(var=var, fixed=fixed)
        

class RecordingCode(BaseRecording):
    """Mix-in class that defines spike code related methods"""
    def codes(self, nids=None, tranges=None, experiments=None, shufflecodes=False):
        """Return a Codes object, a 2D array where each row is a neuron code constrained
        to tranges, or to the tranges of experiments"""
        if nids == None:
            nids = self.get_nids() # sorted nids of all active neurons
        neurons = [ self.n[nid] for nid in nids ] # sorted list of neurons
        if tranges == None:
            tranges = [self.trange] # use whole Recording trange
            if experiments != None:
                tranges = [ e.trange for e in experiments ] # assume a list of Experiments
        codes = Codes(neurons=neurons, tranges=tranges, shufflecodes=shufflecodes)
        codes.calc()
        return codes
    '''
    # unused
    def spikecorr(self, nid1, nid2, tranges=None):
        """Calculate the correlation coefficient of the codes of two neurons"""
        code1 = self.n[nid1].code(tranges=tranges)
        code2 = self.n[nid2].code(tranges=tranges)
        return corrcoef(code1.c, code2.c)
    '''
    def sc(self, tranges=None, width=None, tres=None, weights=None, shift=0,
           shiftcorrect=0, experiments=None, nids=None, R=None):
        """Return a SpikeCorr object"""
        if weights in ['synch', 'desynch']:
            r, t = self.lfp.si(plot=False)
            if weights == 'synch':
                weights = r, t
            else: # weights == 'desynch'
                weights = 1-r, t
        sc = SpikeCorr(recording=self, tranges=tranges, width=width, tres=tres,
                       weights=weights, shift=shift, shiftcorrect=shiftcorrect,
                       experiments=experiments, nids=nids, R=None)
        # run sc.calc() as late as possible, not here
        return sc


class BaseNetstate(object):
    """Base class of Network state analyses.
    Implements a lot of the analyses on network states found in the 2006 Schneidman paper

    WARNING: not sure if self.tranges, which derives from self.experiments, is being
    used at all yet. See codes() method below.
    """
    def __init__(self, recording, experiments=None, nids=None):
        self.r = recording
        if experiments == None:
            self.tranges = [self.r.trange]
            # or should we check to see if this Recording has a tranges field due to
            # appending Neurons?
        else:
            experiments = tolist(experiments)
            try:
                # is experiments a list of Experiments?
                self.tranges = [ e.trange for e in experiments ]
            except AttributeError:
                # assume experiments is a list of experiment ids:
                self.tranges = [ self.r.e[ei].trange for ei in experiments ]
                # convert to a list of Experiments
                experiments = [ self.r.e[ei] for ei in experiments ]
        self.e = experiments # save list of Experiments (could potentially be None)
        self.neurons = self.r.n
        self.nneurons = len(self.neurons)
        if nids == None:
            self.nidswasNone = True
            nids = self.neurons.keys() # get all neuron indices in this Recording
            nids.sort() # make sure they're sorted
        else:
            self.nidswasNone = False
        self.cs = self.codes(nids=nids) # generate and save the Codes object for all the nids
        # if you need to retrieve the nids, you can get them from self.cs.nids. Leave
        # self.nids open for subclasses to use for their own purposes

    def codes(self, nids=None, shufflecodes=False):
        """Returns the appropriate Codes object, depending on the recording
        and experiments defined for this Netstate object"""
        ## TODO: codes are not currently constrained to when stimuli are on the screen,
        ## although this shouldn't be a big deal most of the time...
        return self.r.codes(nids=nids, experiments=self.e, shufflecodes=shufflecodes)

    def get_wordts(self, nids=None, mids=None):
        """Returns word times, ie the times of the left bin edges for which all the
        neurons in the mids in this Netstate object have a 1 in them, and all
        the rest have a 0 in them. nids lists the total population of neuron ids"""
        if nids == None:
            nids = self.cs.nids
        mids = toiter(mids)
        for mid in mids:
            assert mid in nids # make sure mids is a subset of nids
        cs = self.codes(nids=nids) # make a new codes object using the nids population
        nids2niis = cs.nids2niis
        notmids = [ nid for nid in nids if nid not in mids ] # nids not in mids
        # take product down all rows, only synchronous events across all mids cells will
        # survive, boolean array:
        mids_high = cs.c[nids2niis(mids)].prod(axis=0) == 1
        notmids_low = cs.c[nids2niis(notmids)].sum(axis=0) == 0 # boolean array
        # indices where mids are 1 and all the others are 0:
        i = (mids_high * notmids_low).nonzero()[0]
        return cs.t[i] # return the times at those indices

    def get_wordtsms(self, nids=None, mids=None):
        """Returns word times to the nearest msec, with the on bits specified in mids.
        nids lists the total population of neuron ids"""
        return np.int32(np.round(self.get_wordts(nids=nids, mids=mids) / 1e3))

    def get_intcodes(self, nids=None, shufflecodes=False):
        """Given neuron indices (ordered LSB to MSB top to bottom), returns an array of the
        integer representation of the neuronal population binary code for each time bin"""
        uns = get_ipython().user_ns
        assert uns['CODEKIND'] == 'binary'
        if nids == None:
            # randomly sample CODEWORDLEN bits of the nids
            nids = random.sample(self.cs.nids, uns['CODEWORDLEN'])
        return binarray2int(self.codes(nids=nids, shufflecodes=shufflecodes).c)

    def intcodesPDF(self, nids=None):
        """Returns the observed pdf across all possible population binary code words,
        labelled according to their integer representation"""
        uns = get_ipython().user_ns
        if nids == None:
            # randomly sample CODEWORDLEN bits of the nids
            nids = random.sample(self.cs.nids, uns['CODEWORDLEN'])
        intcodes = self.get_intcodes(nids=nids)
        nbits = len(nids)
        p, bins = pmf(intcodes, bins=np.arange(2**nbits))
        return p, bins

    def intcodesFPDF(self, nids=None):
        """the F stands for factorial. Returns the probability of getting each population
        binary code word, assuming independence between neurons, taking into account each
        neuron's spike (and no spike) probability"""
        uns = get_ipython().user_ns
        if nids == None:
            # randomly sample CODEWORDLEN bits of the nids
            nids = random.sample(self.cs.nids, uns['CODEWORDLEN'])
        nbits = len(nids)
        intcodes = np.arange(2**nbits)
        # this is like dict comprehension, pretty awesome!:
        #neurons = dict( (ni, self.neurons[ni]) for ni in nids )
        codes = self.codes(nids=nids)
        spikeps = [] # list spike probabilities for all neurons
        for neuroncode in codes.c: # for each neuron, ie each row
            # calc the average p of getting a spike for this neuron, within any time bin
            spikeps.append(neuroncode.mean())
        # convert to an nbits*1 array, make sure it's explicitly treated as a 2D array that
        # can be transposed, or something
        spikeps = np.array(spikeps, ndmin=2)
        nospikeps = 1 - spikeps
        #print 'spikesps: ', spikeps.__repr__()
        #print 'nospikesps: ', nospikeps.__repr__()
        binarytable = core.getbinarytable(nbits)
        # 2D array of probs of having a 1 in the right place for all possible population
        # code words:
        pon = binarytable * spikeps.transpose()
        # 2D array of probs of having a 0 in the right place for all possible population
        # code words:
        poff = (1 - binarytable) * nospikeps.transpose()
        #print 'pon', pon.__repr__()
        #print 'poff', poff.__repr__()
        # add the 2D arrays, each has zero p values where the other has non-zero p values:
        x = pon + poff
        #print 'x', x.__repr__()
        # take the product along the 0th axis (the columns) to get the prob of each
        # population code word
        intcodeps = x.prod(axis=0)
        return intcodeps, intcodes

    def ising(self, nids=None, R=None, algorithm='CG'):
        """Returns a maximum entropy Ising model that takes into account pairwise
        correlations within neuron codes. R = (R0, R1) torus. Algorithm can be 'CG', 'BFGS',
        'LBFGSB', 'Powell', or 'Nelder-Mead'"""
        uns = get_ipython().user_ns
        if nids == None:
            nids = self.cs.nids[0:uns['CODEWORDLEN']]
        #print 'nids:', nids.__repr__()
        if R:
            assert len(R) == 2 and R[0] < R[1] # should be R = (R0, R1) torus
        codes = self.codes(nids=nids)

        #c = codes.c
        # convert values in codes object from [0, 1] to [-1, 1] by mutliplying by 2 and
        # subtracting 1
        c = codes.c.copy() # don't modify the original
        c = c*2 - 1 # this should be safe to do cuz c is a 2D array of signed int8 values
        #print 'c:', c.__repr__()
        means = [ row.mean() for row in c ] # iterate over rows of codes in c
        nrows = c.shape[0]
        pairmeans = []
        for i in range(0, nrows):
            for j in range(i+1, nrows):
                if R == None or (R[0] < core.dist(self.r.n[nids[i]].pos,
                                                  self.r.n[nids[j]].pos) < R[1]):
                    # take a pair of rows, find the mean of their elementwise product:
                    pairmeans.append((c[i]*c[j]).mean())
                else:
                    # pair are outside the torus, ignore their pairmeans:
                    pairmeans.append(None)
        ising = core.Ising(means=means, pairmeans=pairmeans, algorithm=algorithm)
        return ising


class NetstateIsingHist(BaseNetstate):
    """Netstate Ising parameter histograms. See Schneidman 2006 Fig 3b"""
    def calc(self, ngroups=5, algorithm='CG'):
        """Collects hi and Jij parameter values computed from ising models
        of ngroups subgroups of cells of size nbits"""
        uns = get_ipython().user_ns
        self.nbits = uns['CODEWORDLEN']
        self.ngroups = ngroups
        self.algorithm = algorithm

        self.ims = [] # holds Ising Model objects
        self.his = []
        self.Jijs = []

        for groupi in range(self.ngroups): # for each group of nbits cells
            nids = random.sample(self.cs.nids, self.nbits) # randomly sample nbits of nids
            im = self.ising(nids=nids, algorithm=algorithm) # returns a maxent Ising model
            self.ims.append(im)
            self.his.append(im.hi)
            self.Jijs.append(im.Jij)
        return self

    def plot(self, nbins=50, hirange=(-2.5, 2.5), Jijrange=(-1.1, 1.1)):
        """Plots hi and Jij histograms in separate figures"""
        try: self.his, self.Jij
        except AttributeError: self.calc()

        # histogram them in linear space
        hibins  = np.linspace(start=hirange[0], stop=hirange[1], num=nbins, endpoint=True)
        Jijbins = np.linspace(start=Jijrange[0], stop=Jijrange[1], num=nbins, endpoint=True)
        nhi  = np.histogram(self.his, bins=hibins, density=True)[0]
        nJij = np.histogram(self.Jijs, bins=Jijbins, density=True)[0]

        # plot the hi histogram
        f1 = pl.figure()
        a1 = f1.add_subplot(111)
        a1.hold(True)
        a1.bar(left=hibins, height=nhi, width=hibins[1]-hibins[0], color='g', edgecolor='g')
        gcfm().window.setWindowTitle(lastcmd())
        a1.set_title('hi histogram\n%s, nbits=%d, ngroups=%d, algorithm=%s'
                     % (lastcmd(), self.nbits, self.ngroups, self.algorithm))
        a1.set_ylabel('probability density')
        a1.set_xlabel('hi')
        a1.set_xlim(hirange)

        # plot the Jij histogram
        f2 = pl.figure()
        a2 = f2.add_subplot(111)
        a2.hold(True)
        a2.bar(left=Jijbins, height=nJij, width=Jijbins[1]-Jijbins[0],
               color='m', edgecolor='m')
        gcfm().window.setWindowTitle(lastcmd())
        a2.set_title('Jij histogram\n%s, nbits=%d, ngroups=%d, algorithm=%s'
                     % (lastcmd(), self.nbits, self.ngroups, self.algorithm))
        a2.set_ylabel('probability density')
        a2.set_xlabel('Jij')
        a2.set_xlim(Jijrange)

        f1.tight_layout(pad=0.3) # crop figure to contents
        f2.tight_layout(pad=0.3) # crop figure to contents
        self.f = {1:f1, 2:f2}
        self.a = {1:a1, 2:a2}
        return self


class NetstateNspikingPMF(BaseNetstate):
    """Netstate PMF of number of cells spiking in the same bin. See 2006 Schneidman fig 1e"""
    def calc(self):
        """Calcs the PMF of observing n cells spiking in the same time bin,
        as well as the PMF for indep cells (shuffled codes)"""
        uns = get_ipython().user_ns
        nbits = uns['CODEWORDLEN']
        if self.nidswasNone:
            self.nids = random.sample(self.cs.nids, nbits) # randomly sample nbits of the nids
            self.nids.sort()
            self.nbits = nbits
        else:
            self.nids = self.cs.nids
            self.nbits = len(self.nids)
        self.words = {}
        self.nspiking = {}
        self.pnspiking = {}
        self.bins = {}
        bins = np.arange(self.nneurons+1)
        for shufflecodes in (False, True):
            self.words[shufflecodes] = self.get_intcodes(nids=self.nids,
                                                         shufflecodes=shufflecodes)
            # collect observances of the number of cells spiking for each pop code time bin.
            # Convert the word at each time bin to binary, count the number of 1s in it.
            # np.binary_repr() is a bit faster than using core.bin()
            self.nspiking[shufflecodes] = [ np.binary_repr(word).count('1')
                                            for word in self.words[shufflecodes] ]
            # want all probs to add to 1, not their area, so use pmf:
            self.pnspiking[shufflecodes], self.bins[shufflecodes] = (
                pmf(self.nspiking[shufflecodes], bins=bins))

        assert (self.bins[False] == self.bins[True]).all() # paranoid, just checking
        # since they're identical, get rid of the dict and just keep one:
        self.bins = self.bins[False]
        assert approx(self.pnspiking[False].sum(), 1.0), ('total observed probs: %f' %
                      self.pnspiking[False].sum())
        assert approx(self.pnspiking[True].sum(), 1.0), ('total indep probs: %f' %
                      self.pnspiking[True].sum())

        return self

    def plot(self, xlim=(-0.5, 15.5), ylim=(10**-6, 10**0)):
        """Plots nspikingPMF, for both observed and shuffled (forcing independence) codes"""
        try: self.pnspiking, self.bins
        except AttributeError: self.calc()

        f = pl.figure()
        a = f.add_subplot(111)
        a.hold(True)
        a.plot(self.bins, self.pnspiking[False], 'r.-')
        a.plot(self.bins, self.pnspiking[True], 'b.-')
        titlestr = ''#'PMF of observing n cells spiking in the same time bin'
        titlestr += '\n%s' % lastcmd()
        if self.nidswasNone:
            titlestr += '\nnids: %r' % self.nids
        a.set_title(titlestr)
        a.legend(('observed', 'indep (shuffled)'))
        a.set_yscale('log')
        if xlim:
            a.set_xlim(xlim)
        if ylim:
            a.set_ylim(ylim)
        gcfm().window.setWindowTitle(lastcmd())
        a.set_xlabel('number of spiking cells in a bin')
        a.set_ylabel('probability')

        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        self.a = a
        return self


class NetstateScatter(BaseNetstate):
    """Netstate scatter analysis object. See Schneidman Figures 1f and 2a"""
    def calc(self, model='both', R=None, shufflecodes=False, algorithm='CG'):
        """Calculates the expected probabilities, assuming a model in ['indep', 'ising',
        'both'], of all possible population codes vs their observed probabilities. R = (R0,
        R1) torus. self's nids are treated in LSB to MSB order"""
        uns = get_ipython().user_ns
        self.nbits = uns['CODEWORDLEN']
        self.model = model
        self.R = R
        if R:
            assert len(R) == 2 and R[0] < R[1] # should be R = (R0, R1) torus
        self.shufflecodes = shufflecodes
        self.algorithm = algorithm

        if self.nidswasNone: # nids weren't specified in __init__
            # randomly sample nbits of the nids:
            self.nids = random.sample(self.cs.nids, self.nbits)
            self.nids.sort()
        else: # nids were specified in __init__
            self.nids = self.cs.nids
            self.nbits = min(len(self.nids), self.nbits) # make sure nbits isn't > len(nids)

        self.intcodes = self.get_intcodes(nids=self.nids, shufflecodes=self.shufflecodes)
        self.pobserved, self.observedwords = pmf(self.intcodes, bins=np.arange(2**self.nbits))

        if self.model == 'indep':
            # expected, assuming independence:
            self.pexpected, self.expectedwords = self.intcodesFPDF(nids=self.nids)
        elif self.model == 'ising':
            # get a maxent Ising model:
            ising = self.ising(nids=self.nids, R=self.R, algorithm=self.algorithm)
            self.pexpected = ising.p # expected, assuming maxent Ising model
            self.expectedwords = ising.intsamplespace
        elif self.model == 'both':
            # get a maxent Ising model:
            ising = self.ising(nids=self.nids, R=self.R, algorithm=self.algorithm)
            self.pexpected = ising.p # expected, assuming maxent Ising model
            self.expectedwords = ising.intsamplespace
            # expected, assuming independence:
            self.pindepexpected = self.intcodesFPDF(nids=self.nids)[0]
        else:
            raise ValueError('Unknown model %r' % self.model)
        # make sure we're comparing apples to apples:
        assert (self.observedwords == self.expectedwords).all()
        return self

    def plot(self, figsize=(7.5, 6.5), model='both', scale='freq',
             xlim=(10**-4, 10**2), ylim=(10**-11, 10**2),
             yticks=(10**-11, 10**-9, 10**-7, 10**-5, 10**-3, 10**-1, 10**1),
             color=False):
        """Scatterplots the expected probabilities of all possible population codes (y axis)
        vs their observed probabilities (x axis). nids are in LSB to MSB order"""
        try: self.pobserved, self.pexpected
        except AttributeError: self.calc(model=model)

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        lo = min(xlim[0], ylim[0])
        hi = max(xlim[1], ylim[1])
        a.plot((lo, hi), (lo, hi), 'k-') # plot a y=x line
        a.hold(True)

        ## TODO: add legend instead of text colorguide
        ## TODO: make tooltips work again, old wx code disabled for now:
        # create a long tooltip with newline to get around bug where newlines aren't
        # recognized on subsequent self.tooltip.SetTip() calls
        #self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100))
        #self.tooltip.Enable(False) # leave disabled for now
        #self.tooltip.SetDelay(0) # set popup delay in ms
        #gcfm().canvas.SetToolTip(self.tooltip) # connect the tooltip to the canvas
        # connect the mpl event to the action:
        #f.canvas.mpl_connect('motion_notify_event', self._onmotion)

        # pylab.scatter(pobserved, pexpected), followed by setting the x and y axes to log
        # scale freezes the figure and runs 100% cpu:
        #gca().set_xscale('log')
        #gca().set_yscale('log')
        # use loglog() instead

        # colour each scatter point according to how many 1s are in the population code word
        # it represents. This is done very nastily, could use a cleanup:
        tres = get_ipython().user_ns['CODETRES']
        if scale == 'freq':
            norm = tres / 1e6 # convert scale to pattern freq in Hz
        elif scale == 'prob':
            norm = 1 # leave scale as pattern probabilities
        else:
            raise ValueError('Unknown scale %r' % scale)
        self.norm = norm

        if color:
            inds = []
            for nspikes in range(0, 5):
                inds.append([])
                [ inds[nspikes].append(i) for i in range(0, 2**self.nbits) if
                  core.bin(i).count('1') == nspikes ]
            # make local copies that are safe to modify for colour plotting and stuff
            pobserved = self.pobserved.copy()
            pexpected = self.pexpected.copy()
            pobserved1 = pobserved[inds[1]]; pexpected1 = pexpected[inds[1]]
            pobserved2 = pobserved[inds[2]]; pexpected2 = pexpected[inds[2]]
            pobserved3 = pobserved[inds[3]]; pexpected3 = pexpected[inds[3]]
            pobserved4 = pobserved[inds[4]]; pexpected4 = pexpected[inds[4]]
            pobserved[inds[1]], pexpected[inds[1]] = None, None # remove all these
            pobserved[inds[2]], pexpected[inds[2]] = None, None
            pobserved[inds[3]], pexpected[inds[3]] = None, None
            pobserved[inds[4]], pexpected[inds[4]] = None, None

        colorguide = ''
        if self.model == 'both': # plot the indep model too, and plot it first
            a.loglog(self.pobserved/norm, self.pindepexpected/norm, '.', color='blue', ms=10)
            colorguide = ('blue, red = indep, ising\n')
        # plot whichever model was specified
        if color:
            # plot what's left in black:
            a.loglog(pobserved/norm, pexpected/norm, '.', color='black')
            a.loglog(pobserved4/norm, pexpected4/norm, '.', color='magenta')
            a.loglog(pobserved3/norm, pexpected3/norm, '.', color='blue')
            a.loglog(pobserved2/norm, pexpected2/norm, '.', color=(0, 1, 0))
            a.loglog(pobserved1/norm, pexpected1/norm, '.', color='red')
            colorguide = '    red: 1 spike patterns\n' + \
                         '  green: 2 spike patterns\n' + \
                         '   blue: 3 spike patterns\n' + \
                         'magenta: 4 spike patterns\n' + \
                         '  black: other patterns  \n'
        else:
            a.loglog(self.pobserved/norm, self.pexpected/norm, '.', color='red', ms=10)
        '''
        a.plot(pobserved/norm, pexpected/norm, 'k.')
        '''
        gcfm().window.setWindowTitle(lastcmd())
        missingcodeis = (self.pobserved == 0).nonzero()[0]
        nmissing = len(missingcodeis)
        percentmissing = nmissing / float(2**self.nbits) * 100
        '''
        missingcodetext = ''
        if nmissing != 0:
            missingcodes = self.observedwords[missingcodeis]
            pexpectedmissing = self.pexpected[missingcodeis]
            maxpi = pexpectedmissing.argmax()
            maxp = pexpectedmissing[maxpi]
            maxpcode = self.expectedwords[missingcodeis[maxpi]]
            missingcodetext += ('\n nmissingcodes: %d, maxpmissingcode: (%r, pexpected=%.3g)'
                                % (nmissing, core.bin(maxpcode, minbits=self.nbits), maxp))
        '''
        a.set_title(lastcmd())
        if scale == 'freq':
            labelend = 'state frequency (Hz)'
        elif scale == 'prob':
            labelend = 'state probability'
        a.set_xlabel('observed ' + labelend)
        a.set_ylabel('predicted ' + labelend)
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        if yticks:
            a.set_yticks(yticks)
        if self.model =='both':
            DJSstring = ('(%.4f, %.4f)' % (core.DJS(self.pobserved, self.pindepexpected),
                         core.DJS(self.pobserved, self.pexpected)))
        else:
            DJSstring = '%.4f' % core.DJS(self.pobserved, self.pexpected)

        # add stuff to top left of plot:
        a.text(0.01, 0.99, 'nids = %s\n'
                           '%s\n'
                           'dt = %d min'
                           % (self.nids, self.r.name, intround(self.r.dtmin)),
                           transform = a.transAxes,
                           horizontalalignment='left',
                           verticalalignment='top')
        # add stuff to bottom right of plot:
        uns = get_ipython().user_ns
        a.text(0.99, 0.01, '%s'
                           'DJS = %s\n'
                           '%.1f%% missing\n'
                           'tres = %d ms\n'
                           'phase = %d deg\n'
                           'R = %r um\n'
                           'minrate = %.2f Hz'
                           % (colorguide, DJSstring, percentmissing, uns['CODETRES']//1000,
                              uns['CODEPHASE'], self.R, uns['MINRATE'],),
                           transform = a.transAxes,
                           horizontalalignment='right',
                           verticalalignment='bottom')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        self.a = a
        return self

    def _onmotion(self, event):
        """Called during mouse motion over scatterplot figure. Pops up the corresponding
        population code word and its int representation when hovering over a neuron scatter
        point"""
        tres = get_ipython().user_ns['CODETRES']
        if event.xdata != None and event.ydata != None: # if mouse is inside the axes
            i  = approx(event.xdata, self.pobserved/self.norm, rtol=1e-1, atol=0).nonzero()[0] # find for what indices (if any) xdata == pobserved
            ii = approx(event.ydata, self.pexpected[i]/self.norm, rtol=1e-1, atol=0).nonzero()[0] # for those above, find for what index (if any) ydata == pexpected
            codeis = i[ii]
            if codeis.size > 0:
                #tip += 'i: %r' % i
                #tip += '\nii: %r' % ii
                #tip += '\ncodeis: %r' % codeis
                intcodes = self.observedwords[codeis] # get the int rep for those indices from either self.observedwords[i] or self.expectedwords[i], doesn't matter which since they should be identical
                codes = [ core.bin(intcode, minbits=self.nbits) for intcode in intcodes ]
                tip =  'codes: %s' % repr(codes).replace('\'', '')
                tip += '\nintcodes: %r' % list(intcodes)
                activenids = [ list(np.asarray(self.nids)[::-1][charfind(code, '1')]) for code in codes ]
                tip += '\nactivenids: %r' % activenids
                tip += '\npattern counts: %r' % [ (self.intcodes == intcode).sum() for intcode in intcodes ]
                tip += '\npattern freqs (Hz): %s' % repr([ '%.3g' % (p / tres * 1e6) for p in self.pobserved[codeis] ]).replace('\'', '')
                tip += '\n(pobserved, pexpected): %s' % repr(zip([ '%.3g' % val for val in self.pobserved[codeis] ], [ '%.3g' % val for val in self.pexpected[codeis] ])).replace('\'', '')
                tip += '\npobserved / pexpected: %s' % repr([ '%.3g' % (o/float(e)) for o, e in zip(self.pobserved[codeis], self.pexpected[codeis]) ]).replace('\'', '')
                tip += '\npexpected / pobserved: %s' % repr([ '%.3g' % (e/float(o)) for o, e in zip(self.pobserved[codeis], self.pexpected[codeis]) ]).replace('\'', '')
                self.tooltip.SetTip(tip) # update the tooltip
                self.tooltip.Enable(True) # make sure it's enabled
            else:
                self.tooltip.Enable(False) # disable the tooltip
        else: # mouse is outside the axes
            self.tooltip.Enable(False) # disable the tooltip


class NetstateI2vsIN(BaseNetstate):
    """Netstate I2/IN vs IN (fraction of pairwise correlated entropy vs all correlated
    entropy) analysis. See Schneidman fig 2c"""
    def calc(self, N=10, ngroups=15):
        """Computes I2/IN vs IN, for ngroups of cells. This shows what fraction of
        network correlation is accounted for by the maxent pairwise model. Plotting
        Icond-indep is different from S1 (see below), and sounds annoying and not worth it
        (see methods in Schneidman 2006)"""
        self.N = N
        self.ngroups = ngroups

        self.nidss = nCrsamples(objects=self.neurons.keys(),
                                r=self.N, # pick N neurons at random
                                nsamples=self.ngroups) # do it ngroups times
        I2s = []
        INs = []
        tres = get_ipython().user_ns['CODETRES']
        for groupi, nids in enumerate(self.nidss):
            p1 = np.asarray(self.intcodesFPDF(nids=nids)[0]) # indep model
            p2 = self.ising(nids=nids).p # expected, assuming maxent Ising model
            pN = np.asarray(self.intcodesPDF(nids=nids)[0]) # observed word probs
            S1 = entropy_no_sing(p1) # ignore any singularities
            S2 = entropy_no_sing(p2)
            SN = entropy_no_sing(pN)
            IN = S1 - SN
            I2 = S1 - S2
            I2s.append(I2 / tres * 1e6) # convert to bits/sec
            INs.append(IN / tres * 1e6)
            print 'groupi',
        print('\n')
        self.I2s = np.asarray(I2s)
        self.INs = np.asarray(INs)
        self.I2divIN = self.I2s / self.INs

        return self

    def plot(self, xlim=(0.0, None), ylim=(0.0, 1.0)):
        """Plots I2/IN vs IN"""

        try: self.I2s, self.INs, self.I2divIN
        except AttributeError: self.calc()

        f = pl.figure()
        gcfm().window.setWindowTitle(lastcmd())
        a = f.add_subplot(111)
        a.plot(self.INs, self.I2divIN, 'r.')
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xlabel('IN (bits / sec)')
        a.set_ylabel('I2 / IN')
        a.set_title('%s' % lastcmd())
        # add mean and std to bottom right:
        a.text(0.99, 0.01, 'mean=%.3f, std=%.3f' %
               (self.I2divIN.mean(), self.I2divIN.std()),
               transform=a.transAxes,
               horizontalalignment='right',
               verticalalignment='bottom')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        self.a = a
        return self


class NetstateDJSHist(BaseNetstate):
    """Jensen-Shannon histogram analysis. See Schneidman 2006 figure 2b"""
    def calc(self, ngroups=5, models=['ising', 'indep'], R=None, shufflecodes=False,
             algorithm='CG'):
        """Calculates Jensen-Shannon divergences and their ratios
        for ngroups random groups of cells, each of length nbits. R = (R0, R1) torus"""
        uns = get_ipython().user_ns
        self.nbits = uns['CODEWORDLEN']
        self.ngroups = ngroups
        self.models = models
        self.R = R
        if R:
            assert len(R) == 2 and R[0] < R[1] # should be R = (R0, R1) torus
        self.shufflecodes = shufflecodes
        self.algorithm = algorithm

        self.nidss = [] # list of lists, each sublist is a group of neuron indices
        # hold the Jensen-Shannon divergences for different models and different groups of
        # neurons:
        self.DJSs = {}
        for model in  self.models:
            # init a dict with the model names as keys, and empty lists as values
            self.DJSs[model] = []
        for groupi in range(self.ngroups): # for each group of nbits cells
            # randomly sample nbits of the Netstate Codes' nids attrib:
            nids = random.sample(self.cs.nids, self.nbits)
            self.nidss.append(nids)
            for modeli, model in enumerate(self.models): # for each model, use the same nids
                so = NetstateScatter(recording=self.r, experiments=self.e, nids=nids)
                so.calc(model=model, R=self.R, shufflecodes=self.shufflecodes,
                        algorithm=self.algorithm)
                self.DJSs[model].append(core.DJS(so.pobserved, so.pexpected))
            print '%d' % groupi,
        print('\n')
        # now find the DJSratios between the two models, for each group of neurons
        # do it only if there's 2 models, otherwise it's indeterminate which two to take
        # ratio of
        if len(self.models) == 2:
            # 2nd model as a ratio of the 1st
            self.DJSratios = (np.asarray(self.DJSs.values()[1]) /
                              np.asarray(self.DJSs.values()[0]))

        return self

    def plot(self, ngroups=5, logrange=(-3.667, -0.333), nbins=50, publication=False):
        """Plots histogram DJSs and DJSratios in logspace"""
        try: self.nidss, self.DJSs
        except AttributeError: self.calc(ngroups=ngroups)
        x = np.logspace(start=logrange[0], stop=logrange[1], num=nbins, endpoint=True,
                        base=10.0)
        n = {} # stores a list of the bin heights in a separate key for each model
        for model in self.models:
            n[model] = np.histogram(self.DJSs[model], bins=x, density=False)[0]
        color = {'indep': 'blue', 'ising': 'red'} # maps from model name to colour
        # each bar will have a different width, convert to list so you can append
        barwidths = list(np.diff(x))
        # need to add one more entry to barwidth to the end to get nbins of them:
        barwidths.append(0) # don't display the last one
        logbinwidth = (logrange[1]-logrange[0]) / float(nbins)
        #barwidths.append(10**(logrange[1]+logbinwidth) - x[-1]) # should be exactly correct

        # plot DJSs of all models on the same axes
        f1 = pl.figure()
        a1 = f1.add_subplot(111)
        #a1.hold(True)
        bars = {}
        heights = {}
        for model in self.models:
            heights[model] = n[model] / float(self.ngroups * logbinwidth) # density
            bars[model] = a1.bar(left=x, height=heights[model], width=barwidths,
                                 color=color[model], edgecolor=color[model])
        # need to set scale of x axis AFTER bars have been plotted, otherwise
        # autoscale_view() call in bar() raises a ValueError for log scale:
        a1.set_xscale('log', basex=10)
        a1.set_xlim(xmin=10**logrange[0], xmax=10**logrange[1])
        gcfm().window.setWindowTitle(lastcmd())
        a1.set_title('%s' % lastcmd())
        if publication:
            a1.set_xticklabels(['', '0.001', '0.01', '0.1', '']) ## TODO: terrible hack!
            a1.set_ylim(ymin=0, ymax=4)
            a1.set_yticks((0, 2, 4))
            for label in a1.get_xticklabels():
                label.set_size(30)
            for label in a1.get_yticklabels():
                label.set_size(30)
            a1.legend([ bars[model][0] for model in self.models ],
                      ['pairwise', 'independent'], loc='upper right',
                      # grab first bar for each model, label it with model name:
                      prop=mpl.font_manager.FontProperties(size=20) )
        else:
            a1.set_xlabel('DJS (bits)')
            #a1.set_ylabel('number of groups of %d cells' % self.nbits)
            a1.set_ylabel('probability density (1 / log10(DJS))')
            a1.legend([ bars[model][0] for model in self.models ],
                      ['pairwise', 'independent'], loc='upper right')

        # add stuff to top left of plot:
        a1.text(0.01, 0.99, '%s\n'
                            'dt = %d min'
                            % (self.r.name, intround(self.r.dtmin)),
                            transform = a1.transAxes,
                            horizontalalignment='left',
                            verticalalignment='top')
        '''
        # plot DJSratios
        if len(self.models) == 2:
            f2 = pl.figure()
            a2 = f2.add_subplot(111)
            # bin heights for the DJSratios
            nratios = np.histogram(self.DJSratios, bins=x, density=False)[0]
            a2.bar(left=x, height=nratios, width=barwidths, color='g', edgecolor='g')
            # need to set scale of x axis AFTER bars have been plotted, otherwise
            # autoscale_view() call in bar() raises a ValueError for log scale
            a2.set_xscale('log', basex=10)
            gcfm().window.setWindowTitle(lastcmd())
            a2.set_title('Jensen-Shannon divergence ratios histogram\n%s' % lastcmd())
            a2.set_ylabel('number of groups of %d cells' % self.nbits)
            a2.set_xlabel('DJS ratio (%s / %s)' % (self.models[1], self.models[0]))
        '''
        f1.tight_layout(pad=0.3) # crop figure to contents
        #f2.tight_layout(pad=0.3) # crop figure to contents
        #self.f = {1:f1, 2:f2}
        #self.a = {1:a1, 2:a2}
        self.f = {1:f1}
        self.a = {1:a1}
        return self


class NetstateS1INvsN(BaseNetstate):
    """Analysis of uncorrelated entropy and reduction by correlated entropy for increasing
    network size N"""
    def calc(self, minN=4, maxN=15, maxnsamples=10):
        """Calculates the average independent (uncorrelated) cell entropy S1
        and average network multi-information IN (IN = S1 - SN) vs network size N.
        IN is how much the correlated entropy reduces the total entropy of the system.
        For each network size up to maxN, averages S1 and IN over maxnsamples (or less if
        that many aren't possible) number of groups at each value of N"""
        self.minN = minN
        self.maxN = maxN
        self.maxnsamples = maxnsamples

        self.S1ss = [] # as f'n of N
        self.INss = []
        self.N = range(self.minN, self.maxN+1) # network sizes from minN up to maxN
        #tstart = time.clock()
        # nsamples as a f'n of N. For each value of N, take up to maxnsamples of all the
        # other neurons, if that many are even possible
        self.nsamples = [ min(nCr(self.nneurons, r), self.maxnsamples) for r in self.N ]
        tres = get_ipython().user_ns['CODETRES']
        for ni, n in enumerate(self.N): # for all network sizes
            # get a list of lists of neuron indices
            nidss = nCrsamples(objects=self.neurons.keys(),
                               r=n, # pick n neurons
                               nsamples=self.nsamples[ni] ) # do it at most maxnsamples times
            S1s = []
            INs = []
            for nidsi, nids in enumerate(nidss):
                #t2 = time.clock()
                p1 = np.asarray(self.intcodesFPDF(nids=nids)[0]) # indep model
                pN = np.asarray(self.intcodesPDF(nids=nids)[0]) # observed word probs
                #print 'calcing ps took: %f sec' % (time.clock()-t2)
                S1 = entropy_no_sing(p1) # ignore any singularities
                SN = entropy_no_sing(pN)
                # better be, indep model assumes the least structure:
                assert S1 > SN or approx(S1, SN), 'S1 is %.20f, SN is %.20f' % (S1, SN)
                IN = S1 - SN
                #print S1, SN, IN
                S1s.append(S1 / tres * 1e6) # convert to bits/sec
                INs.append(IN / tres * 1e6)
            self.S1ss.append(S1s)
            self.INss.append(INs)
        print('\n')
        self.S1mean = [ np.asarray(S1s).mean() for S1s in self.S1ss ]
        self.S1std = [ np.asarray(S1s).std() for S1s in self.S1ss ]
        self.S1sem = np.asarray(self.S1std) / sqrt(np.asarray(self.nsamples))
        self.INmean = [ np.asarray(INs).mean() for INs in self.INss ]
        self.INstd = [ np.asarray(INs).std() for INs in self.INss ]
        self.INsem = np.asarray(self.INstd) / sqrt(np.asarray(self.nsamples))

        return self

    def plot(self, xlim=(1e0, 1e3), ylim=(1e-2, 1e4)):
        """Plots the average independent (uncorrelated) cell entropy S1
        and average network multi-information IN (IN = S1 - SN) vs network size N."""

        try: self.S1ss
        except AttributeError: self.calc()

        f = pl.figure()
        gcfm().window.setWindowTitle(lastcmd())
        a = f.add_subplot(111)
        a.hold(True)
        # plot all the samples before plotting the means with errorbars:
        for n, S1s in zip(self.N, self.S1ss):
            a.plot([n]*len(S1s), S1s, '_', markersize=4, color='lightblue')
        for n, INs in zip(self.N, self.INss):
            a.plot([n]*len(INs), INs, '_', markersize=4, color='pink')
        S1line = a.errorbar(self.N, self.S1mean, yerr=self.S1sem, fmt='b.')[0]
        INline = a.errorbar(self.N, self.INmean, yerr=self.INsem, fmt='r.')[0]
        # do least squares polynomial fit in log10 space
        # returns slope and y intercept:
        mS1, bS1 = sp.polyfit(log10(self.N), log10(self.S1mean), 1)
        mIN, bIN = sp.polyfit(log10(self.N), log10(self.INmean), 1)
        xintersect = (bIN - bS1) / (mS1 - mIN)
        x = np.array([-1, 3]) # define x in log10 space, this is really [0.1, 1000]
        self.yS1 = mS1*x + bS1 # y = mx + b
        self.yIN = mIN*x + bIN
        # take their power to make up for both the x and y scales being log
        a.plot(10.0**x, 10.0**self.yS1, 'b-')
        a.plot(10.0**x, 10.0**self.yIN, 'r-')
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xlabel('Number of cells')
        a.set_ylabel('bits / sec')
        a.set_title('S1 & IN vs N\n%s' % lastcmd())
        a.legend((S1line, INline), ('S1, slope=%.3f' % mS1, 'IN, slope=%.3f' % mIN),
                 loc='lower right')
        # add text box to upper right corner of axes
        a.text(0.99, 0.98, 'Nc=%d' % np.round(10**xintersect),
                           transform = a.transAxes,
                           horizontalalignment = 'right',
                           verticalalignment = 'top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        self.a = a
        return self


class NetstateNNplus1(BaseNetstate):
    """Analysis of amount of mutual information between N cells and the N+1th cell"""
    def calc(self, Nplus1s=None, maxN=15, maxnsamples=10):
        """Calculates Schneidman Figure 5b. Averages over as many as maxnsamples different
        groups of N cells for each N+1th cell in Nplus1s, all done for different values of N
        up to maxN"""
        # list of all indices of neurons that will be treated as the N+1th neuron:
        if Nplus1s == None:
            Nplus1s = self.cs.nids
        else:
            Nplus1s = toiter(Nplus1s)
        nNplus1s = len(Nplus1s)
        dims = (maxN, nNplus1s, maxnsamples)
        mask = np.zeros(dims) # this will be converted to an array of Falses
        # masked array that holds the mutual info between N and N+1th cells, as a ratio of
        # the N+1th cell's entropy. Index like: IdivS[ni, Nplus1i, samplei], ie group size,
        # N+1th cell you're comparing to, and number of samples of size N taken from the
        # possible combs:
        IdivS = np.ma.array(mask, mask=mask, fill_value=666)
        # cell group size, excluding the N+1th neuron. This will be the x axis in the plot:
        self.N = range(1, maxN+1)
        #self.N=[15]#.reverse() # for fun and pleasure
        # take up to maxnsamples of all the other neurons, if that many even exist (for the
        # lower N values, the constraint may end up being the total number of possible
        # combinations of cells), for each N+1th cell. Taken from nNplus1s-1 cuz you always
        # have to exclude an N+1th neurons:
        nsamples = [ min(maxnsamples, nCr(nNplus1s-1, r)) for r in self.N ]
        for ni, n in enumerate(self.N): # for all group sizes
            # mask out the sampleis that are out of range for this value of N, if any:
            IdivS.mask[ni, :, nsamples[ni]::] = True
        maximum = nNplus1s*sum(nsamples)

        # get the binary array for the whole population, then index into it appropriately in
        # the sample loop, find the corresponding integer codes, and feed it to MIbinarray,
        # so you don't have to unnecessarily re-generate it on every iteration:
        nids = self.cs.nids
        nids2niis = self.cs.nids2niis
        for ni, n in enumerate(self.N):
            for Nplus1i, Nplus1 in enumerate(Nplus1s): # for each N+1th neuron to compare to
                mii = nids2niis(Nplus1)
                nidscopy = copy(nids) # make a copy of neuron indices
                nidscopy.remove(Nplus1) # keep just the indices of all the other neurons
                # nsamples random unique choices of n items from nidscopy:
                samples = nCrsamples(nidscopy, n, nsamples[ni])
                # collect nsamples different combinations of the N other cells:
                for samplei, sample in enumerate(samples):
                    # most of the time (for n>1), sample will be a sequence of nids. Build
                    # an array of niis out of it to use as indices into the binary code
                    # array. Sometimes (for n=1) sample will be a scalar, hence the need to
                    # push it through toiter()
                    niis = np.array([ nids2niis(s) for s in toiter(sample) ])
                    IdivS[ni, Nplus1i, samplei] = (
                        MIbinarrays(Nbinarray=self.cs.c[niis],
                                    Mbinarray=self.cs.c[mii]).IdivS) # do it
        # reshape such that you collapse all Nplus1s and samples into a single dimension
        # (columns). The N are still in the rows:
        self.IdivS = IdivS.reshape(maxN, nNplus1s*maxnsamples)
        #logIdivS = log10(IdivS)
        # average over all Nplus1s and all samples. Values that are masked are ignored:
        self.IdivSmeans = self.IdivS.mean(axis=1)
        self.IdivSstds = self.IdivS.std(axis=1) # find stdev for the same
        assert self.IdivSmeans.shape == (maxN,)
        assert self.IdivSstds.shape == (maxN,)
        self.IdivSsems = self.IdivSstds / sqrt(np.asarray(nsamples)*nNplus1s)

        return self

    def plot(self, maxN=15, maxnsamples=10, xlim=(10**np.log10(0.9), 1e3), ylim=(1e-3, 1e1)):
        """Plots the figure with error bars"""
        try: self.IdivS
        except AttributeError: self.calc(maxN=maxN, maxnsamples=maxnsamples)

        f = pl.figure()
        gcfm().window.setWindowTitle(lastcmd())
        a = f.add_subplot(111)
        a.hold(True)
        for n, row in zip(self.N, self.IdivS): # underplot the samples for each value of N
            a.plot([n]*len(row), row, '_', markersize=4, color='deepskyblue')
        # plot the means and sems:
        a.errorbar(self.N, self.IdivSmeans, yerr=self.IdivSsems, fmt='b.')
        # do some linear regression in log10 space:
        m, b = sp.polyfit(log10(self.N), log10(self.IdivSmeans), 1) # slope and y intercept
        x = np.array([log10(0.9), 3]) # define x in log10 space, this is really [0.9, 1000]
        y = m*x + b
        xintersect = (0-b) / m # intersection point of regression line with y=1=10**0 line
        # raise them to the power to make up for the fact that both the x and y scales
        # will be log:
        plot(10.0**x, 10.0**y, 'b-')
        plot(10.0**x, [1e0]*2, 'r--') # plot horizontal line at y=1
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xlabel('Number of cells')
        a.set_ylabel('mutualinfo(N, N+1th) / entropy(N+1th)')
        a.set_title('fraction of info that N cells provide about the N+1th cell\n%s'
                    % lastcmd())
        # add text box to upper right corner of axes:
        a.text(0.99, 0.98, 'Nc=%d' % np.round(10**xintersect),
            transform = a.transAxes,
            horizontalalignment = 'right',
            verticalalignment = 'top')
        # add slope of fit line to bottom right:
        a.text(0.99, 0.01, 'slope=%.3f' % m,
            transform = a.transAxes,
            horizontalalignment = 'right',
            verticalalignment = 'bottom')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        self.a = a
        return self
        '''
        # plot the distributions of IdivS
        for ni, n in enumerate(self.N):
            f = pl.figure()
            gcfm().window.setWindowTitle('%s IdivS distrib for N=%d' % (lastcmd(), n))

            # indexes the non-masked entries in IdivS, for this ni:
            notmaskedis = self.IdivS[ni].mask==False

            a1 = f.add_subplot(211) # axes with linear bins
            heights, bins = np.histogram(self.IdivS[ni, notmaskedis], bins=arange(0, 1, 0.02))
            barwidth = bins[1]-bins[0]
            a1.bar(left=bins, height=heights, width=barwidth, bottom=0, color='k')
            #a1.set_xlabel('mutualinfo(N, N+1th) / entropy(N+1th)')
            a1.set_ylabel('count')
            a1.set_title('IdivS distrib for N=%d' % n)

            a2 = f.add_subplot(212) # axes with log bins
            start = log10(0.001)
            stop = log10(1)
            bins = np.logspace(start=start, stop=stop, num=50, endpoint=True, base=10.0)
            heights, bins = np.histogram(self.IdivS[ni, notmaskedis], bins=bins)
            # each bar will have a different width, convert to list so you can append:
            barwidth = list(diff(bins))
            # need to add one more entry to barwidth to the end to get nbins of them:
            #barwidth.append(barwidth[-1]) # not exactly correct
            logbinwidth = (log10(stop)-log10(stop)) / float(len(bins))
            barwidth.append(10**(log10(stop)+logbinwidth) - stop) # should be exactly correct
            a2.bar(left=bins, height=heights, width=barwidth, bottom=0, color='k')
            a2.set_xscale('log')
            a2.set_xlabel('mutualinfo(N, N+1th) / entropy(N+1th)')
            a2.set_ylabel('count')
        '''

class NetstateCheckcells(BaseNetstate):
    """Analysis of how activity rates of each cell in the population vary with
    the overall amount of activity in the rest of the population"""
    def _calc(self, ni=None, othernids=None, shufflecodes=False):
        """Calculates the joint pdf of cell ni activity and the number of cells in
        othernids being active at the same time. ni should not be in othernids"""
        assert ni not in othernids
        nids2niis = self.cs.nids2niis
        nii = nids2niis(ni)
        otherniis = nids2niis(othernids)
        nothers = len(othernids)

        # 0s and 1s, this picks out the row in the binary code array that corresponds to ni:
        nicode = self.cs.c[nii]
        othercodes = self.cs.c[otherniis]
        if shufflecodes:
            nicode = np.asarray(shuffle(nicode))
            othercodes = np.asarray(shuffle(othercodes))
        nothersactive = othercodes.sum(axis=0) # anywhere from 0 up to and including nothers

        # build up joint pdf of the nicode and nothersactive:
        xedges = np.array([0, 1, 2]) # 2 is needed as the rightmost bin edge for histogram2d
        # anywhere from 0 up to and including nothers, plus nothers+1 as the
        # rightmost bin edge:
        yedges = np.arange(nothers+2)
        bins = [xedges, yedges]
        # generate joint pdf, nicode are in the rows, nothersactive are in the columns,
        # leave it unnormalized, just counts:
        jpdf, xedgesout, yedgesout = np.histogram2d(nicode, nothersactive, bins,
                                                    density=False)
        # now, normalize each column separately, so that say, for nothersactive==5,
        # p(checkcell==0)+p(checkcell==1) == 1.0
        jpdf = np.float64(jpdf) # convert to floats, updated entries are trunc'd to ints
        for coli in range(jpdf.shape[-1]):
            # save the normalized column back to the jpdf
            jpdf[:, coli] = normalize(jpdf[:, coli])
        return jpdf

    def calc(self, nids=None, othernids=None, nothers=None, nsamples=10, shufflecodes=False):
        """Calcs the probability of each cell (in nids) being active vs. the number of
        other active cells (in the Recording) at that time. For each ni, calcs an average over
        nsamples, each being a different sample of nothers from othernids.
        See Schneidman figure 5c"""
        if nids == None:
            self.nids = self.cs.nids
        else:
            self.nids = toiter(nids)
        if othernids == None:
            self.othernids = self.cs.nids
        else:
            self.othernids = toiter(othernids)
        if nothers == None:
            less = 0
            while nCr(len(self.othernids), len(self.othernids)-less) < nsamples:
                less += 1
            # -1 to remove ni, -less again to allow for at least nsamples combos of othernids:
            nothers = len(self.othernids) - 1 - less
        self.N = np.arange(nothers+1)

        try: self.jpdfss
        except AttributeError:
            # init dicts to store jpdfs and other stuff in
            self.jpdfss = {}
            self.jpdfmeans = {}
            self.jpdfstds = {}
            self.jpdfsems = {}

        for ni in self.nids:
            try:
                self.jpdfss[ni]
            except KeyError:
                othernids = copy(self.othernids) # don't modify the original
                try:
                    # all the other possible nids, excluding the current ni
                    othernids.remove(ni)
                except ValueError: # ni isn't in othernids, nothing to remove
                    pass
                # get nsamples unique random samples of length nothers from othernids:
                othernidss = nCrsamples(objects=othernids, r=nothers, nsamples=nsamples)
                jpdfs = []
                for othernids in othernidss: # collect jpdfs across all random samples
                    jpdf = self._calc(ni=ni, othernids=othernids, shufflecodes=shufflecodes)
                    jpdfs.append(jpdf)
                jpdfs = np.asarray(jpdfs) # this is an nsamples x 2 x (1+nothers) matrix
                self.jpdfss[ni] = jpdfs
                # find the mean jpdf across all nsamples jpdfs:
                self.jpdfmeans[ni] = jpdfs.mean(axis=0)
                # find the stdev across all nsamples jpdfs:
                self.jpdfstds[ni] = jpdfs.std(axis=0)
                self.jpdfsems[ni] = self.jpdfstds[ni] / sqrt(nsamples)
        return self

    def plot(self, nids=None, nothers=None, nsamples=10):
        """Plots the desired neurons so you can see if they behave like check cells"""
        try: self.jpdfss
        except AttributeError: self.calc(nids=nids, nothers=nothers, nsamples=nsamples)
        try: self.f
        except AttributeError: self.f = {}
        try: self.a
        except AttributeError: self.a = {}
        if nids == None:
            nids = self.nids
        else:
            nids = toiter(nids)
        for ni in nids:
            f = pl.figure()
            gcfm().window.setWindowTitle('%s for ni=%d' % (lastcmd(), ni))
            a = f.add_subplot(111)
            a.hold(True)
            # plot all the samples first
            for jpdf in self.jpdfss[ni]: # iter over the hyperrows
                # marginal pdf of getting a 1 for the check cell:
                a.plot(self.N, jpdf[1], '_', markersize=4, color='grey')
            # plot the stdevs, means, and sems of marginal pdf of getting a 1 for the
            # check cell:
            #a.errorbar(self.N, self.jpdfmeans[ni][1], yerr=self.jpdfstds[ni][1],
            #           fmt=None, capsize=0, ecolor='grey')
            a.errorbar(self.N, self.jpdfmeans[ni][1], yerr=self.jpdfsems[ni][1], fmt='k.-')
            a.set_ylim(ymin=0, ymax=1)

            titlestr = '%s\nni=%d' % (lastcmd(), ni)
            titlestr += ', nsamples=%d' % nsamples
            a.set_title(titlestr)
            a.set_xlabel('Number of other active cells')
            a.set_ylabel('Probability of cell ni being active')

            f.tight_layout(pad=0.3) # crop figure to contents
            self.f[ni] = f
            self.a[ni] = a
        return self


class NetstateTriggeredAverage(BaseNetstate):
    """Analysis that reverse correlates the occurence of a specific netstate
    to the stimuli in the Experiments in this Recording to build up a netstate
    triggered average"""
    def cut(self, trange):
        """Cuts network state word times according to trange"""
        lo, hi = self.wordts.searchsorted([trange[0], trange[1]]) # returns indices where tstart and tend would fit in wordts
        if trange[1] == self.wordts[min(hi, len(self.wordts)-1)]: # if tend matches a word time (protect from going out of index bounds when checking)
            hi += 1 # inc to include a word time if it happens to exactly equal tend. This gives us end inclusion
            hi = min(hi, len(self.wordts)) # limit hi to max slice index (==max value index + 1)
        cutwordts = self.wordts[lo:hi] # slice it
        return cutwordts

    def calc(self, intcode=None, nt=9, ti0=-4):
        """Calculate the network state triggered average for word intcode, using nt revcorr timepoints,
        starting at revcorr timepoint index ti0.

        For now, this uses the Codes object created across the entire Recording
        """
        self.intcode = intcode
        self.nt = nt # number of revcorr timepoints
        self.ti0 = ti0
        self.tis = range(ti0, ti0+nt, 1) # revcorr timepoint indices, can be -ve. these will be multiplied by the movie frame time
        words = binarray2int(self.cs.c)
        i = (words == intcode)
        assert i.any(), 'netstate intcode %d never occured'
        self.wordts = self.cs.t[i] # netstate intcode word times

        if self.e == None:
            self.e = self.r.e # if no specific experiments were specified to revcorr to in __init__, revcorr to all of them
        self.frames = {} # accumulates movie frames separately for each timepoint across all Experiments' movies
        for ti in self.tis:
            self.frames[ti] = []
        tstart = time.clock()
        pd = wx.ProgressDialog(title='NSTA progress: loading movies, collecting frames', message='',
                               maximum=len(self.e), style=1) # create a progress dialog
        for ei, e in enumerate(self.e.values()):
            cont, skip = pd.Update(ei-1, newmsg='experiment %d\nelapsed: %.1fs' % (e.id,
                                   time.clock()-tstart))
            if not cont:
                pd.Destroy()
                return
            # for now, we're using the Codes object created across the entire Recording. It might be slightly more correct to generate a separate codes object for each Experiment. That way, the wordts for would be aligned to the start of each Experiment, as opposed to the start of the Recording, as they are now.
            cutwordts = self.cut(e.trange) # get wordts that were active during this experiment. This isn't really necessary is it??????????????????????????? Might speed things up for the next searchsorted call, since cutwordts is shorter than wordts
            rcdini = e.din[:, 0].searchsorted(cutwordts) - 1 # revcorr dini. Find where the cutwordts times fall in the din, dec so you get indices that point to the most recent din value for each cutwordt
            #self.din = e.din[rcdini, 1] # get the din (frame indices) at the rcdini
            # for now, only do revcorr if experiment.stims has only one entry. Stims no longer exists in dimstim >= 0.16 anyway
            assert len(e.stims) == 1
            movie = e.stims[0] # this will have to change for dimstim >= 0.16
            movie.load() # ensure movie is loaded
            data = movie.data
            try:
                self.width
            except AttributeError: # init stuff
                self.width = movie.data.shape[-1]
                self.height = movie.data.shape[-2] # dims are nframes, height, width
                self.sweeptimeMsec = movie.sweeptimeMsec # changes in dimstim >= 0.16
                self.ndinperframe = intround(movie.sweeptimeMsec / float(e.REFRESHTIME / 1000.)) # changes in dimstim >= 0.16
                self.regionwidthDeg = movie.regionwidthDeg
                self.regionheightDeg = movie.regionheightDeg
                self.origDeg = e.origDeg
            # assert that all movies are the same size and in the same spot. that way you can just accumulate frames in a single array
            assert self.width == movie.data.shape[-1] # dims are nframes, height, width
            assert self.height == movie.data.shape[-2]
            assert self.sweeptimeMsec == movie.sweeptimeMsec # changes in dimstim >= 0.16
            assert self.regionwidthDeg == movie.regionwidthDeg
            assert self.regionheightDeg == movie.regionheightDeg
            assert self.origDeg == e.origDeg
            for ti in self.tis:
                shiftedrcdini = rcdini - ti*self.ndinperframe # this can unintentionally introduce -ve valued indices at the left boundary, or out of range values at right boundary
                shiftedrcdini = shiftedrcdini[shiftedrcdini >= 0] # remove any -ve valued indices. Is this the most efficient way to do this?
                shiftedrcdini = shiftedrcdini[shiftedrcdini <= len(e.din)-1] # remove any out of range values
                frameis = e.din[shiftedrcdini, 1] # get the din values (frame indices) at the rcdini for this timepoint
                # in Cat 15, we erroneously duplicated the first frame of the mseq movies at the end, giving us one more frame (0 to 65535 for mseq32) than we should have had (0 to 65534 for mseq32). We're now using the correct movies, but the din for Cat 15 mseq experiments still have those erroneous frame indices (65535 and 16383 for mseq32 and mseq16 respectively), so we'll just ignore them for revcorr purposes.
                if movie.oname == 'mseq32':
                    frameis = frameis[frameis != 65535] # remove all occurences of 65535
                elif movie.oname == 'mseq16':
                    frameis = frameis[frameis != 16383] # remove all occurences of 16383
                frames = data.take(frameis.astype(np.int32), axis=0) # collect the relevant frames for this timepoint, take is much faster than direct indexing, but have to typecast indices to int32, maybe cuz this machine is 32bit?
                self.frames[ti].append(frames)
            pd.Destroy()


        # now that we're outside the experiment loop, get the mean for each timepoint across all experiments
        self.rf = zeros([self.nt, self.height, self.width], dtype=np.float64) # 3D matrix to store the NSTA at each timepoint. rf == 'receptive field'
        self.t = [ ti*intround(self.sweeptimeMsec) for ti in self.tis ]
        pd = wx.ProgressDialog(title='NSTA progress: averaging frames', message='',
                               maximum=self.nt, style=1) # create a progress dialog
        for tii, ti in enumerate(self.tis):
            cont, skip = pd.Update(tii-1, newmsg='timepoint: %dms\nelapsed: %.1fs' % (self.t[tii], time.clock()-tstart))
            if not cont:
                pd.Destroy()
                return
            self.frames[ti] = np.concatenate(tuple(self.frames[ti])) # need to concatenate all lists for this ti into a single array
            self.rf[tii] = mean_accum(self.frames[ti]) # this is much faster than frames.mean()
            #self.rf[ti] = mean_accum2(data, frameis)
        pd.Destroy()
        return self

    def plot(self, intcode=None, nt=10, ti0=-4, interp='nearest', normed=True, scale=2.0):
        """Plots the spatiotemporal RF as bitmaps in a wx.Frame"""
        try:
            self.rf
        except AttributeError:
            self.calc(intcode=intcode, nt=nt, ti0=ti0)
        rf = self.rf.copy() # create a copy to manipulate for display purposes, (nt, width, height)
        if normed: # normalize across the timepoints for this RevCorr
            norm = mpl.colors.normalize(vmin=rf.min(), vmax=rf.max(), clip=True) # create a single normalization object to map luminance to the range [0,1]
            rf = norm(rf) # normalize the rf the same way across all timepoints
        else: # don't normalize across timepoints, leave each one to autoscale
            for ti in range(self.nt):
                norm = mpl.colors.normalize(vmin=None, vmax=None, clip=True) # create a normalization object to map luminance to the range [0,1], autoscale
                rf[ti] = norm(rf[ti]) # normalize the rf separately at each timepoint
        cmap = mpl.cm.jet # get a colormap object
        rf = cmap(rf)[::, ::, ::, 0:3] # convert luminance to RGB via the colormap, throw away alpha channel (not used for now in ReceptiveFieldFrame)
        rf = rf * 255 # scale up to 8 bit values
        rf = rf.round().astype(np.uint8) # downcast from float to uint8 for feeding to ReceptiveFieldFrame
        self.rfframe = NetstateReceptiveFieldFrame(title=lastcmd(), rfs=[rf], intcodes=self.intcode, t=self.t, scale=scale)
        self.rfframe.Show()
        return self


class RecordingNetstate(BaseRecording):
    """Mix-in class that defines Netstate related Recording methods"""
    def ns_(self, experiments=None, nids=None):
        """Returns a BaseNetstate object"""
        return BaseNetstate(recording=self, experiments=experiments, nids=nids)

    def ns_isinghist(self, experiments=None, nids=None):
        """Returns a NetstateIsingHist object"""
        return NetstateIsingHist(recording=self, experiments=experiments)

    def ns_nspikingpmf(self, experiments=None, nids=None):
        """Returns a NetstateNspikingPMF object"""
        return NetstateNspikingPMF(recording=self, experiments=experiments, nids=nids)

    def ns_scatter(self, experiments=None, nids=None, R=None):
        """Returns a NetstateScatter object"""
        ns_so = NetstateScatter(recording=self, experiments=experiments, nids=nids)
        ns_so.calc(R=R)
        return ns_so

    def ns_i2vsin(self, experiments=None, nids=None):
        """Returns a NetstateI2vsIN object"""
        return NetstateI2vsIN(recording=self, experiments=experiments, nids=nids)

    def ns_djshist(self, experiments=None, nids=None, ngroups=5, R=None):
        """Returns a NetstateDJSHist object"""
        ns_djso = NetstateDJSHist(recording=self, experiments=experiments, nids=nids)
        ns_djso.calc(ngroups=ngroups, R=R)
        return ns_djso

    def ns_s1invsn(self, experiments=None, nids=None):
        """Returns a NetstateS1INvsN object"""
        return NetstateS1INvsN(recording=self, experiments=experiments, nids=nids)

    def ns_nnplus1(self, experiments=None, nids=None):
        """Returns a NetstateNNplus1 object"""
        return NetstateNNplus1(recording=self, experiments=experiments, nids=nids)

    def ns_checkcells(self, experiments=None, nids=None):
        """Returns a NetstateCheckcells object"""
        return NetstateCheckcells(recording=self, experiments=experiments, nids=nids)

    def ns_ta(self, experiments=None, nids=None):
        """Returns a NetstateTriggeredAverage object"""
        return NetstateTriggeredAverage(recording=self, experiments=experiments, nids=nids)


class Recording(RecordingRevCorr,
                RecordingRaster,
                RecordingCode,
                RecordingNetstate,
                BaseRecording):
    """Inherits all the Recording classes into a single Recording class"""
    pass
