"""Defines the Track class"""

from __future__ import division

import os
import StringIO

import numpy as np

import pylab as pl
from pylab import get_current_fig_manager as gcfm

import core
from core import dictattr, TAB, td2usec, lastcmd, intround
from recording import Recording
from sort import TrackSort


class Track(object):
    """A track can have multiple recordings"""
    def __init__(self, path, animal=None):
        self.level = 2 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # string buffer to print tree hierarchy to
        self.path = path
        self.animal = animal
        if animal != None:
            # update parent animal's track dict, in case self wasn't loaded by its parent
            animal.tr[self.id] = self
        self.r = dictattr() # store recordings in a dictionary with attrib access

    def get_name(self):
        return os.path.split(self.path)[-1]
    
    name = property(get_name)

    def get_absname(self):
        """Absolute name including parent Animal, kind of like absolute path, but more
        abbreviated, as one would enter it at the IPython prompt"""
        return '.'.join((self.animal.name, self.name))

    absname = property(get_absname)

    def get_id(self):
        # replace first occurrence of 'tr' with nothing, keep the rest:
        id = self.name.lower().replace('tr', '', 1)
        if not id:
            raise NameError('Badly formatted track name: %s' % self.name)
        '''
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, as in '7c', leave it as a string
        '''
        return id
        
    id = property(get_id)

    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        if self.animal != None:
            self.animal.writetree(string)

    def load(self):
        treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)
        # collect recording names: 1st char of each name must be a digit, that's all:
        rnames = [ name for name in os.listdir(self.path)
                   if os.path.isdir(os.path.join(self.path, name))
                   and name[0].isdigit() ]
        rnames.sort() # alphabetical order
        dt = 0 # calculate total track duration by summing durations of all recordings
        for rname in rnames:
            path = os.path.join(self.path, rname)
            recording = Recording(path, track=self)
            recording.load()
            self.r[recording.id] = recording
            self.__setattr__('r' + str(recording.id), recording) # add shortcut attrib
            dt += recording.dt
        self.rnames = rnames # easy way to print out all recording names
        self.dt = dt
        self.dtsec = self.dt / 1e6
        self.dtmin = self.dtsec / 60
        self.dthour = self.dtmin / 60

        # create a TrackSort with TrackNeurons:
        self.sort = TrackSort(self)
        self.sort.load()
        # one way of calculating self.trange:
        #tranges = np.asarray([ n.trange for n in self.alln.values() ])
        #self.trange = min(tranges[:, 0]), max(tranges[:, 1])
        # better way of calculating self.trange:
        rids = sorted(self.r.keys()) # all recording ids in self
        r0 = self.r[rids[0]]
        r1 = self.r[rids[-1]]
        assert r0.datetime == self.datetime
        self.trange = r0.td+r0.trange[0], r1.td+r1.trange[1]

        self.calc_meanrates()

        # pttype better be the same for all member recordings:
        pttype = self.r[rids[0]].pttype # init to pttype of first recording
        for rid in rids[1:]:
            r = self.r[rid]
            # if recording doesn't have a pttype, it's probably from an old .spk file,
            # so don't bother doing this test:
            if hasattr(r, 'pttype') and pttype != r.pttype:
                raise ValueError("inconsistent polytrode types %r and %r in track %s"
                                 % (pttype, r.pttype, self.id))

    def calc_meanrates(self):
        """Calculate mean firing rates of all TrackNeurons in this track"""
        TRACKNEURONPERIOD = get_ipython().user_ns['TRACKNEURONPERIOD']
        if TRACKNEURONPERIOD == 'track':
            # calc tn.meanrate using entire track duration:
            for tn in self.alln.values():
                tn.meanrate = tn.nspikes / self.dtsec
        elif TRACKNEURONPERIOD == 'trange':
            # calc tn.meanrate using duration between its first and last spike:
            for tn in self.alln.values():
                if tn.dtsec == 0:
                    tn.meanrate = 0.0
                else:
                    tn.meanrate = tn.nspikes / tn.dtsec
        else:
            raise ValueError("invalid value for TRACKNEURONPERIOD: %r" % TRACKNEURONPERIOD)

    def get_meanrates(self):
        """Return mean firing rates of all TrackNeurons in this track"""
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

    def get_nids(self, rids=None):
        """Return nids of active neurons common to all recordings specified in rids.
        Otherwise, return all active nids in all recordings. Active neurons in a recording
        are those with at least MINRATE mean spike rate during the recording"""
        if rids == None: # return all nids in all recordings
            rids = list(self.r.keys())
            return np.unique(np.hstack([ self.r[rid].n.keys() for rid in rids ]))
        else: # return intersection of nids of specified recordings
            nids = [ self.r[rid].n.keys() for rid in rids ]
            return core.intersect1d(nids, assume_unique=True)

    def get_allnids(self, rids=None):
        """Return nids of all neurons (active and quiet) common to all recordings
        specified in rids, ie return the intersection. If rids==None, return the union
        of all nids in the track instead"""
        if rids == None:
            rids = sorted(self.r.keys())
            allnids = np.hstack([ self.r[rid].alln.keys() for rid in rids ])
            return np.unique(allnids)
        else:
            allnids = [ self.r[rid].alln.keys() for rid in rids ]
            return core.intersect1d(allnids, assume_unique=True)

    def pospdfrec(self, neurons=None, rids=None, dim='y', nbins=10, a=None, figsize=(7.5, 6.5)):
        """Plot PDF of cell positions ('x' or 'y') along the polytrode
        for each recording to get an idea of how cells are distributed in space,
        and how that changes from one recording to the next"""
        if a == None:
            f = pl.figure(figsize=figsize)
            a = f.add_subplot(111)
        else: # add to existing axes
            a.hold(True)
            f = pl.gcf()

        if rids == None:
            rids = self.r.keys()
            rids.sort()
        for rid in rids:
            r = self.r[rid]
            r.pospdf(neurons=neurons, dim=dim, nbins=nbins, a=a, stats=False, figsize=figsize)

        return a

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

        # use CLUSTERCOLOURDICT for familiarity with len 10 1-based id to colour mapping
        #color = CLUSTERCOLOURDICT[int(self.id)]
        color = 'k'

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
            a.text(0.99, 0.99, 'mean = %.3f\n'
                               'median = %.3f\n'
                               'mode = %.3f\n'
                               'stdev = %.3f\n'
                               'minrate = %.2f Hz\n'
                               'nneurons = %d\n'
                               'dt = %d min'
                               % (mean, median, mode, stdev,
                                  uns['MINRATE'], self.nneurons, intround(self.dtmin)),
                               transform = a.transAxes,
                               horizontalalignment='right',
                               verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        f.canvas.draw() # this is needed if a != None when passed as arg
        return a

    def scstim(self, method='weighted mean', width=None, tres=None, figsize=(7.5, 6.5)):
        """Scatter plot some summary statistic of spike correlations of each recording vs what
        stimulus group each recording falls into. width and tres dictate tranges to split
        recordings up into, if any"""

        ## TODO: for each pair of recordings, find common subset of active neurons and calculate
        ## pairwise corrs for each recording in that pair using just those neurons

        ## TODO: maybe limit to visually responsive cells

        uns = get_ipython().user_ns
        blankmseqrids = uns['BLANKMSEQRIDS'][self.absname]
        movdriftrids = uns['MOVDRIFTRIDS'][self.absname]

        blankmseqcorrs = []
        movdriftcorrs = []
        for rid in (blankmseqrids + movdriftrids):
            r = self.r[rid]
            print('%s: %s' % (r.absname, r.name))
            spikecorr = r.sc(width=width, tres=tres)
            sc = spikecorr.sct(method=method)[0]
            sc = sc[0] # pull out the spike correlation values that span all laminae
            if rid in blankmseqrids:
                blankmseqcorrs.append(sc)
            else:
                movdriftcorrs.append(sc)
        blankmseqcorrs = np.hstack(blankmseqcorrs)
        movdriftcorrs = np.hstack(movdriftcorrs)
        # repeat each element in blankmseqcorrs len(movdriftcorrs) times:
        x = np.repeat(blankmseqcorrs, len(movdriftcorrs))
        # tile movdriftcorrs len(blankmseqcorrs) times:
        y = np.tile(movdriftcorrs, len(blankmseqcorrs))

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        lim = min([x.min(), y.min(), 0]), max([x.max(), y.max()])
        a.plot(lim, lim, c='e', ls='--', marker=None) # y=x line
        a.plot(x, y, 'k.')
        #a.set_xlim(lim)
        #a.set_ylim(lim)
        a.set_xlabel('%s spike correlations: blankscreen and mseq' % method)
        a.set_ylabel('%s spike correlations: movie and drift bar' % method)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents
        f.show()

    def scsistim(self, method='weighted mean', width=None, tres=None, figsize=(7.5, 6.5)):
        """Scatter plot some summary statistic of spike correlations of each recording vs
        synchrony index SI. Colour each point according to stimulus type. width and tres
        dictate tranges to split recordings up into, if any"""

        ## TODO: for each pair of recordings, find common subset of active neurons and calculate
        ## pairwise corrs for each recording in that pair using just those neurons

        ## TODO: maybe limit to visually responsive cells

        uns = get_ipython().user_ns
        if width == None:
            width = uns['SIWIDTH'] # want powers of two for efficient FFT
        if tres == None:
            tres = width
        blankmseqrids = uns['BLANKMSEQRIDS'][self.absname]
        movdriftrids = uns['MOVDRIFTRIDS'][self.absname]

        blankmseq_scs, blankmseq_sis = [], []
        movdrift_scs, movdrift_sis = [], []
        for rid in (blankmseqrids + movdriftrids):
            r = self.r[rid]
            print('%s: %s' % (r.absname, r.name))
            spikecorr = r.sc(width=width, tres=tres)
            sc, si = spikecorr.si(method=method, plot=False) # calls sc.calc() and sc.si()
            sc = sc[0] # pull out the spike correlation values that span all laminae
            if rid in blankmseqrids:
                blankmseq_scs.append(sc)
                blankmseq_sis.append(si)
            else:
                movdrift_scs.append(sc)
                movdrift_sis.append(si)

        blankmseq_scs = np.hstack(blankmseq_scs)
        blankmseq_sis = np.hstack(blankmseq_sis)
        movdrift_scs = np.hstack(movdrift_scs)
        movdrift_sis = np.hstack(movdrift_sis)
        
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.plot(blankmseq_scs, blankmseq_sis, 'k.')
        a.plot(movdrift_scs, movdrift_sis, 'r.')
        a.set_ylim(0, 1)
        a.set_xlabel('%s spike correlations' % method)
        a.set_ylabel('synchrony index')
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents
