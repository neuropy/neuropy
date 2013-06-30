"""Defines the Sort class"""

from __future__ import division

import os
import StringIO
import datetime

import numpy as np

import core
from core import dictattr, rstrip, eof, TAB, PTCSHeader, SPKHeader, EPOCH, td2usec
from neuron import Neuron, TrackNeuron


class Sort(object):
    """A sort is a single spike extraction. Generally, there is one sort per recording,
    and sorts of the same name within the same track were extracted in the same spike
    sorting session"""
    def __init__(self, path, id=None, recording=None):
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.path = path
        self.id = id
        self.r = recording
        self.alln = {} # dict to store all Neurons

    def get_n(self):
        """Return dict of neurons that meet MINRATE"""
        n = {}
        MINRATE = get_ipython().user_ns['MINRATE']
        for neuron in self.alln.values():
            if neuron.meanrate >= MINRATE:
                n[neuron.id] = neuron
        return n

    n = property(get_n)

    def get_qn(self):
        """Return dict of quiet neurons, ie those that fail to meet MINRATE"""
        qn = {}
        MINRATE = get_ipython().user_ns['MINRATE']
        for neuron in self.alln.values():
            if neuron.meanrate < MINRATE:
                qn[neuron.id] = neuron
        return qn

    qn = property(get_qn)

    name = property(lambda self: os.path.split(self.path)[-1])
    nneurons = property(lambda self: len(self.n))
    nqneurons = property(lambda self: len(self.qn))
    nallneurons = property(lambda self: len(self.alln))
    nspikes = property(lambda self: self.header.nspikes)
    # .ptcs specific properties:
    # datetime object, calculated from header.datetime days since EPOCH"""
    datetime = property(lambda self: EPOCH + datetime.timedelta(days=self.header.datetime))
    pttype = property(lambda self: self.header.pttype)
    chanpos = property(lambda self: self.header.chanpos)

    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)
    
    def load(self):
        treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)
        
        if os.path.isfile(self.path):
            ext = os.path.splitext(self.path)[1]
            assert ext == '.ptcs'
            # it's a single .ptcs file
            self.loadptcs()
        elif os.path.isdir(self.path):
            # it's a directory of .spk files
            self.loadspk()
        else:
            raise RuntimeError

    def loadptcs(self):
        """Load neurons from a single .ptcs file"""
        self.header = PTCSHeader()
        with open(self.path, 'rb') as f:
            self.header.read(f)
            for i in range(self.header.nneurons):
                neuron = Neuron(self.path, sort=self)
                neuron.loadptcs(f, self.header)
                self.alln[neuron.id] = neuron # save it
            assert eof(f), 'File %s has unexpected length' % self.path

    def loadspk(self):
        """Load neurons from multiple .spk files"""
        self.header = SPKHeader(self.path)
        for spkfname in self.header.spkfnames:
            path = os.path.join(self.path, spkfname)
            neuron = Neuron(path, sort=self)
            self.header.read(neuron)
            self.alln[neuron.id] = neuron # save it


class TrackSort(object):
    """A kind of fake sort that holds a concatenation of neurons from all of a track's
    recordings. These neurons are stored as TrackNeurons"""
    def __init__(self, track=None):
        self.level = 3 # level in the hierarchy, just below Track
        self.path = track.path
        self.id = track.id
        self.tr = track
        self.alln = {} # dict to store all Neurons
        self.nspikes = None
        self.datetime = None
        self.pttype = None
        self.chanpos = None

    def get_n(self):
        """Return dict of neurons that meet MINRATE"""
        n = {}
        MINRATE = get_ipython().user_ns['MINRATE']
        for neuron in self.alln.values():
            if neuron.meanrate >= MINRATE:
                n[neuron.id] = neuron
        return n

    n = property(get_n)

    def get_qn(self):
        """Return dict of quiet neurons, ie those that fail to meet MINRATE"""
        qn = {}
        MINRATE = get_ipython().user_ns['MINRATE']
        for neuron in self.alln.values():
            if neuron.meanrate < MINRATE:
                qn[neuron.id] = neuron
        return qn

    qn = property(get_qn)

    nneurons = property(lambda self: len(self.n))
    nqneurons = property(lambda self: len(self.qn))
    nallneurons = property(lambda self: len(self.alln))

    def load(self):
        """Load TrackNeurons by concatenating spikes from neurons from all recordings"""
        tr = self.tr
        rids = sorted(tr.r.keys()) # all recording ids in tr
        recs = [ tr.r[rid] for rid in rids ]
        # copy some attribs from first sort, should be the same for all of them:
        sort = recs[0].sort
        datetime0 = sort.datetime # start of acquisition (t=0) of first recording
        self.datetime = datetime0
        self.pttype = sort.pttype
        self.chanpos = sort.chanpos
        # get the union of all nids in recs:
        nids = tr.get_allnids()
        spikes = {}
        for nid in nids:
            spikes[nid] = [] # init each value to empty list
        alln = {} # dict of first Neurons encountered across recordings
        for rec in recs:
            # store time delta between start of track and start of rec:
            rec.td = td2usec(rec.sort.datetime - datetime0) # (us)
            rec.tdsec = rec.td / 1e6
            rec.tdmin = rec.tdsec / 60
            rec.tdhour = rec.tdmin / 60
            # for each neuron in this recording append appropriately offset spikes
            # array to entry in spikes dict:
            for n in rec.alln.values():
                spikes[n.id].append(n.spikes + rec.td)
                # for each nid, store the first neuron encountered when iterating over
                # recordings;
                if n.id not in alln:
                    alln[n.id] = n

        nspikes = 0 # add them up
        for nid in nids:
            spikes[nid] = np.hstack(spikes[nid]) # concatenate each nid's spikes arrays:
            assert (np.sort(spikes[nid]) == spikes[nid]).all() # should come out sorted
            # replace Neuron with TrackNeuron:
            n = alln[nid]
            tn = TrackNeuron(self)
            # point TrackNeuron attribs to relevant Neuron attribs, probably not copies:
            tn.id = nid
            tn.descr = n.descr
            tn.pos = n.pos
            # ptcs version 3 added sigma:
            try: tn.sigma = n.sigma
            except AttributeError: pass
            tn.nchans = n.nchans
            tn.chans = n.chans
            tn.maxchan = n.maxchan
            tn.wavedata = n.wavedata
            tn.wavestd = n.wavestd
            # assign spikes and calc static attribs:
            tn.spikes = spikes[nid]
            tn.nspikes = len(tn.spikes)
            tn.trange = tn.spikes[0], tn.spikes[-1]
            tn.dt = tn.trange[1] - tn.trange[0]
            tn.dtsec = tn.dt / 1e6
            tn.dtmin = tn.dtsec / 60
            tn.dthour = tn.dtmin / 60
            
            alln[nid] = tn # replace
            nspikes += tn.nspikes

        self.nspikes = nspikes
        self.alln = alln # save it
