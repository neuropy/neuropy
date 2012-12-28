"""Defines the Track class"""

from __future__ import division

import os
import StringIO

import numpy as np

import pylab as pl

import core
from core import dictattr, TAB, td2usec, pad0s
from recording import Recording


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
        # collect recording names: 1st char in dirname must be a digit, that's all:
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname[0].isdigit() ]
        dirnames.sort() # alphabetical order
        for dirname in dirnames:
            path = os.path.join(self.path, dirname)
            recording = Recording(path, track=self)
            recording.load()
            self.r[recording.id] = recording
            self.__setattr__('r' + str(recording.id), recording) # add shortcut attrib
        self.rnames = dirnames # easy way to print out all recording names

        rids = np.sort(self.r.keys()) # all recording ids in self
        if len(rids) > 0:
            # calculate total track duration, this is slightly different from what you get
            # from the source .srf files, because exact recording duration is not exported,
            # only experiment din and spike times, which are used to generate rec.trange
            r0 = self.r[rids[0]]
            r1 = self.r[rids[-1]]
            if hasattr(r0, 'datetime') and hasattr(r1, 'datetime'):
                self.dt = td2usec(r1.datetime - r0.datetime) - r0.trange[0] + r1.trange[1]
                self.dtsec = self.dt / 1e6
                self.dtmin = self.dtsec / 60
                self.dthour = self.dtmin / 60

    def get_nids(self, rids=None):
        """Return nids of normal (active) neurons common to all recordings specified in
        rids. Active neurons in a recording are those with at least MINRATE mean
        spike rate during the recording"""
        if rids == None:
            rids = self.r.keys() # all recording ids in self
            return np.unique(np.hstack([ self.r[rid].n.keys() for rid in rids ]))
        else:
            # allow int rids:
            rids = [ pad0s(rid, ndigits=2) for rid in rids ]
            nids = [ self.r[rid].n.keys() for rid in rids ]
            return core.intersect1d(nids, assume_unique=True)

    def get_allnids(self, rids=None):
        """Return nids of all neurons (active and quiet) common to all recordings
        specified in rids"""
        if rids == None:
            rids = self.r.keys() # all recording ids in self
            return np.unique(np.hstack([ self.r[rid].alln.keys() for rid in rids ]))
        else:
            allnids = [ self.r[rid].alln.keys() for rid in rids ]
            return core.intersect1d(allnids, assume_unique=True)

    def get_nspikes(self, rids=None):
        """Return total number of spikes in recordings specified by rids"""
        if rids == None:
            rids = self.r.keys() # all recording ids in self
        nspikes = [ self.r[rid].nspikes for rid in rids ]
        return sum(nspikes)

    def get_meanrates(self):
        """Return mean firing rates of all neurons across all recordings.
        Neurons are counted as many times as they have spikes in a given recording.
        Maybe this should be weighted by the duration of each recording"""
        meanrates = []
        for r in self.r.values():
            meanrates.append([n.meanrate for n in r.alln.values()])
        return np.hstack(meanrates)

    allnids = property(get_allnids)
    nallneurons = property(lambda self: len(self.allnids))
    nspikes = property(get_nspikes)
    meanrates = property(get_meanrates)

    def meanratehist(self, bins=None):
        f = pl.figure()
        if bins == None:
            bins = np.arange(0, 1, 0.01)
        pl.hist(self.meanrates, bins=bins)
