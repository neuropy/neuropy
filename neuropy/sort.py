"""Defines the Sort class"""

from __future__ import division

import os
import StringIO
import datetime

import core
from core import dictattr, rstrip, eof, TAB, PTCSHeader, SPKHeader, EPOCH
from neuron import Neuron


class Sort(object):
    """A sort is a single spike extraction. Generally, sorts of the same name within
    the same track were extracted in the same spike sorting session"""
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
        """Return dict of (quiet) neurons, ie those that fail to meet MINRATE"""
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
