"""Defines the Sort class"""

from __future__ import division

import os
import StringIO
import datetime

from core import dictattr, rstrip, eof, TAB, PTCSHeader, EPOCH
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
        self.n = dictattr() # store Neurons in a dictionary with attrib access
        #self.cn = dictattr() # store ConstrainedNeurons in a dictionary with attrib acces

    def get_name(self):
        return os.path.split(self.path)[-1]

    name = property(get_name)

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
                self.n[neuron.id] = neuron # save it
            assert eof(f), 'File %s has unexpected length' % self.path

    def get_datetime(self):
        """Return datetime object, calculated from header.datetime days since EPOCH"""
        return EPOCH + datetime.timedelta(days=self.header.datetime)

    datetime = property(get_datetime)

    def loadspk(self):
        """Load neurons from multiple .spk files"""
        fnames = [ fname for fname in os.listdir(self.path)
                   if os.path.isfile(os.path.join(self.path, fname))
                   and fname.endswith('.spk') ] # spike filenames
        
        # Look for neuron2pos.py file, which contains a dict mapping from neuron id to (x, y) position
        if 'neuron2pos.py' in os.listdir(self.path):
            os.chdir(self.path)
            from neuron2pos import neuron2pos
        
        for fname in fnames:
            path = os.path.join(self.path, fname)
            neuron = Neuron(path, sort=self)
            neuron.loadspk() # load the neuron
            self.n[neuron.id] = neuron # save it
            #self.__setattr__('n' + str(neuron.id), neuron) # add shortcut attrib
            # repeat for ConstrainedNeurons
            #cneuron = ConstrainedNeuron(path, sort=self)
            #cneuron.load()
            #self.cn[cneuron.id] = cneuron
            #self.__setattr__('cn' + str(cneuron.id), cneuron) # add shortcut attrib
            
            try: # binding neuron (x, y) positions
                neuron.record.xpos, neuron.record.ypos = neuron2pos[neuron.id]
                #cneuron.pos = neuron2pos[cneuron.id]
            except NameError: # there was no neuron2pos file
                pass
