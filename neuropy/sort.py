"""Defines the Sort class"""

import os
import StringIO

from core import dictattr, rstrip, TAB
from neuron import Neuron, ConstrainedNeuron


class Sort(object):
    """A sort is a single spike extraction. Generally, Sorts of the same name within
    the same track were extracted in the same spike sorting session"""
    def __init__(self, path, id=None, recording=None):
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.path = path
        self.id = id
        self.r = recording
        self.n = dictattr() # store Neurons in a dictionary with attrib access
        self.cn = dictattr() # store ConstrainedNeurons in a dictionary with attrib acces

    def get_name(self):
        dirname = os.path.split(self.path)[-1]
        return rstrip(dirname, '.sort')

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
        fnames = [ fname for fname in os.listdir(self.path)
                   if os.path.isfile(os.path.join(self.path, fname))
                   and fname.endswith('.spk') ] # spike filenames
        '''
        # Look for neuron2pos.py file, which contains a dict mapping from neuron id to (x, y) position
        if 'neuron2pos.py' in os.listdir(self.path):
            os.chdir(self.path)
            from neuron2pos import neuron2pos
        '''
        for fname in fnames:
            path = os.path.join(self.path, fname)
            neuron = Neuron(path, sort=self)
            neuron.load() # load the neuron
            self.n[neuron.id] = neuron # save it
            #self.__setattr__('n' + str(neuron.id), neuron) # add shortcut attrib
            # repeat for ConstrainedNeurons
            cneuron = ConstrainedNeuron(path, sort=self)
            cneuron.load()
            self.cn[cneuron.id] = cneuron
            #self.__setattr__('cn' + str(cneuron.id), cneuron) # add shortcut attrib
            '''
            try: # binding neuron (x, y) positions
                neuron.pos = neuron2pos[neuron.id]
                cneuron.pos = neuron2pos[cneuron.id]
            except NameError: # there was no neuron2pos file
                pass
            '''
