"""Defines the Rip class"""

#print 'importing Rip'

from Core import *

class Rip(object):
    """A Rip is a single spike extraction. Generally, Rips of the same name within the same Track
    were generated with the same spike template, though of course Rips in different Tracks must
    be generated from different templates, even if the Rips have the same name. In the context of a
    Model, a Rip is a set of spike times generated with a certain set of modelling parameters."""

    from Recording import Recording

    def __init__(self, id=None, name=None, parent=Recording):
        self.level = 4 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.r = parent() # init parent Recording object
        except TypeError: # parent is an instance, not a class
            self.r = parent # save parent Recording object
        if name is None:
            raise ValueError, 'rip name can\'t be None'
        # rips don't have ids, at least for now. Just names
        self.id = id # not really used by the Rip class, just there for user's info
        self.name = name
        self.path = self.r.path + self.name + '.rip' + SLASH # have to add .rip extension to rip name to get its actual folder name
        self.n = {} # store Neurons in a dictionary
    '''
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)
    '''
    def load(self):

        from Neuron import Neuron

        #treestr = self.level*TAB + self.name + '/'
        #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        neuronNames = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.spk') ] # returns spike filenames without their .spk extension
        for neuronName in neuronNames:
            neuron = Neuron(id=None, name=neuronName, parent=self) # make an instance using just the neuron name (let it figure out the neuron id)
            neuron.load() # load the neuron
            self.n[neuron.id] = neuron # save it
        # then, maybe add something that loads info about the rip, say from some file describing the template used, and all the thresholds, exported to the same folder by SURF
        # maybe also load the template used for the rip, perhaps also stored in the same folder

