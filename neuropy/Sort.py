"""Defines the Sort class"""

from Core import *

class Sort(object):
    """A Sort is a single spike extraction. Generally, Sorts of the same name within the same Track
    were generated with the same spike template, though of course Sorts in different Tracks must
    be generated from different templates, even if the Sorts have the same name"""

    from Recording import Recording

    def __init__(self, id=None, name=None, parent=Recording):
        self.level = 4 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.r = parent() # init parent Recording object
        except TypeError: # parent is an instance, not a class
            self.r = parent # save parent Recording object
        if name is None:
            raise ValueError, 'sort name can\'t be None'
        # sorts don't have ids, at least for now. Just names
        self.id = id # not really used by the Sort class, just there for user's info
        self.name = name
        self.path = os.path.join(self.r.path, self.name) + '.sort' # have to add .sort extension to sort name to get its actual folder name
        self.n = dictattr() # store Neurons in a dictionary with attrib access
        self.cn = dictattr() # store ConstrainedNeurons in a dictionary with attrib acces
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

        from Neuron import Neuron, ConstrainedNeuron

        #treestr = self.level*TAB + self.name + '/'
        #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        neuronNames = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(self.path)
                        if os.path.isfile(os.path.join(self.path, fname))
                        and fname.endswith('.spk') ] # returns spike filenames without their .spk extension

        # Look for neuron2pos.py file, which contains a dict mapping from neuron id to (x, y) position
        if 'neuron2pos.py' in os.listdir(self.path):
            os.chdir(self.path)
            from neuron2pos import neuron2pos

        for neuronName in neuronNames:
            neuron = Neuron(id=None, name=neuronName, parent=self) # make an instance using just the neuron name (let it figure out the neuron id)
            neuron.load() # load the neuron
            self.n[neuron.id] = neuron # save it
            self.__setattr__('n' + str(neuron.id), neuron) # add shortcut attrib
            # repeat for ConstrainedNeurons
            cneuron = ConstrainedNeuron(id=None, name=neuronName, parent=self)
            cneuron.load()
            self.cn[cneuron.id] = cneuron
            self.__setattr__('cn' + str(cneuron.id), cneuron) # add shortcut attrib
            try: # binding neuron (x, y) positions
                neuron.pos = neuron2pos[neuron.id]
                cneuron.pos = neuron2pos[cneuron.id]
            except NameError: # there was no neuron2pos file
                pass

        # then, maybe add something that loads info about the sort, say from some file describing the template used, and all the thresholds, exported to the same folder by SURF
        # maybe also load the template file used for the sort, perhaps also stored in the same folder
