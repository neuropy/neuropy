"""Defines the System class"""

#print 'importing System'

from Core import *
from Core import _model # ensure it's imported, in spite of leading _

class System(object):
    """A model System can have multiple modelling Runs"""
    def __init__(self, name=DEFAULTSYSTEMNAME, parent=_model):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.m = parent # save parent Model object
        self.name = name
        self.path = self.m.path + self.name + SLASH
        self.m.s[self.name] = self # add this System to its parent's dict of Systems, in case this System wasn't loaded by its parent
        self.r = {} # store Runs in a dictionary
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.m.writetree(string)
    # doesn't need a id2name or name2id method, since there are no system ids
    def load(self):
        if not os.path.isdir(self.path):
            raise NameError, 'Cannot find System(%r), path %r does not exist' % (self.name, self.path)
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        runNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname[0:2].isdigit() and dirname.count(' - ') == 1 ] # 1st 2 chars in dirname must be digits, must contain exactly 1 occurrence of ' - '
        for runName in runNames:
            run = Run(id=None, name=runName, parent=self) # make an instance using just the runName (let it figure out the run id)
            run.load() # load the Run
            self.r[run.id] = run # save it
        #if len(self.r) == 1:
        #   self.r = self.r.values[0] # pull it out of the dictionary
