"""Defines the Track class"""

#print 'importing Track'

from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class Track(object):
    """A Track can have multiple Recordings"""
    def __init__(self, id=DEFAULTTRACKID, name=None, parent=None):

        from Cat import Cat

        self.level = 2 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        if parent == None:
            try:
                self.c = _data.c[DEFAULTCATID] # see if the default Cat has already been init'd
            except KeyError:
                self.c = Cat() # init the default Cat...
                _data.c[self.c.id] = self.c  # ...and add it to the default Data object's list of Cats
        else:
            self.c = parent # save parent Cat object
        if id is not None:
            name = self.id2name(self.c.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'Track id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.c.path + self.name + SLASH
        self.c.t[self.id] = self # add/overwrite this Track to its parent's dict of Tracks, in case this Track wasn't loaded by its parent
        self.r = {} # store Recordings in a dictionary
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.c.writetree(string)
    def id2name(self, path, id):
        name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith('Track '+str(id)) ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Track id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        id = name.replace('Track ','',1) # replace first occurrence of 'Track ' with nothing, keep the rest
        if not id:
            raise NameError, 'Badly formatted Track name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self):

        from Recording import Recording

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        recordingNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname[0:2].isdigit() and dirname.count(' - ') == 1 ] # 1st 2 chars in dirname must be digits, must contain exactly 1 occurrence of ' - '
        for recordingName in recordingNames:
            recording = Recording(id=None, name=recordingName, parent=self) # make an instance using just the recording name (let it figure out the recording id)
            recording.load() # load the Recording
            self.r[recording.id] = recording # save it
        #if len(self.r) == 1:
        #   self.r = self.r.values[0] # pull it out of the dictionary
