"""Defines the Cat class"""

#print 'importing Cat'

from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class Cat(object):
    """A Cat can have multiple Tracks"""
    def __init__(self, id=DEFAULTCATID, name=None, parent=_data):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.d = parent # save the parent Data object
        if id is not None:
            name = self.id2name(self.d.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'Cat id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.d.path + self.name + SLASH
        self.d.c[self.id] = self # add/overwrite this Cat to its parent's dict of Cats, in case this Cat wasn't loaded by its parent
        self.t = {} # store Tracks in a dictionary
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.d.writetree(string)
    def id2name(self, path, id):
        if len(str(id)) == 1: # if id is only 1 digit long
            id = '0'+str(id) # add a leading zero
        name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith('Cat '+str(id)) ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Cat id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        id = name.replace('Cat ','',1) # replace first occurrence of 'Cat ' with nothing, keep the rest
        if not id:
            raise NameError, 'Badly formatted Cat name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self):

        from Track import Track

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        trackNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Track ') ]
        for trackName in trackNames:
            track = Track(id=None, name=trackName, parent=self) # make an instance using just the track name (let it figure out the track id)
            track.load() # load the Track
            self.t[track.id] = track # save it
        #if len(self.t) == 1:
        #   self.t = self.t.values[0] # pull it out of the dictionary
