"""Defines the Animal class"""

#print 'importing Animal'

from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class Animal(object):
    """An Animal can have multiple Tracks.
    Animals are identified by their unique id (e.g. ptc15)"""
    def __init__(self, id=None, parent=_data):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.d = parent # save the parent Data object
        self.id = id
        self.path = os.path.join(self.d.path, self.id)
        self.d.a[self.id] = self # add/overwrite this Animal to its parent's dict of Animals, in case this Animal wasn't loaded by its parent
        self.t = dictattr() # store Tracks in a dictionary with attrib access
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.d.writetree(string)
    def load(self):

        from Track import Track

        treestr = self.level*TAB + self.id + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname.lower().startswith('tr') ] # collect all track folder names for this animal
        for dirname in dirnames:
            track = Track(id=None, name=dirname, parent=self) # make an instance using just the track name (let it figure out the track id)
            track.load() # load the Track
            self.t[track.id] = track # save it
