"""Defines the Animal, Cat, and Rat classes"""

#print 'importing Animal'

from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class Animal(object):
    """An abstract Animal class, not meant to be instantiated directly.
    An Animal can have multiple Tracks. Animals are identified by their full name
    for clarity, since ids could be confusing (would id==5 mean Cat 5 or Rat 5?)"""
    def __init__(self, name=None, parent=None):
        self.level = 1 # level in the hierarchy
        self.treebuf = cStringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.d = parent # save the parent Data object
        self.name = name
        self.path = os.path.join(self.d.path, self.name)
        self.d.a[self.name] = self # add/overwrite this Animal to its parent's dict of Animals, in case this Animal wasn't loaded by its parent
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

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname.startswith('Track ') ]
        for dirname in dirnames:
            track = Track(id=None, name=dirname, parent=self) # make an instance using just the track name (let it figure out the track id)
            track.load() # load the Track
            self.t[track.id] = track # save it


class Cat(Animal):
    """This Animal is a Cat"""
    def __init__(self, id=CATID, parent=_data):
        #id = pad0s(id, ndigits=2) # returns a string
        name = 'Cat ' + str(id)
        self.id = int(id) # save it as an int (this will again remove any leading zeros)
        super(Cat, self).__init__(name=name, parent=parent)


class Rat(Animal):
    """This Animal is a Rat"""
    def __init__(self, id=RATID, parent=_data):
        #id = pad0s(id, ndigits=2) # returns a string
        name = 'Rat ' + str(id)
        self.id = int(id) # save it as an int (this will again remove any leading zeros)
        super(Rat, self).__init__(name=name, parent=parent)
