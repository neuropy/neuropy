"""Defines the Track class"""

#print 'importing Track'

from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class Track(object):
    """A Track can have multiple Recordings"""
    def __init__(self, id=TRACKID, name=None, parent=None):

        from Animal import Cat, Rat

        self.level = 2 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        if parent == None:
            try:
                self.a = _data.a[ANIMALNAME] # see if the default Animal has already been init'd
            except KeyError: # it hasn't, init the default animal, save it in self.a, and add it to _data.a
                if SPECIES == 'Cat':
                    defaultanimal = Cat(id=CATID, parent=_data)
                elif SPECIES == 'Rat':
                    defaultanimal = Rat(id=RATID, parent=_data)
                self.a = defaultanimal
                _data.a[defaultanimal.name] = defaultanimal
        else:
            self.a = parent # save parent Animal object
        if id is not None:
            name = self.id2name(self.a.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'Track id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = os.path.join(self.a.path, self.name)
        self.a.t[self.id] = self # add/overwrite this Track to its parent Animal's dict of Tracks, in case this Track wasn't loaded by its parent
        self.r = dictattr() # store Recordings in a dictionary with attrib access
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.a.writetree(string)
    def id2name(self, path, id):
        name = [ dirname for dirname in os.listdir(path)
                 if os.path.isdir(os.path.join(path, dirname))
                 and dirname.startswith('Track %s' % id) ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Track id: %s' % id
        else:
            return name[0] # pull the string out of the list
    def name2id(self, name):
        id = name.replace('Track ', '', 1) # replace first occurrence of 'Track ' with nothing, keep the rest
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
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname[0].isdigit() ] # 1st char in dirname must be a digit, that's all
        for dirname in dirnames:
            recording = Recording(id=None, name=dirname, parent=self) # make an instance using just the recording name (let it figure out the recording id)
            recording.load() # load the Recording
            self.r[recording.id] = recording # save it
