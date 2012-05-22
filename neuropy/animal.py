"""Defines the Animal class"""

from __future__ import division

import os
import StringIO

from core import dictattr, TAB
from track import Track


class Animal(object):
    """An animal can have multiple tracks.
    Animals are identified by their globally unique name (e.g. ptc15)"""
    def __init__(self, path):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # string buffer to print tree hierarchy to
        self.path = path
        self.tr = dictattr() # store tracks in a dictionary with attrib access

    def get_name(self):
        return os.path.split(self.path)[-1]
    
    name = property(get_name)

    def get_id(self):
        return self.name

    id = property(get_id)

    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self, string):
        """Write to self's tree buffer"""
        self.treebuf.write(string)

    def load(self):
        treestr = self.level*TAB + self.id + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)
        # all track folder names for this animal:
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname.lower().startswith('tr') ]
        dirnames.sort() # alphabetical order
        for dirname in dirnames:
            path = os.path.join(self.path, dirname)
            track = Track(path, animal=self)
            track.load()
            self.tr[track.id] = track
            self.__setattr__('tr' + str(track.id), track) # add shortcut attrib
