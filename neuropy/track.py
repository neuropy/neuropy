"""Defines the Track class"""

from __future__ import division

import os
import StringIO

from core import dictattr, TAB
from recording import Recording


class Track(object):
    """A track can have multiple recordings"""
    def __init__(self, path, animal=None):
        self.level = 2 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # string buffer to print tree hierarchy to
        self.path = path
        self.animal = animal
        if animal != None:
            # update parent animal's track dict, in case self wasn't loaded by its parent
            animal.tr[self.id] = self
        self.r = dictattr() # store recordings in a dictionary with attrib access

    def get_name(self):
        return os.path.split(self.path)[-1]
    
    name = property(get_name)

    def get_id(self):
        # replace first occurrence of 'tr' with nothing, keep the rest:
        id = self.name.lower().replace('tr', '', 1)
        if not id:
            raise NameError('Badly formatted track name: %s' % self.name)
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, as in '7c', leave it as a string
        return id
        
    id = property(get_id)

    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        if self.animal != None:
            self.animal.writetree(string)

    def load(self):
        treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname[0].isdigit() ] # collect recording names: 1st char in dirname must be a digit, that's all
        for dirname in dirnames:
            path = os.path.join(self.path, dirname)
            recording = Recording(path, track=self)
            recording.load()
            self.r[recording.id] = recording
            self.__setattr__('r' + str(recording.id), recording) # add shortcut attrib
