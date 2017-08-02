"""Defines the Animal class"""

import os
from io import StringIO

from core import dictattr, tolist, TAB
from track import Track


class Animal(object):
    """An animal can have multiple tracks.
    Animals are identified by their globally unique name (e.g. ptc15)"""
    def __init__(self, path):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO() # string buffer to print tree hierarchy to
        self.path = path
        self.tr = dictattr() # store tracks in a dictionary with attrib access

    def get_longname(self):
        return os.path.split(self.path)[-1]
    
    longname = property(get_longname)

    def get_name(self):
        """Return shorter names for blab mice"""
        long2shortstrain = {'BL6': 'bl',
                            'PVCre': 'pvc',
                            'Ntsr1-Cre': 'nts'}
        fields = self.longname.split('_') # split name into fields separated by underscores
        if len(fields) == 1: # no underscores, something like 'ptc22'
            return fields[0]
        if len(fields) == 2: # 1 underscore, something like 'BL6_0348'
            strain, num = fields
            year = ''
        elif len(fields) == 3: # 2 underscores, something like 'BL6_2017_0001'
            strain, year, num = fields
        else:
            raise ValueError("Don't know how to parse long name %r" % self.longname)
        strain = long2shortstrain[strain] # filter strain, strict capitalization
        if year:
            assert year.isnumeric()
            year = str(int(year) - 2000) # year wrt 2000
            # ensure year is exactly 2 digits:
            year= year.zfill(2)
            assert len(year) == 2
        assert num.isnumeric()
        num = str(int(num)) # drop any leading 0s
        num = num.zfill(2) # ensure animal number is at least 2 digits
        return ''.join((strain, year, num))

    name = property(get_name)

    def get_id(self):
        return self.name

    id = property(get_id)

    def tree(self):
        """Print tree hierarchy"""
        print(self.treebuf.getvalue(), end='')

    def writetree(self, string):
        """Write to self's tree buffer"""
        self.treebuf.write(string)

    def load(self, tracknames=None):
        treestr = self.level*TAB + self.id + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)
        if tracknames != None:
            tracknames = tolist(tracknames)
            dirnames = tracknames
        else:
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

    def get_rnames(self):
        rnames = []
        trids = sorted(self.tr)
        for trid in trids:
            rnames.append(self.tr[trid].rnames)
        return rnames

    rnames = property(get_rnames)
