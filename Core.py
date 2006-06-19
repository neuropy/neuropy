r"""Core neuropy functions and classes

                         object hierarchies:

                                level:

                Data              0                Model
                  |                                  |
                 Cat              1               System
                  |                                  |
                Track             2                  |
                  |                                  |
              Recording           3                 Run
             /         \                           /   \
       Experiment      Rip        4       Experiment    Rip
            |           |                      |         |
          Movie       Neuron      5          Movie     Neuron
"""

DEFAULTDATAPATH = 'C:/data/' # the convention in neuropy is that all 'path' var names have a trailing slash
DEFAULTMODELPATH = 'C:/model/'
DEFAULTCATID    = 15
DEFAULTSYSTEMNAME = 'Cat 15'
DEFAULTTRACKID  = '7c'
RIPKEYWORDS = ['best'] # a Rip with one of these keywords (listed in decreasing priority) will be loaded as the default Rip for its Recording/Run
SLASH = '/' # use forward slashes instead of having to use double backslashes
TAB = '    ' # 4 spaces

DEFAULTMOVIEPATH = 'C:/pub/Movies/'
DEFAULTMOVIENAME = 'mseq32.m'

import os
import types
from pprint import pprint
import struct
import re
import StringIO
import sys

import numpy as np
import pylab as pl
import scipy.signal as sig
import Dimstim.Movies
from Dimstim.Core import buildSweepTable

# Rips should really have ids to make them easier to reference to: r[83].rip[0] instead of r[83].rip['conservative spikes'] - this means adding id prefixes to rip folder names (or maybe suffixes: 'conservative spikes.0.rip', 'liberal spikes.1.rip', etc...). Prefixes would be better cuz they'd force sorting by id in explorer (which uses alphabetical order) - ids should be 0-based of course
# worry about conversion of ids to strings: some may be only 1 digit and may have a leading zero!
# maybe make two load() f'ns for Experiment and Neuron: one from files, and a future one from a database
# make a save() f'n that pickles the object (including any of its results, like its STA, tuning curve points, etc)? - just use IPython's %store

'''
def str2(data):
    if type(data) is types.IntTypes:
        s = str(data)
        if len(s) == 1:
            s = '0'+s # add a leading zero for single digits
'''
def txtdin2binarydin(fin, fout):
    """Converts a csv text .din file to an int64 binary .din file"""
    fi = file(fin, 'r') # open the din file for reading in text mode
    fo = file(fout,'wb') # for writing in binary mode
    for line in fi:
        line = line.split(',')
        '''
        # for old NVS display, converts from NVS condition numbers (which increment with repeats) to Dimstim sweepis (which don't)
        nruns = 18
        line[1] = int(line[1]) % nruns
        '''
        fo.write( struct.pack('@qq',int(line[0]),int(line[1])) ) # read both values in as a C long longs, using the system's native ('@') byte order
    fi.close()
    fo.close()
    print 'Converted ascii din: ', fin, ' to binary din: ', fout

def renameSpikeFiles(path, newname):
    """Renames all .spk files in path to newname, retaining their '_t##.spk' ending"""
    for fname in os.listdir(path):
        if fname.endswith('.spk'):
            i=fname.find('_t')
            if i!=-1:
                newfname = newname+fname[i::]
                print newfname
                os.rename(path+SLASH+fname, path+SLASH+newfname)

def warn(msg,level=2,exit_val=1):
    """Standard warning printer. Gives formatting consistency. Stolen from IPython.genutils"""
    if level>0:
        header = ['','','WARNING: ','ERROR: ','FATAL ERROR: ']
        print >> sys.stderr, '%s%s' % (header[level],msg)
        if level == 4:
            print >> sys.stderr,'Exiting.\n'
            sys.exit(exit_val)
'''
def warn(msg):
    import warnings
    warnings.warn(msg, category=RuntimeWarning, stacklevel=2)
'''
def unique(inseq):
    """Return unique items from a 1-dimensional sequence. Stolen from numpy.unique(), modified to return list instead of array"""
    # Dictionary setting is quite fast.
    outseq = {}
    for item in inseq:
        outseq[item] = None
    return list(outseq.keys())
'''
def unique(objlist):
    """Returns the input list minus any repeated objects it may have had. Also defined in Dimstim"""
    return list(set(objlist)) # this requires Python >= 2.4
'''
'''
def unique(objlist):
    """Does in-place removal of non-unique objects in a list of objects"""
    for (i,obj1) in enumerate(objlist):
        for (j,obj2) in enumerate(objlist):
            if i != j and obj1 == obj2:
                del objlist[j]
'''
def iterable(y):
    """Check if the input is iterable, stolen from numpy.iterable()"""
    try: iter(y)
    except: return 0
    return 1
'''
def tolist(obj):
    """Takes either scalar or sequence input and returns a list,
    useful when you want to iterate over an object (like in a for loop),
    and you don't want to have to do type checking or handle exceptions
    when the object isn't a sequence"""
    try: # assume obj is a sequence
        return list(obj) # converts any sequence to a list
    except TypeError: # obj is probably a scalar
        return [obj] # converts any scalar to a list
'''
def approx(a, b, rtol=1.e-14, atol=1.e-14):
    """Returns a boolean array describing which components of a and b are equal
    subject to given tolerances. The relative error rtol must be positive and << 1.0
    The absolute error atol comes into play for those elements of y that are very
    small or zero; it says how small x must be also. Copied and modified from
    numpy.allclose()"""
    x = np.array(a, copy=False)
    y = np.array(b, copy=False)
    print x.shape
    print y.shape
    return np.less(np.absolute(x-y), atol + rtol * np.absolute(y))

def histogramSorted(sorteda, bins=10, range=None):
    """Builds a histogram, stolen from numpy.histogram(), modified to assume sorted input"""
    a = np.asarray(sorteda).ravel()
    if not iterable(bins):
        if range is None:
            range = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins, endpoint=False)
    #n = np.sort(a).searchsorted(bins)
    n = a.searchsorted(bins)
    n = np.concatenate([n, [len(a)]]) # this adds a bin that includes overflow points
    n = n[1:]-n[:-1] # subtracts a shifted version of itself
    #if normed:
    #   db = bins[1] - bins[0]
    #   return 1.0/(a.size*db) * n, bins # this seems a bit weird
    #else:
    return n, bins

def sah(t, y, ts, keep=False):
    """Resample using sample and hold. Returns resampled values at ts given the original points (t,y)
    such that the resampled values are just the most recent value in y (think of a staircase with non-uniform steps).
    Assumes that t is sorted. t and ts arrays should be of the same data type. Contributed by Robert Kern."""
    i = np.searchsorted(t, ts) - 1 # find where ts falls in t, dec so you get indices that point to the most recent value in y
    i = np.where(i < 0, 0, i) # handle the cases where ts is smaller than the first point.
    '''this has an issue of not keeping the original data point where ts == t'''
    if keep:
        # The following ensures that the original data point is kept when ts == t, doesn't really work if the shortest ISI is less than tres in ts
        di = np.diff(i).nonzero()[0] # find changes in i, nonzero() method returns a tuple, pick the result for the first dim with [0] index
        si = approx(t[1::], ts[di]) # check at those change indices if t ~= ts (ignoring potential floating point representational inaccuracies). If so, inc i at that point so you keep y at that point.
        #print i
        i[di[si]] += 1
        #print i
    return y[i]

def corr(x,y):
    """Returns correlation of signals x and y. This should be equivalent to np.corrcoef(),
    but that one doesn't seem to work for signals with zeros in them. Check how std() works exactly"""
    x = np.array(x)
    y = np.array(y)
    return ((x * y).mean() - x.mean() * y.mean()) / (x.std() * y.std())

class Data(object): # use 'new-style' classes
    """Data can have multiple Cats"""
    def __init__(self, dataPath=DEFAULTDATAPATH):
        self.level = 0 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.name = 'Data'
        self.path = dataPath
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        # Data has no parent to write to
    def load(self):
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.c = {} # store Cats in a dictionary
        catNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Cat ') ] # os.listdir() returns all dirs AND files
        for catName in catNames:
            cat = Cat(id=None, name=catName, parent=self) # make an instance using just the catName (let it figure out the cat id)
            cat.load() # load the Cat
            self.c[cat.id] = cat # save it
        #if len(self.c) == 1:
        #   self.c = self.c.values[0] # pull it out of the dictionary


class Model(Data):
    """Abstract model data class. Model can have multiple modelling Runs"""
    def __init__(self, modelPath=DEFAULTMODELPATH):
        self.level = 0 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.name = 'Model'
        self.path = modelPath
    def load(self):
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.s = {} # store model Systems in a dictionary
        systemNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) ] # os.listdir() returns all dirs AND files
        for systemName in systemNames:
            system = System(name=systemName, parent=self) # make an instance using the systemName
            system.load() # load the System
            self.s[system.name] = system # save it
        #if len(self.s) == 1:
        #   self.s = self.s.values[0] # pull it out of the dictionary


class Cat(object):
    """A Cat can have multiple Tracks"""
    def __init__(self, id=DEFAULTCATID, name=None, parent=Data):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.d = parent() # init parent Data object
        except TypeError: # parent is an instance, not a class
            self.d = parent # save parent Data object
        if id is not None:
            name = self.id2name(self.d.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'cat id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.d.path + self.name + SLASH
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
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
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.t = {} # store Tracks in a dictionary
        trackNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Track ') ]
        for trackName in trackNames:
            track = Track(id=None, name=trackName, parent=self) # make an instance using just the track name (let it figure out the track id)
            track.load() # load the Track
            self.t[track.id] = track # save it
        #if len(self.t) == 1:
        #   self.t = self.t.values[0] # pull it out of the dictionary


class System(object):
    """A model System can have multiple modelling Runs"""
    def __init__(self, name=DEFAULTSYSTEMNAME, parent=Model):
        self.level = 1 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.m = parent() # init parent Model object
        except TypeError: # parent is an instance, not a class
            self.m = parent # save parent Model object
        self.name = name
        self.path = self.m.path + self.name + SLASH
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.m.writetree(string)
    # doesn't need a id2name or name2id method, since there are no system ids
    def load(self):
        if not os.path.isdir(self.path):
            raise NameError, 'Cannot find System(%s), path %s does not exist' % (repr(self.name), repr(self.path))
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.r = {} # store Runs in a dictionary
        runNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname[0:2].isdigit() and dirname.count(' - ') == 1 ] # 1st 2 chars in dirname must be digits, must contain exactly 1 occurrence of ' - '
        for runName in runNames:
            run = Run(id=None, name=runName, parent=self) # make an instance using just the runName (let it figure out the run id)
            run.load() # load the Run
            self.r[run.id] = run # save it
        #if len(self.r) == 1:
        #   self.r = self.r.values[0] # pull it out of the dictionary


class Track(object):
    """A Track can have multiple Recordings"""
    def __init__(self, id=DEFAULTTRACKID, name=None, parent=Cat):
        self.level = 2 # level in the hierarchy
        try:
            self.c = parent() # init parent Cat object
        except TypeError: # parent is an instance, not a class
            self.c = parent # save parent Cat object
        if id is not None:
            name = self.id2name(self.c.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'track id and name can\'t both be None'
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.id = id
        self.name = name
        self.path = self.c.path + self.name + SLASH
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
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
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.r = {} # store Recordings in a dictionary
        recordingNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname[0:2].isdigit() and dirname.count(' - ') == 1 ] # 1st 2 chars in dirname must be digits, must contain exactly 1 occurrence of ' - '
        for recordingName in recordingNames:
            recording = Recording(id=None, name=recordingName, parent=self) # make an instance using just the recording name (let it figure out the recording id)
            recording.load() # load the Recording
            self.r[recording.id] = recording # save it
        #if len(self.r) == 1:
        #   self.r = self.r.values[0] # pull it out of the dictionary


class Recording(object):
    """A Recording corresponds to a single SURF file, ie everything recorded between when
    the user hits record and when the user hits stop and closes the SURF file, including any
    pauses in between Experiments within that Recording. A Recording can have multiple Experiments,
    and multiple spike extractions, called Rips"""
    def __init__(self, id=None, name=None, parent=Track):
        self.level = 3 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.t = parent() # init parent Track object
        except TypeError: # parent is an instance, not a class
            self.t = parent # save parent Track object
        if id is not None:
            name = self.id2name(self.t.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'recording id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.t.path + self.name + SLASH
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.t.writetree(string)
    def id2name(self, path, id):
        if len(str(id)) == 1: # if id is only 1 digit long
            id = '0'+str(id) # add a leading zero
        name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith(str(id)+' - ') ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Recording id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        try:
            id = name[0:name.index(' - ')] # everything before the first ' - ', index() raises ValueError if it can't be found
        except ValueError:
            raise ValueError, 'Badly formatted Recording name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self):
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.e = {} # store Experiments in a dictionary
        experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.din') ] # returns din filenames without their .din extension
        for (experimentid, experimentName) in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
            experiment = Experiment(id=experimentid, name=experimentName, parent=self) # pass both the id and the name
            experiment.load() # load the Experiment
            self.e[experiment.id] = experiment # save it
        #if len(self.e) == 1:
        #   self.e = self.e.values[0] # pull it out of the dictionary
        self.rip = {} # store Rips in a dictionary
        ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
        defaultRipNames = [ ripName for ripName in ripNames for ripkeyword in RIPKEYWORDS if ripName.count(ripkeyword) ]
        if len(defaultRipNames) < 1:
            warn('Couldn\'t find a default Rip for Recording(%s)' % self.id)
        if len(defaultRipNames) > 1: # This could just be a warning instead of an exception, but really, some folder renaming is in order
            raise RuntimeError, 'More than one Rip folder in Recording(%s) has a default keyword: %s' %(self.id, defaultRipNames)
        for (ripid, ripName) in enumerate(ripNames): # ripids will be according to alphabetical order of ripNames
            rip = Rip(id=ripid, name=ripName, parent=self) # pass both the id and the name
            rip.load() # load the Rip
            self.rip[rip.name] = rip # save it
            # make the Neurons from the default Rip (if it exists in the Recording path) available in the Recording, so you can access them via r.n[nid] instead of having to do r.rip[name].n[nid]. Make them just another pointer to the data in r.rip[ripName].n
            for ripkeyword in RIPKEYWORDS[::-1]: # reverse the keywords so first one gets processed last
                if rip.name.count(ripkeyword): # if the keyword is in the ripName
                    self.n = self.rip[rip.name].n # make it the default Rip
        #if len(self.rip) == 1:
        #   self.rip = self.rip.values[0] # pull it out of the dictionary


class Run(object):
    """A Run corresponds to a single modelling run. A Run can have multiple Experiments
    and modelling Rips (a set of spike times generated with a certain set of modelling parameters)."""
    def __init__(self, id=None, name=None, parent=System):
        self.level = 3 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.s = parent() # init parent System object
        except TypeError: # parent is an instance, not a class
            self.s = parent # save parent System object
        if id is not None:
            name = self.id2name(self.s.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'recording id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.s.path + self.name + SLASH
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.s.writetree(string)
    def id2name(self, path, id):
        if len(str(id)) == 1: # if id is only 1 digit long
            id = '0'+str(id) # add a leading zero
        name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith(str(id)+' - ') ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Run id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        try:
            id = name[0:name.index(' - ')] # everything before the first ' - ', index() raises ValueError if it can't be found
        except ValueError:
            raise ValueError, 'Badly formatted Run name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self):
        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.e = {} # store Experiments in a dictionary
        experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.din') ] # returns din filenames without their .din extension
        for (experimentid, experimentName) in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
            experiment = Experiment(id=experimentid, name=experimentName, parent=self) # pass both the id and the name
            experiment.load() # load the Experiment
            self.e[experiment.id] = experiment # save it
        #if len(self.e) == 1:
        #   self.e = self.e.values[0] # pull it out of the dictionary
        self.rip = {} # store Rips in a dictionary
        ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
        defaultRipNames = [ ripName for ripName in ripNames for ripkeyword in RIPKEYWORDS if ripName.count(ripkeyword) ]
        if len(defaultRipNames) < 1:
            warn('Couldn\'t find a default Rip for Run(%s)' % self.id)
        if len(defaultRipNames) > 1: # This could just be a warning instead of an exception, but really, some folder renaming is in order
            raise RuntimeError, 'More than one Rip folder in Run(%s) has a default keyword: %s' %(self.id, defaultRipNames)
        for (ripid, ripName) in enumerate(ripNames): # ripids will be according to alphabetical order of ripNames
            rip = Rip(id=ripid, name=ripName, parent=self) # pass both the id and the name
            rip.load() # load the Rip
            self.rip[rip.name] = rip # save it
            # make the Neurons from the default Rip (if it exists in the Run path) available in the Run, so you can access them via r.n[nid] instead of having to do r.rip[name].n[nid]. Make them just another pointer to the data in r.rip[ripName].n
            for ripkeyword in RIPKEYWORDS[::-1]: # reverse the keywords so first one gets processed last
                if rip.name.count(ripkeyword): # if the keyword is in the ripName
                    self.n = self.rip[rip.name].n # make it the default Rip
        #if len(self.rip) == 1:
        #   self.rip = self.rip.values[0] # pull it out of the dictionary


class Experiment(object):
    """An Experiment corresponds to a single contiguous VisionEgg stimulus session.
    It contains information about the stimulus during that session, including
    the DIN values, the text header, and any Movies that were involved"""
    def __init__(self, id=None, name=None, parent=Recording): # Experiment IDs are 1-based in the .din filenames, at least for now. They should be renamed to 0-based. Here, they're treated as 0-based
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.r = parent() # init parent Recording object
        except TypeError: # parent is an instance, not a class
            self.r = parent # save parent Recording object
        if name is None:
            raise ValueError, 'experiment name can\'t be None'
        self.id = id # not really used by the Experiment class, just there for user's info
        self.name = name
        self.path = self.r.path
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)
    # doesn't need a id2name or name2id method, neither can really be derived from the other in an easy way (although could use re), the id is just chronological (which is also alphabetical) order, at least for now
    def load(self):
        f = file(self.path + self.name + '.din', 'rb') # open the din file for reading in binary mode
        self.din = np.fromfile(f, dtype=np.int64).reshape(-1,2) # reshape to nrows x 2 columns
        f.close()
        try:
            f = file(self.path + self.name + '.textheader', 'r') # open the textheader file for reading
            self.textheader = f.read() # read it all in
            f.close()
        except IOError:
            if type(self.r) is Recording: # parent is a Recording, which normally have textheaders in their Experiments
                warn('Error reading: <%s.textheader>, text header not loaded' % self.name)
            elif type(self.r) is Run: # parent is a Run, which don't have textheaders in their Experiments, so don't print a warning
                pass
            else:
                raise RuntimeError, 'parent is invalid type: %s %s' % (type(self.r), Run)
            self.textheader = '' # set to empty

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

        if self.textheader: # if it isn't empty
            names1 = locals().copy() # namespace before execing the textheader
            exec(self.textheader) # execute the textheader as Python code, maybe add some checks here to prevent changes to filesystem from accidental code, unload the os or sys modules or something?
            names2 = locals().copy() # namespace after
            newnames = [ n2 for n2 in names2 if not names1.has_key(n2) and n2 != 'names1' ] # names that were added to the namespace, excluding the 'names1' name itself
            for newname in newnames:
                self.__setattr__(newname, eval(newname)) # for each variable that was defined in the textheader, bind it as an attribute of this Experiment

            try:
                self.stims = unique(self.playlist) # self.stims is a non-repeating list of object oriented stim objects (Movie is the only possible kind right now) in this Experiment
            except AttributeError: # this was a simple non object-oriented stim, has no playlist
                self.stims = []
            #if len(self.stims) == 1:
            #   self.stims = self.stims[0] # get rid of the list
            '''If you inited a stim object(s) (like a movie) while execing the textheader, you didn't have a chance to pass this exp as the parent in the init. So just set the attribute manually:'''
            for s in self.stims:
                s.e = self
                try: # this'll probably only apply to Movies stim, cuz others won't have fnames
                    if s.name == None:
                        s.name = s.fname # fname should've been defined when loading in the textheader
                except AttributeError: # probably not a Movie stim
                    pass
                # Search self.moviepath string (from textheader) for 'Movies' word (preferably case insensitive). Everything after that is the relative path to your base movies folder. Eg, if self.moviepath = 'C:\\Desktop\\Movies\\reliability\\e\\single\\', then set self.relpath = 'reliability\\e\\single\\'
                try:
                    s.moviepath = s.moviepath.replace('\\','/') # replace annoying double backslashes with single forward slashes, which seem to work
                    s.relpath = s.moviepath[ s.moviepath.index('Movies/')+len('Movies/') :: ]
                    s.path = path + s.relpath
                except AttributeError: # this Movie was manually inited, not loaded from a textheader. s.moviepath doesn't exist, use s.path instead. Or it might not even be a Movie
                    pass
                '''
                # Also, if you inited a stim that needs to be loaded (like a movie), maybe you should also load it now (this wasn't done when execing the textheader)
                try:
                    s.load()
                except:
                    pass
                '''
            # Generate the sweeptable here, no need to load if from files anymore...
            #self.sweeptable = {[]} # dictionary of lists, ie sweeptable={'ori',[0,45,90],'sfreq',[1,1,1]}
            # so you index into it with self.sweeptable['var'][sweepi]
            # vars = self.sweeptable.keys()

            if self.stims: # this Experiment has object-oriented stim(s)
                for s in self.stims:
                    varvals={} # init a dictionary that will contain variable values
                    for var in s.varlist:
                        varvals[var]=eval('s.'+var) # generate a dictionary with var:val entries, to pass to buildSweepTable
                    s.sweepTable = buildSweepTable(s.varlist, varvals, s.nruns, s.shuffleRuns, s.blankSweep, s.shuffleBlankSweeps, makeSweepTableText=0)[0] # passing varlist by reference, dim indices end up being modified
            else: # this is a simple stim (not object oriented)
                varvals={} # init a dictionary that will contain variable values
                for var in self.varlist:
                    varvals[var]=eval('self.'+var) # generate a dictionary with var:val entries, to pass to buildSweepTable
                self.sweepTable = buildSweepTable(self.varlist, varvals, self.nruns, self.shuffleRuns, self.blankSweep, self.shuffleBlankSweeps, makeSweepTableText=0)[0] # passing varlist by reference, dim indices end up being modified

            '''
            # Old code for creating a sweeptable file (used by Matlab and NVS):
            sweeptabletext = sweeptabletext.replace('[','') # get rid of brackets and ' in first line, these demarcate dimensions, but aren't needed in matlab
            sweeptabletext = sweeptabletext.replace(']','')
            sweeptabletext = sweeptabletext.replace('\'','')
            sweeptabletext = sweeptabletext.replace(',','') # also, replace any commas or spaces (in the first line) with tabs for delimiting
            sweeptabletext = sweeptabletext.replace(' ','\t')

            fname = string.replace(sys.argv[0],sys.path[0]+'\\','') # name of file that launched Python
            fname = fname.replace('.textheader','') # remove .textheader part of filename
            #fname = fname.replace('.py','') # remove .py part of filename
            fname += '.sweeptable' # add .sweeptable extension

            fullpathfname = sys.path[0]+'\\'+fname

            print 'Writing to file:', fullpathfname

            f = file(fullpathfname,'w')
            f.write(sweeptabletext)
            f.close()
            '''

            for defaultm in [MSEQ32, MSEQ16]: # check if this Experiment uses specific default movies
                for s in self.stims: # for all stims inited by the textheader
                    if s.name == defaultm.name and isinstance(s,Movie):
                        if defaultm.data == None: # see if this default movie has yet to be loaded
                            defaultm.load() # load this default movie
                        s.data = defaultm.data # point this movie's data to default movie data
                        treestr = s.level*TAB + s.name
                        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

        try:
            self.REFRESHTIME = int(round(1/float(self.REFRESHRATE)*1000000)) # in us, keep 'em integers
        except AttributeError:
            self.REFRESHTIME = self.din[1,0] - self.din[0,0] # use the time difference between the first two din instead
        #self.buildsweepranges()

    def buildsweepranges(self):
        self.sweepranges = []

    def code(self, neuron, **kwargs):
        """Returns a Neuron.Code object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
        try:
            return neuron.code(trange=trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].code(trange=trange, **kwargs) # neuron is probably a Neuron id

    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code objects"""
        code1 = self.code(neuron1, **kwargs)
        code2 = self.code(neuron2, **kwargs)
        return corr(code1.c, code2.c)

    class CodeCorrPDF(object):
        def __init__(self, experiment, range=(-1, 1), nbins=100, normed=True, **kwargs):
            self.e = experiment
            self.r = self.e.r
            self.range = range
            self.nbins = nbins
            self.normed = normed
            self.kwargs = kwargs
        def __eq__(self, other):
            selfd = self.__dict__.copy()
            otherd = other.__dict__.copy()
            # Delete their n and c attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
            [ d.__delitem__(key) for d in [selfd, otherd] for key in ['n', 'c'] if d.has_key(key) ]
            if type(self) == type(other) and selfd == otherd:
                return True
            else:
                return False
        def calc(self):
            neurons = self.r.n.keys()
            nneurons = len(neurons)
            # this is too slow! use np's builtin corrcoef somehow
            corrs = [ self.e.codecorr(neurons[ni1], neurons[ni2], **self.kwargs) for ni1 in range(0,nneurons) for ni2 in range(ni1,nneurons) if ni1 != ni2 ]
            self.n, self.c = np.histogram(corrs, bins=self.nbins, normed=self.normed)
        def plot(self):
            pl.figure()
            barwidth = (range[1] - range[0]) / float(self.nbins)
            #pl.hist(self.n, bins=self.r, normed=0, bottom=0, width=None, hold=False) # doesn't seem to work
            pl.bar(left=self.c, height=self.n, width=barwidth, bottom=0, color='k', yerr=None, xerr=None, ecolor='k', capsize=3)
            pl.xlim
            pl.title('neuron code corrleation pdf - experiment %d - %s' % (self.e.id, self.e.name))
            if self.normed:
                pl.ylabel('probability')
            else:
                pl.ylabel('count')
            pl.xlabel('correlation coefficient')

    def codecorrpdf(self, **kwargs):
        """Returns an existing CodeCorrPDF object, or creates a new one if necessary"""
        try:
            self.ccpdfs
        except AttributeError: # doesn't exist yet
            self.ccpdfs = [] # create a list that'll hold CodeCorrPDF objects
        co = self.CodeCorrPDF(experiment=self, **kwargs) # init a new one
        for ccpdf in self.ccpdfs:
            if co == ccpdf: # need to define special == method for class CodeCorrPDF()
                return ccpdf # returns the first object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self.ccpdfs
        co.calc() # no matching object was found, calculate it
        self.ccpdfs.append(co) # add it to the object list
        return co

    def rate(self, neuron, **kwargs):
        """Returns a Neuron.Rate object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
        try:
            return neuron.rate(trange=trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].rate(trange=trange, **kwargs) # neuron is probably a Neuron id

    def ratepdf(self):
        pass

class Rip(object):
    """A Rip is a single spike extraction. Generally, Rips of the same name within the same Track
    were generated with the same spike template, though of course Rips in different Tracks must
    be generated from different templates, even if the Rips have the same name. In the context of a
    Model, a Rip is a set of spike times generated with a certain set of modelling parameters."""
    def __init__(self, id=None, name=None, parent=Recording):
        self.level = 4 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.r = parent() # init parent Recording object
        except TypeError: # parent is an instance, not a class
            self.r = parent # save parent Recording object
        if name is None:
            raise ValueError, 'rip name can\'t be None'
        # rips don't have ids, at least for now. Just names
        self.id = id # not really used by the Rip class, just there for user's info
        self.name = name
        self.path = self.r.path + self.name + '.rip' + SLASH # have to add .rip extension to rip name to get its actual folder name
    '''
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)
    '''
    def load(self):
        #treestr = self.level*TAB + self.name + '/'
        #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        self.n = {} # store Neurons in a dictionary
        neuronNames = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.spk') ] # returns spike filenames without their .spk extension
        for neuronName in neuronNames:
            neuron = Neuron(id=None, name=neuronName, parent=self) # make an instance using just the neuron name (let it figure out the neuron id)
            neuron.load() # load the neuron
            self.n[neuron.id] = neuron # save it
        # then, maybe add something that loads info about the rip, say from some file describing the template used, and all the thresholds, exported to the same folder by SURF
        # maybe also load the template used for the rip, perhaps also stored in the same folder


class Neuron(object):
    """A Neuron object\'s spike data spans all the Experiments within a Recording.
    If different Recordings have Rips with the same name, you can assume that the
    same spike template was used for all of those Recordings, and that therefore
    the neuron ids are the same"""
    def __init__(self, id=None, name=None, parent=Rip): # neuron names don't include the '.spk' ending, although neuron filenames do
        self.level = 5 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.rip = parent() # init parent Rip object
        except TypeError: # parent is an instance, not a class
            self.rip = parent # save parent Rip object
        if id is not None:
            name = self.id2name(self.rip.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'neuron id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.rip.path
    '''
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.rip.writetree(string)
    '''
    def id2name(self, path, id):
        #if len(str(id)) == 1: # if id is only 1 digit long
        #    id = '0'+str(id) # add a leading zero
        name = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(path) if os.path.isfile(path+fname) and \
               ( fname.find('_t'+str(id)+'.spk')!=-1 or fname.find('_t0'+str(id)+'.spk')!=-1 or fname.find('_t00'+str(id)+'.spk')!=-1 ) ] # have to deal with leading zero ids, go up to 3 digit ids, should really use a re to do this properly...
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Neuron id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        try:
            id = name[name.rindex('_t')+2::] # everything from just after the last '_t' to the end of the neuron name, index() raises ValueError if it can't be found
        except ValueError:
            raise ValueError, 'Badly formatted Neuron name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self): # or loadspikes()?
        f = file(self.path + self.name + '.spk', 'rb') # open the spike file for reading in binary mode
        self.spikes = np.fromfile(f, dtype=np.int64) # read in all spike times in us
        f.close()
        self.nspikes = len(self.spikes)
        #self.results = {} # a dictionary to store results in
        #treestr = self.level*TAB + self.name + '/'
        #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

    def cut(self, *args):
        """Returns all of the Neuron's spike times where tstart <= spikes <= tend"""
        if len(args) == 0: # passed nothing
            tstart = self.spikes[0]
            tend = self.spikes[-1]
        elif len(args) == 1: # passed a tuple
            tstart = args[0][0]
            tend = args[0][1]
        elif len(args) == 2: # passed tstart and tend as separate args
            tstart = args[0]
            tend = args[1]
        else:
            raise RuntimeError, 'Too many arguments'
        if tstart == 0: # shorthand for "from first spike" - would be problematic if a spike existed at t=0
            tstart = self.spikes[0]
        if tend == -1: # shorthand for "to last spike" - would be problematic if a spike existed at t=-1
            tend = self.spikes[-1]
        '''
        # this is what we're trying to do:
        return self.spikes[ (self.spikes >= tstart) & (self.spikes <= tend) ]
        self.searchsorted(values) method does it faster. It returns an index where the value would fit in self. The index is such that self[index-1] < value <= self[index]. In this formula self[self.size]=inf and self[-1]= -inf
        Another possibility could be to use a masked array instead?
        '''
        lo, hi = self.spikes.searchsorted([tstart, tend]) # returns indices where tstart and tend would fit in spikes
        if tend  == self.spikes[min(hi, self.nspikes-1)]: # if tend matches a spike (protect from going out of index bounds when checking)
            hi += 1 # inc to include a spike if it happens to exactly equal tend. This gives us end inclusion
            hi = min(hi, self.nspikes) # limit hi to max slice index (==max value index + 1)
        return self.spikes[ lo : hi ] # slice it

    def cutrel(self, *args):
        """Cuts Neuron spike times and returns them relative to tstart"""
        if len(args) == 0: # passed nothing
            tstart = self.spikes[0]
            tend = self.spikes[-1]
        elif len(args) == 1: # passed a tuple
            tstart = args[0][0]
            tend = args[0][1]
        elif len(args) == 2: # passed tstart and tend as separate args
            tstart = args[0]
            tend = args[1]
        else:
            raise RuntimeError, 'Too many arguments'
        #if tstart == 0: # shorthand for "from first spike"
        #   tstart = self.spikes[0]
        if tend == -1: # shorthand for "to last spike"
            tend = self.spikes[-1]
        tstart = np.int64(round(tstart)) # let's keep all the returned spikes as integers, us is more than accurate enough
        return self.cut(tstart, tend) - tstart

    def code(self, kind='binary', **kwargs):
        """Returns an existing Code object, or creates a new one if necessary"""
        try:
            self.codes
        except AttributeError: # self.codes doesn't exist yet
            self.codes = [] # create a list that'll hold Code objects
        if kind == 'binary':
            co = Neuron.BinaryCode(neuron=self, **kwargs) # init a new BinaryCode object
        else:
            raise ValueError, 'Unknown kind: %s' % repr(self.kind)
        for code in self.codes:
            if co == code: # need to define special == method for class Code()
                return code # returns the first Code object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self.codes
        co.calc() # no matching Rate was found, calculate it
        self.codes.append(co) # add it to the Code object list
        return co

    class Code(object):
        """Abstract spike code class"""
        def __init__(self, neuron=None, trange=None):
            self.neuron = neuron
            if trange == None:
                trange = self.neuron.spikes[0], self.neuron.spikes[-1]
            self.trange = trange
        def __eq__(self, other):
            selfd = self.__dict__.copy()
            otherd = other.__dict__.copy()
            # Delete their c and t attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
            [ d.__delitem__(key) for d in [selfd, otherd] for key in ['c', 't'] if d.has_key(key) ]
            if type(self) == type(other) and selfd == otherd:
                return True
            else:
                return False
        def plot(self):
            # plot some kind of long grid of white and black elements?
            pl.figure()

    class BinaryCode(Code):
        """Quantizes the spike train into a binary signal with a given time resolution"""
        def __init__(self, neuron=None, trange=None, tres=20000):
            super(Neuron.BinaryCode, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'binary'
            self.tres = tres
        def calc(self):
            # make the start of the timepoints be an even multiple of self.tres. Round down to the nearest multiple. Do the same for the end of the timepoints, but round up. This way, timepoints will line up for different code objects
            tstart = self.trange[0] - (self.trange[0] % self.tres)
            tend   = self.trange[1] - (self.trange[1] % self.tres)
            self.t = np.arange( tstart, tend+self.tres, self.tres ) # t sequence demarcates left bin edges, add tres to trange[1] to make t end inclusive
            s = self.neuron.cut(self.trange) # spike times
            self.c = np.zeros(len(self.t)) # binary code signal
            self.c[np.unique(self.t.searchsorted(s)) - 1] = 1 # dec index by 1 so that you get indices that point to the most recent bin edge. For each bin that has 1 or more spikes in it, set its bit to 1
        def plot(self):
            super(Neuron.BinaryCode, self).plot()
            pl.title('neuron %d - binary spike code' % self.neuron.id)

    def rate(self, kind='nisi', **kwargs):
        """Returns an existing Rate object, or creates a new one if necessary"""
        try:
            self.rates
        except AttributeError: # self.rates doesn't exist yet
            self.rates = [] # create a list that'll hold Rate objects
        if kind == 'bin':
            ro = Neuron.BinRate(neuron=self, **kwargs) # init a new BinRate object
        elif kind == 'nisi':
            ro = Neuron.nISIRate(neuron=self, **kwargs) # init a new nISIRate object
        elif kind == 'wnisi':
            ro = Neuron.WnISIRate(neuron=self, **kwargs) # init a new WnISIRate object
        elif kind == 'idp':
            ro = Neuron.IDPRate(neuron=self, **kwargs) # init a new IDPRate object
        elif kind == 'gauss':
            ro = Neuron.GaussRate(neuron=self, **kwargs) # init a new GaussRate object
        elif kind == 'rect':
            ro = Neuron.RectRate(neuron=self, **kwargs) # init a new RectRate object
        else:
            raise ValueError, 'Unknown kind: %s' % repr(self.kind)
        for rate in self.rates:
            if ro == rate: # need to define special == method for class Rate()
                return rate # returns the first Rate object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self.rates
        ro.calc() # no matching Rate was found, calculate it
        self.rates.append(ro) # add it to the Rate object list
        return ro

    class Rate(object):
        """Abstract firing rate class"""
        def __init__(self, neuron=None, trange=None):
            self.neuron = neuron
            if trange == None:
                trange = self.neuron.spikes[0], self.neuron.spikes[-1]
            self.trange = trange
        def __eq__(self, other):
            selfd = self.__dict__.copy()
            otherd = other.__dict__.copy()
            # Delete their r and t attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
            [ d.__delitem__(key) for d in [selfd, otherd] for key in ['r', 't', 'rawr', 'rawt'] if d.has_key(key) ]
            if type(self) == type(other) and selfd == otherd:
                return True
            else:
                return False
        def plot(self):
            pl.figure()
            pl.plot(self.t, self.r)
            # diagnostic for comparing interpolated to raw nisi rate:
            #pl.plot(self.rawt, self.rawr, 'r+')
            pl.xlabel('t')
            pl.ylabel('spike rate')

    class BinRate(Rate):
        """Uses simple binning to calculate firing rate"""
        def __init__(self, neuron=None, trange=None, tres=100000):
            super(Neuron.BinRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'bin'
            self.tres = tres
        def calc(self):
            # make the start of the timepoints be an even multiple of self.tres. Round down to the nearest multiple. Do the same for the end of the timepoints. This way, timepoints will line up for different code objects
            tstart = self.trange[0] - (self.trange[0] % self.tres)
            tend   = self.trange[1] - (self.trange[1] % self.tres)
            t = np.arange( tstart, tend+self.tres, self.tres ) # t sequence demarcates left bin edges, add tres to trange[1] to make t end inclusive
            s = self.neuron.cut(trange) # spike times
            self.r, self.t = histogramSorted(self.neuron.spikes, bins=t) # assumes spikes are in chrono order
            self.r = self.r / float(self.tres) * 1000000 # spikes/sec
        def plot(self):
            super(Neuron.BinRate, self).plot()
            pl.title('neuron %d - binned spike rate' % self.neuron.id)

    class nISIRate(Rate):
        """Uses nisi inter spike intervals to calculate firing rate, with a fixed number of spikes nisi+1 per bin. nisi == 1 is the ISI rate"""
        def __init__(self, neuron=None, trange=None, nisi=3, tres=50000, interp='sah'):
            super(Neuron.nISIRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'nisi'
            self.nisi = nisi
            self.tres = tres
            self.interp = interp
        def calc(self):
            s = self.neuron.cut(self.trange) # spike times
            #n0 = self.n-1 # 0-based n
            # compare s to a shifted version of itself:
            diff = s[self.nisi::] - s[:-self.nisi:] # (n0 to end, single steps) - (beginning to n0 from end, single steps)
            r = float(self.nisi+1) / diff * 1000000 # spikes/sec
            t = s[self.nisi::] # for the corresponding timepoints, pick the spike time at the end of the group of n spikes to keep it causal
            #self.rawr = r
            #self.rawt = t
            if not self.interp:
                self.t = t
                self.r = r
                return
            # make the start of our interpolated timepoints be an even multiple of self.tres. Round down to the nearest multiple. This way, the timepoints will line up, even if different Rates have different starting points, like neuron.rate() vs experiment.rate()
            tstart = t[0] - (t[0] % self.tres)
            # should we have tend = t[-1] + self.tres ?
            self.t = np.arange(tstart, t[-1], self.tres) # new set of timepoints to interpolate over
            if self.interp == 'sah':
                self.r = sah(t, r, self.t, keep=False)
            elif self.interp == 'linear':
                f = sig.interpolate.interp1d(t, r, kind='linear') # returns an interpolation f'n
                self.r = f(self.t) # interpolate over the new timepoints
                # or maybe try sig.resample() instead
            else:
                raise ValueError, 'unknown interpolation method: %s' % self.interp

        def plot(self):
            super(Neuron.nISIRate, self).plot()
            pl.title('neuron %d - %d-inter-spike-interval spike rate, %s interpolation' % (self.neuron.id, self.nisi, self.interp))

    class WnISIRate(Rate):
        """Uses a weighted sum of various n inter spike intervals to calculate rate"""
        pass
        #def plot(self)

    class IDPRate(Rate):
        """Instantaneous discharge probability. See Pauluis and Baker, 2000"""
        def __init__(self, neuron=None, IDP_a=4, tres=50000):
            super(Neuron.IDPRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'idp'
            self.tres = tres
        def calc(self):

            # Step 1: find sudden changes in ISI distribution
            Is = diff(self.neuron.spikes) # intervals
            # y=gdtrc(a,b,x) returns the integral from x to infinity of the gamma probability density function.  SEE gdtr, gdtri
            b = IDP_a
            increaseOnsets = [] # stores the interval indices of significant sudden increases in ISI rate
            decreaseOnsets = [] # stores the interval indices of significant sudden decreases in ISI rate
            for (n,I) in enumerate(Is):
                try:
                    IDP_mu = Is(n+1) # the next interval
                    a = IDP_a / float(IDP_mu)
                    p1 = gdtrc(a,b,I) # integrate from I to infinity
                    IDP_mu = Is(n+2) # the interval after that
                    a = IDP_a / float(IDP_mu)
                    p2 = gdtrc(a,b,I) # integrate from I to infinity
                    if p1 < pthresh and p2 < pthresh: # then this nth interval was a sudden increase in ISI firing rate
                        inreaseOnsets.append(n)
                except IndexError:
                    pass
                try:
                    IDP_mu = Is(n-1) # the previous interval
                    a = IDP_a / float(IDP_mu)
                    p1 = gdtrc(a,b,I) # integrate from I to infinity
                    IDP_mu = Is(n-2) # the interval before that
                    a = IDP_a / float(IDP_mu)
                    p2 = gdtrc(a,b,I) # integrate from I to infinity
                    if p1 < pthresh and p2 < pthresh: # then this nth interval was a sudden increase in ISI firing rate
                        dereaseOnsets.append(n)
                except IndexError:
                    pass

            # Step 2: extend back half an ISI from sudden rate increase, and forward half an ISI from sudden rate decrease

            # Step 3: sample our new ISI rate with extended intervals with sufficient resolution
            isi = self.neuron.rate(kind='nisi', nisi=1, interp=None) # ISI rate object with with no interpolation
            isi.r
            isi.t
            isi = sah(isi)

            # use np.piecewise() ?
        pass
        #def plot(self)

    class GaussRate(Rate):
        """Uses a sliding Gaussian window to calculate firing rate"""
        def __init__(self, neuron=None, trange=None, width=200000):
            super(Neuron.GaussRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'gauss'
            self.width = width
        def calc(self):
            pass
        def plot(self):
            super(Neuron.GaussRate, self).plot()
            pl.title('Gaussian sliding window spike rate')

    class RectRate(Rate):
        """Uses a sliding rectangular window to calculate firing rate"""
        def __init__(self, neuron=None, trange=None, width=200000):
            super(Neuron.RectRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'rect'
            self.width = width
        def calc(self):
            pass
        def plot(self):
            super(Neuron.RectRate, self).plot()
            pl.title('rectangular sliding window spike rate')

    def ratepdf(self, **kwargs):
        """Returns an existing RatePDF object, or creates a new one if necessary"""
        try:
            self.ratepdfs
        except AttributeError: # self.ratepdfs doesn't exist yet
            self.ratepdfs = [] # create a list that'll hold RatePDF objects
        rpdf = Neuron.RatePDF(neuron=self, **kwargs)
        for ratepdf in self.ratepdfs:
            if rpdf == ratepdf: # need to define special == method for class RatePDF()
                return ratepdf # returns the first RatePDF object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self.ratepdfs
        rpdf.calc() # no matching RatePDF object was found, calculate it
        self.ratepdfs.append(rpdf) # add it to the RatePDF object list
        return rpdf

    class RatePDF(object):
        """Firing rate probability distribution function object"""
        def __init__(self, neuron=None, rrange=(0, 200), nbins=100, scale='log', normed=False, **kwargs): # rrange == rate range, ie limits of pdf x axis; nbins == number of rate bins
            self.neuron = neuron
            self.rrange = rrange
            self.nbins = nbins
            self.scale = scale
            self.normed = normed
            self.rate = neuron.rate(**kwargs) # returns a Rate object of some kind
        def __eq__(self, other):
            selfd = self.__dict__.copy()
            otherd = other.__dict__.copy()
            # Delete their n and r and logrrange attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
            [ d.__delitem__(key) for d in [selfd, otherd] for key in ['n', 'r', 'logrrange'] if d.has_key(key) ]
            if type(self) == type(other) and selfd == otherd:
                return True
            else:
                return False
        def calc(self):
            if self.scale == 'log':
                safe = list(self.rrange) # convert from immutable tuple to mutable list
                # prevent from taking log(0):
                for (i,r) in enumerate(safe):
                    if r == 0:
                        safe[i] = 0.1 # set to 0.1 Hz
                self.logrrange = np.log10(tuple(safe)) # convert back to tuple, is now safe to take log
                # r sequence demarcates left rate bin edges
                r = np.logspace(start=self.logrrange[0], stop=self.logrrange[1], num=self.nbins, endpoint=True, base=10.0)
            elif self.scale == 'linear':
                r = np.linspace(start=self.rrange[0], stop=self.rrange[1], num=self.nbins, endpoint=True)
            else:
                raise ValueError, 'Unknown scale: %s' % repr(scale)
            self.n, self.r = np.histogram(self.rate.r, bins=r, normed=True) # don't use histogram()'s normed, not sure what the hell it does
            #if self.normed:
            #    self.n = self.n / float(np.sum(self.n)) # do own normalization
        def plot(self):
            pl.figure()
            if self.scale == 'log':
                pl.axes().set_xscale(self.scale, basex=10)
                barwidth = list(np.diff(self.r)) # each bar will have a different width, convert to list so you can append
                # need to add one more entry to barwidth to the end to get nbins of them:
                #barwidth.append(barwidth[-1]) # not exactly correct
                logbinwidth = (self.logrrange[1]-self.logrrange[0]) / float(self.nbins)
                barwidth.append(10**(self.logrrange[1]+logbinwidth)-self.r[-1]) # should be exactly correct
            elif self.scale == 'linear':
                pl.axes().set_xscale(self.scale)
                barwidth = (self.rrange[1]-self.rrange[0]) / float(self.nbins)
            else:
                raise ValueError, 'Unknown scale: %s' % repr(scale)
            #pl.hist(self.n, bins=self.r, normed=0, bottom=0, width=None, hold=False) # doesn't seem to work
            pl.bar(left=self.r, height=self.n, width=barwidth, bottom=0, color='k', yerr=None, xerr=None, ecolor='k', capsize=3)
            pl.title('neuron %d - %s spike rate PDF' % (self.neuron.id, self.rate.kind))
            if self.normed:
                pl.ylabel('probability')
            else:
                pl.ylabel('count')
            pl.xlabel('spike rate')


    def raster(self):
        # use pylab.vlines()
        pass


class Movie(Dimstim.Movies.Movie): # inherit from Dimstim Movie() class (assumes it's new-style)
    """A Movie stimulus object"""
    def __init__(self, name=None, path=DEFAULTMOVIEPATH, parent=None):
        super(Movie, self).__init__() # first run __init__() of inherited Dimstim Movie class
        self.level = 5 # level in the hierarchy
        try:
            self.e = parent() # init parent Experiment object
        except TypeError: # parent is an instance, not a class
            self.e = parent # save parent Experiment object
        self.name = name
        self.path = path
    def load(self):
        # Load movie data
        f = file(self.path + self.name, 'rb') # open the movie file for reading in binary format
        headerstring = f.read(5)
        if headerstring == 'movie': # a header has been added to the start of the file
            (self.ncellswide,) = struct.unpack('H', f.read(2)) # 'H'== unsigned short int == 2 bytes on this PC
            (self.ncellshigh,) = struct.unpack('H', f.read(2))
            (self.nframes,) = struct.unpack('H', f.read(2))
            if self.nframes == 0: # this was used in Cat 15 mseq movies to indicate 2**16 frames
                self.nframes = 2**16
            offset = 11 # header is this long
        else: # there's no header at the start of the file, set the file pointer back to the beginning and use these hard coded values:
            f.seek(0)
            self.ncellswide = self.ncellshigh = 64
            self.nframes = 6000
            offset = 0 # header is this long
        # read in all of the frames
        self.data = np.fromfile(f, np.uint8, count=self.nframes*self.ncellshigh*self.ncellswide).reshape(self.nframes,self.ncellshigh,self.ncellswide) # read it all in
        #self.data = numarray.fromfile(f, np.UInt8, (self.nframes,self.ncellshigh,self.ncellswide)) # read it all in
        leftover = f.read() # check if there are any leftover bytes in the file
        if leftover != '':
            pprint(leftover)
            print self.nframes,self.ncellshigh,self.ncellswide
            raise RuntimeError, 'There are unread bytes in movie file %s. Width, height, or nframes is incorrect in the movie file header.' % repr(self.name)
        #self.data = self.data[::,::-1,::] # flip the movie frames vertically for OpenGL's bottom left origin
        f.close() # close the movie file




# init some typical movies (but don't load 'em til needed). Then, just point to them within the appropriate Experiments
MSEQ32 = Movie(name='mseq32.m', parent=None)
MSEQ16 = Movie(name='mseq16.m', parent=None)
# shouldn't use sparse bar movies anymore, can access VisionEgg directly now, get the framebuffers to directly do STA
#sparsebars = Movie(path='C:/data/Cat 15/Track 7c/72 - track 7c sparseexps/', name='72 - track 7c sparseexps.sparsebars.movie');


# init and load some neuropy objects:
#print 'Initing and loading Recording(92):'
#r92=Recording(92)
#r92.load()
#print 'Initing and loading Track(\'7c\'):'
#t=Track('7c')
#t.load
print 'Initing and loading Recording(71):'
r71=Recording(71)
r71.load()


