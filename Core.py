"""Core neuropy functions and classes"""

print 'importing Core'

DEFAULTDATAPATH = 'C:/data/' # the convention in neuropy is that all 'path' var names have a trailing slash
DEFAULTMODELPATH = 'C:/model/'
DEFAULTCATID = 15
DEFAULTSYSTEMNAME = 'Cat 15'
DEFAULTTRACKID = '7c'
RIPKEYWORDS = ['best'] # a Rip with one of these keywords (listed in decreasing priority) will be loaded as the default Rip for its Recording/Run
SLASH = '/' # use forward slashes instead of having to use double backslashes
TAB = '    ' # 4 spaces

DEFAULTMOVIEPATH = 'C:/pub/Movies/'
DEFAULTMOVIENAME = 'mseq32.m'

DEFAULTCODEWORDLENGTH = 10

import os
import types
from pprint import pprint
import struct
import re
import StringIO
import sys
import random

import numpy as np
import pylab as pl
import matplotlib as mpl
import scipy as sp
import scipy.signal as sig
from numpy import arange, array, asarray, log, log10, rand, randn, zeros, ones, diff, concatenate, concatenate as cat, histogram
from pylab import figure, plot, loglog, hist, bar, barh, xlabel, ylabel, xlim, ylim, title, gcf, gca, get_current_fig_manager as gcfm, axes, axis, hold, imshow

# Nah!: Rips should really have ids to make them easier to reference to: r[83].rip[0] instead of r[83].rip['conservative spikes'] - this means adding id prefixes to rip folder names (or maybe suffixes: 'conservative spikes.0.rip', 'liberal spikes.1.rip', etc...). Prefixes would be better cuz they'd force sorting by id in explorer (which uses alphabetical order) - ids should be 0-based of course
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
    x = array(a, copy=False)
    y = array(b, copy=False)
    print x.shape
    print y.shape
    return np.less(np.absolute(x-y), atol + rtol * np.absolute(y))

def histogramSorted(sorteda, bins=10, range=None):
    """Builds a histogram, stolen from numpy.histogram(), modified to assume sorted input"""
    a = asarray(sorteda).ravel()
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
    n = cat([n, [len(a)]]) # this adds a bin that includes overflow points
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
        di = diff(i).nonzero()[0] # find changes in i, nonzero() method returns a tuple, pick the result for the first dim with [0] index
        si = approx(t[1::], ts[di]) # check at those change indices if t ~= ts (ignoring potential floating point representational inaccuracies). If so, inc i at that point so you keep y at that point.
        #print i
        i[di[si]] += 1
        #print i
    return y[i]

def corr(x,y):
    """Returns correlation of signals x and y. This should be equivalent to np.corrcoef(),
    but that one doesn't seem to work for signals with zeros in them. Check how std() works exactly"""
    x = array(x)
    y = array(y)
    return ((x * y).mean() - x.mean() * y.mean()) / (x.std() * y.std())

def getargstr(obj):
    import inspect
    argstr = apply(inspect.formatargspec, inspect.getargspec(obj))
    if inspect.isfunction(obj):
        pass
    elif inspect.ismethod(obj):
        # stolen from wx.py.introspect.getCallTip:
        temp = argstr.split(',')
        if len(temp) == 1:  # No other arguments.
            argstr = '()'
        elif temp[0][:2] == '(*': # first param is like *args, not self
            pass
        else:  # Drop the first argument.
            argstr = '(' + ','.join(temp[1:]).lstrip()
    else:
        argstr = '()'
    return argstr

def binaryarray2int(bin):
    """Takes a binary array (only 1s and 0s) and returns the base 10 integer representations"""
    #assert type(bin) == type(array)
    nbits = bin.shape[0] # length of the highest (first) dimension (the rows)
    nd = bin.ndim
    multiplier = []
    for i in range(nbits):
        multiplier.append(2**i)
    multiplier = array(multiplier, ndmin=nd).transpose()
    #print multiplier
    x = bin*multiplier
    #print x
    return x.sum(axis=0) # sum over the lowest dimension (the columns)

def getbinarytable(nbits=8):
    """Generates a 2D binary table containing all possible words for nbits, with bits in the rows and words in the columns (msb at bottomest row)"""
    rowlength = 2**nbits
    '''
    x = zeros((nbits, 2**nbits)) # init an array
    for bit in range(nbits):
        pattern = [0]*2**bit
        pattern.extend([1]*2**bit)
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        row = pattern*npatterns
        x[bit]=row
    return x
    '''
    '''
    x = zeros((nbits, 2**nbits), dtype=np.int8) # init an array
    for bit in range(nbits): # one row at a time
        pattern = array(0, dtype=np.int8).repeat(2**bit)
        pattern = cat((pattern, array(1, dtype=np.int8).repeat(2**bit)))
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        row = np.repmat(pattern, 1, npatterns)
        x[bit::,::] = row
    return x
    '''
    # this seems to be the fastest method:
    x = []
    for bit in range(nbits): # one row at a time
        pattern = array(0, dtype=np.int8).repeat(2**bit)
        pattern = cat((pattern, array(1, dtype=np.int8).repeat(2**bit)))
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        row = np.repmat(pattern, 1, npatterns)
        x.append(row)
    return cat(x)

def shuffle(ip):
    """Takes an input list and returns a shuffled (without replacement) copy, its only benefit
    over and above random.sample() is that you don't have to pass a second argument len(ip)
    every time you use it"""
    return random.sample(ip, len(ip))

def randomize(ip):
    """Takes an input list and returns a randomized (with replacement) output list of
    the same length, sampled from the input sequence"""
    op = [] # init output list
    for i in range(0, len(ip)):
        op.append(random.choice(ip))
    return op


class Data(object): # use 'new-style' classes
    """Abstract data class. Data can have multiple Cats"""
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
    """Abstract model class. Model can have multiple model Systems"""
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
        from Recording import Recording
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
