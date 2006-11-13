"""Defines base neuropy Data and Model classes, as well as other miscellaneous functions and classes"""

#print 'importing Core'

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

DEFAULTCODETRES = 20000 # us
DEFAULTCODEWORDLENGTH = 10 # in bits

import os
import sys
import time
import types
import __main__
import struct
import re
import StringIO
import random
import math

from copy import copy
from pprint import pprint
printraw = sys.stdout.write # useful for raw printing

import numpy as np
import pylab as pl
import matplotlib as mpl
import scipy as sp
import scipy.signal as sig
import scipy.weave as weave
from numpy.random import rand, randn, randint
from numpy import arange, array, array as ar, asarray, log, log10, zeros, ones, diff, concatenate, concatenate as cat
from pylab import figure, plot, loglog, hist, bar, barh, xlabel, ylabel, xlim, ylim, title, gcf, gca, get_current_fig_manager as gcfm, axes, axis, hold, imshow
import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

mpl.use('WXAgg')
mpl.interactive(True)


class Data(object): # use 'new-style' classes
    """Abstract data class. Data can have multiple Cats"""
    def __init__(self, dataPath=DEFAULTDATAPATH):
        self.level = 0 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.name = 'Data'
        self.path = dataPath
        self.c = {} # store Cats in a dictionary
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        # Data has no parent to write to
    def load(self):

        from Cat import Cat

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        catNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Cat ') ] # os.listdir() returns all dirs AND files
        for catName in catNames:
            cat = Cat(id=None, name=catName, parent=self) # make an instance using just the catName (let it figure out the cat id)
            cat.load() # load the Cat
            self.c[cat.id] = cat # save it
        #if len(self.c) == 1:
        #   self.c = self.c.values[0] # pull it out of the dictionary

_data = Data() # init a default Data object to use as a container for everything that falls under the data object hierarchy


class Model(Data):
    """Abstract model class. Model can have multiple model Systems"""
    def __init__(self, modelPath=DEFAULTMODELPATH):
        self.level = 0 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.name = 'Model'
        self.path = modelPath
        self.s = {} # store model Systems in a dictionary
    def load(self):

        from System import System

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        systemNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) ] # os.listdir() returns all dirs AND files
        for systemName in systemNames:
            system = System(name=systemName, parent=self) # make an instance using the systemName
            system.load() # load the System
            self.s[system.name] = system # save it
        #if len(self.s) == 1:
        #   self.s = self.s.values[0] # pull it out of the dictionary

_model = Model() # init a default Model object to use as a container for everything that falls under the model object hierarchy


def getargstr(obj):
    """Returns object's argument list as a string. Stolen from wx.py package?"""
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

def barefigure(*args, **kwargs):
    """Creates a bare figure with no toolbar or statusbar"""
    figure(*args, **kwargs)
    gcfm().frame.GetStatusBar().Hide()
    gcfm().frame.GetToolBar().Hide()
barefigure.__doc__ += '\n' + figure.__doc__

def lastcmd():
    """Returns a string containing the last command entered at the PyShell prompt.
    Maybe this could be extended to work with other shells too?"""
    try:
        # PyShell's shell.py was hacked to save the last command
        # as an attrib in Shell.push()
        return __main__.shell.lastcmd
    except AttributeError:
        return 'unknown'

class CanvasFrame(wx.Frame):
    """A minimal wx.Frame containing a matplotlib figure"""
    def __init__(self, title='frame', size=(550,350)):
        wx.Frame.__init__(self, None, -1, title=title, size=size)
        self.SetBackgroundColour(wx.NamedColor("WHITE"))
        self.figure = mpl.figure.Figure(figsize=(5,4), dpi=100)
        #self.axes = self.figure.add_subplot(111)
        #t = arange(0.0, 3.0, 0.01)
        #s = sin(2*pi*t)
        #self.axes.plot(t,s)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)

        # Capture the paint message, slows frame down a little, can be commented out
        #wx.EVT_PAINT(self, self.OnPaint)
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        self.SetSizer(self.sizer)
        self.Fit()
    def OnPaint(self, event):
        self.canvas.draw()

def frame(**kwargs):
    """Returns a CanvasFrame object"""
    frame = CanvasFrame(**kwargs)
    frame.Show(True)
    return frame
frame.__doc__ += '\n' + CanvasFrame.__doc__
frame.__doc__ += '\n\n**kwargs:\n' + getargstr(CanvasFrame.__init__)


class ReceptiveFieldFrame(wx.Frame):
    """A wx.Frame for plotting a scrollable 2D grid of receptive fields, with neuron and time labels
    rfs is a list of (nt, width, height) sized receptive fields made up of uint8 RGB data"""
    def __init__(self, parent=None, id=-1, title='ReceptiveFieldFrame', rfs=None, neurons=None, t=None, scale=2.0, **kwargs):
        self.rfs = rfs
        self.neurons = neurons
        self.t = t
        self.title = title
        kwargs['style'] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, parent=parent, id=id, title=title, **kwargs)
        self.panel = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)

        self.bitmaps = {}
        for ni, n in enumerate(self.neurons):
            self.bitmaps[ni] = {}
            for ti, t in enumerate(self.t):
                rf = self.rfs[ni][ti]
                im = wx.ImageFromData(width=rf.shape[0], height=rf.shape[1], data=rf.data) # expose rf as databuffer
                im = im.Scale(width=im.GetWidth()*scale, height=im.GetHeight()*scale)
                self.bitmaps[ni][t] = wx.StaticBitmap(parent=self.panel, bitmap=im.ConvertToBitmap())

        #self.Bind(wx.EVT_PAINT, self.OnPaint)
        #self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
        self.__set_properties()
        self.__do_layout()
    '''
    def OnPaint(self, event):
        #self.canvas.draw()
        event.Skip()
    '''
    def OnMouseWheel(self, event):
        pass
    def __set_properties(self):
        self.SetTitle(self.title)
        self.panel.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.panel.SetScrollRate(10, 10)
    def __do_layout(self):
        sizer_1 = wx.GridSizer(1, 1, 0, 0)
        grid_sizer_1 = wx.FlexGridSizer(rows=len(self.neurons)+1, cols=len(self.t)+1, vgap=2, hgap=2) # add an extra row and column for the text labels
        grid_sizer_1.Add((1, 1), 0, wx.ADJUST_MINSIZE, 0) # spacer in top left corner
        for t in self.t:
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "%sms" % t), 0, wx.ADJUST_MINSIZE, 0) # text row along top
        for ni, n in enumerate(self.neurons):
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "n%d" % n.id), 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ADJUST_MINSIZE, 0) # text down left side
            for t in self.t:
                grid_sizer_1.Add(self.bitmaps[ni][t], 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL, 0)
        self.panel.SetAutoLayout(True)
        self.panel.SetSizer(grid_sizer_1)
        grid_sizer_1.Fit(self.panel)
        #grid_sizer_1.SetSizeHints(self.panel) # prevents the panel from being resized to something smaller than the above fit size
        '''
        # might be a more direct way to set these:
        for rowi in range(1, len(self.ns)+1):
            print 'rowi:', rowi
            grid_sizer_1.AddGrowableRow(rowi)
        for coli in range(1, len(self.ts)+1):
            print 'coli:', coli
            grid_sizer_1.AddGrowableCol(coli)
        '''
        sizer_1.Add(self.panel, 1, wx.ADJUST_MINSIZE|wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        #sizer_1.SetSizeHints(self) # prevents the frame from being resized to something smaller than the above fit size
        self.Layout()

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
        fo.write( struct.pack('@qq', int(line[0]), int(line[1])) ) # read both values in as a C long longs, using the system's native ('@') byte order
    fi.close()
    fo.close()
    print 'Converted ascii din: ' + repr(fin) + ' to binary din: ' + repr(fout)

def convertalltxtdin2binarydin(path=None):
    """Converts all text .csv din files in path (or cwd) to 64 bit binary .din files of the same name"""
    if path == None:
        path = os.getcwd()

    listing = os.listdir(path)
    dinfnames = []

    for fname in listing:
        if fname.endswith('.csv'):
            dinfnames.append(fname.rstrip('.csv')) # text din filenames without the .csv extension

    for dinfname in dinfnames:
        fin = path + '\\' + dinfname + '.csv'
        fout = path + '\\' + dinfname + '.din'
        #os.rename(fout, fin) # rename the csv .din file to .din.txt extension
        txtdin2binarydin(fin, fout) # convert the text .csv file to binary .din file
        #os.remove(fin) # delete the .din.txt file

def renameSpikeFiles(path, newname):
    """Renames all .spk files in path to newname, retaining their '_t##.spk' ending"""
    for fname in os.listdir(path):
        if fname.endswith('.spk'):
            i = fname.find('_t')
            if i != -1:
                newfname = newname+fname[i::]
                print newfname
                os.rename(path+SLASH+fname, path+SLASH+newfname)

def warn(msg, level=2, exit_val=1):
    """Standard warning printer. Gives formatting consistency. Stolen from IPython.genutils"""
    if level > 0:
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
def iterable(x):
    """Check if the input is iterable, stolen from numpy.iterable()"""
    try:
        iter(x)
        return True
    except:
        return False

def toiter(x):
    """Convert to iterable. If input is iterable, returns it. Otherwise returns it in a list"""
    if iterable(x):
        return x
    else:
        return [x]

def tolist(x):
    """Convert to list. If input is a dict, returns its values. If it's already a list, returns it.
    Otherwise, input is returned in a list."""
    if x.__class__ == dict:
        return list(x.values())
    elif x.__class__ == list:
        return x
    else:
        return [x] # stick it in a list

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
    #print x.shape
    #print y.shape
    return np.less(np.absolute(x-y), atol + rtol * np.absolute(y))

def histogram(a, bins=10, range=None, normed=False):
    """Builds a histogram, stolen from numpy.histogram(), modified to allow
    normed='pdf' or normed='pmf' (prob mass function)"""
    a = asarray(a).ravel()
    if not iterable(bins):
        if range is None:
            range = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins, endpoint=False)
    n = np.sort(a).searchsorted(bins)
    n = concatenate([n, [len(a)]])
    n = n[1:]-n[:-1]
    if normed:
        if normed == 'pdf':
            db = bins[1] - bins[0]
            return 1.0/(a.size*db) * n, bins
        elif normed == 'pmf':
            return n/float(sum(n)), bins
    else:
        return n, bins

def histogramSorted(sorteda, bins=10, range=None, normed=False):
    """Builds a histogram, stolen from numpy.histogram(), modified to assume
    sorted input and to allow normed='pdf' or normed='pmf' (prob mass function)"""
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
    n = concatenate([n, [len(a)]]) # this adds a bin that includes overflow points
    n = n[1:]-n[:-1] # subtracts a shifted version of itself
    if normed:
        if normed == 'pdf':
            db = bins[1] - bins[0]
            return 1.0/(a.size*db) * n, bins
        elif normed == 'pmf':
            return n/float(sum(n)), bins
    else:
        return n, bins

def sah(t, y, ts, keep=False):
    """Resample using sample and hold. Returns resampled values at ts given the original points (t,y)
    such that the resampled values are just the most recent value in y (think of a staircase with non-uniform steps).
    Assumes that t is sorted. t and ts arrays should be of the same data type. Contributed by Robert Kern."""
    i = np.searchsorted(t, ts) - 1 # find where ts falls in t, dec so you get indices that point to the most recent value in y
    i = np.where(i < 0, 0, i) # handle the cases where ts is smaller than the first point.
    '''this has an issue of not keeping the original data point where ts == t'''

    ###NOTE: can probably get around having to do this by using searchsorted's new 'side' keyword

    if keep:
        # The following ensures that the original data point is kept when ts == t, doesn't really work if the shortest ISI is less than tres in ts
        di = diff(i).nonzero()[0] # find changes in i, nonzero() method returns a tuple, pick the result for the first dim with [0] index
        si = approx(t[1::], ts[di]) # check at those change indices if t ~= ts (ignoring potential floating point representational inaccuracies). If so, inc i at that point so you keep y at that point.
        #print i
        i[di[si]] += 1
        #print i
    return y[i]

def corrcoef(x, y):
    """Returns the correlation coefficient of signals x and y. This just uses np.corrcoef(),
    but converts to floats first, cuz np.corrcoef() seems to have issues with integer signals,
    especially those with zeros in them."""
    #assert len(x) == len(y), 'arrays need to be of equal length'
    x = np.float64(x)
    y = np.float64(y)
    return np.corrcoef(x, y)[0, 1] # pick one of the 2 entries in the correlation coefficient matrix, on the -ve diagonal (er, the one that goes from bottom left to top right, that's what I mean)
    #return ((x * y).mean() - x.mean() * y.mean()) / (x.std() * y.std()) # this works just fine as well, easier to understand too

def bin(i, minbits=8):
    """Return the binary representation of an integer, or sequence of integers.
    If necessary, will append leading zeros if result is less than minbits long.
    Uses np.binary_repr()"""
    ints = toiter(i) # ensure it's iterable
    sints = []
    for i in ints:
        s = np.binary_repr(i)
        nzerostoadd = minbits - len(s) # OK if this is -ve
        s = '0'*nzerostoadd + s # add enough leading zeros to get requested minbits
        sints.append(s)
    if len(sints) == 1:
        sints = sints[0] # pull it out of the list
    return sints

def binslow(i, minbits=8):
    """Return the binary representation of an integer. If necessary, will append leading zeros
    if result is less than minbits long. Seems like np.binary_repr() is a somewhat faster alternative.
    First 2 lines stolen from Andrew Gaul <andrew@gaul.org> off the web"""
    l = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    s = ''.join(map(lambda x, l=l: l[int(x, 16)], hex(i)[2:]))
    s = s.lstrip('0') # strip s of leading zeros
    nzerostoadd = minbits - len(s)
    s = '0'*nzerostoadd + s # add enough leading zeros to get requested minbits
    return s

def binaryarray2int(bin):
    """Takes a 2D binary array (only 1s and 0s, with rows LSB to MSB from top to bottom)
    and returns the base 10 integer representations of the columns"""
    #assert type(bin) == type(array)
    nbits = bin.shape[0] # length of the highest (first) dimension, ie the number of rows
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

def shuffle(x):
    """Takes an input list x and returns a shuffled (without replacement) copy. Its only benefit
    over and above random.sample() is that you don't have to pass a second argument len(x)
    every time you use it.
    In NumPy, it's better (and faster) to use np.random.shuffle()"""
    return random.sample(x, len(x))

def randomize(x):
    """Takes an input list x and returns a randomized (with replacement) output list of
    the same length, sampled from the input sequence"""
    y = [] # init output list
    for i in range(0, len(x)):
        y.append(random.choice(x))
    return y

def fact(n):
    """Factorial!"""
    assert n.__class__ == int
    assert n >= 0
    if n == 0:
        n = 1 # 0! == 1!
    result = n
    for i in range(1, n):
        result *= i
    return result

def nPr(n, r):
    """n Pick r"""
    return fact(n) / fact(n-r)

def nCr(n, r):
    """n Choose r"""
    return nPr(n, r) / fact(r)

ncr = nCr # convenience f'ns
npr = nPr

'''
class Chooser(object):
    """Container class for combination generator.
    Takes a sequence of objects and the number of objects to return in each combination"""
    def __init__(self, objects, r):
        self.objects = asarray(objects)
        self.n = len(objects)
        self.r = r
        self.i = asarray([0]*self.r) # stores all the current i values for all r levels of nested for loops
    def combgen(self, level=0):
        """Generator that yields all possible combinations of self.objects, without replacement.
        Eg, if self.objects=[0,1,2] and self.r=2, this yields [0,1], [0,2], and [1,2], one at a time.
        A recursive generator is used in order to create the necessary r number of nested for loops"""
        try: # recursive case
            if level == 0: # handles special case for starting index of top level for loop
                starti = 0
            else:
                starti = self.i[level-1] + 1 # start this level's loop at one greater than the previous level's current loop index
            for self.i[level] in range(starti, self.n+1): # not too sure why this is n+1, but it works
                for comb in self.combgen(level=level+1): # iterate over next level's generator
                    yield comb # yield whatever next level (level+1) returns, back up to previous level (level-1)
        except: # base case, we're at the deepest recursion level (innermost for loop)
            yield self.objects[self.i] # use the current index state for all levels to yield a combination of objects

def choose(objects, r):
    co = Chooser(objects, r)
    combs = list(co.combgen())
    return combs
'''
def combgen(objects, r=2, i=None, level=0):
    """Generator that yields, without replacement, all length r possible combinations of objects from a length n sequence.
    Eg, if objects=[0,1,2] and r=2, this yields [0,1], [0,2], and [1,2], one at a time.
    A recursive generator is used in order to create the necessary r number of nested for loops"""
    assert r <= len(objects)
    try: # recursive case
        if i == None:
            i = asarray([0]*r) # stores all the current index values for all r nested for loops
        if level == 0: # handles special case for starting index of top level for loop
            starti = 0
        else:
            starti = i[level-1] + 1 # start this level's loop index at one greater than the previous level's current loop index
        for i[level] in range(starti, len(objects)+1): # not too sure why this is n+1, but it works
            for comb in combgen(objects, r=r, i=i, level=level+1): # iterate over next level's generator
                yield comb # yield whatever the next level (level+1) yields, pass it on up to the previous level (level-1)
    except: # base case, we're at the deepest recursion level (innermost for loop)
        yield asarray(objects)[i] # use the current index state for all levels to yield a combination of objects

'''
# this f'n isn't really needed, just use objlist.sort(key=lambda obj: obj.attrib)
def sortby(objs, attrib, cmp=None, reverse=False):
    """Returns objects list sorted according to the specified object attribute.
    attrib should be passed as a string"""
    objs.sort(key=lambda obj: obj.__getattribute__(attrib), cmp=cmp, reverse=reverse) # sort in-place
    return objs
'''

def mean_accum(data):
    """Takes mean by accumulating over 0th axis in data,
    much faster than numpy's mean() method because it avoids making any copies of the data
    Suggested by Tim Hochberg"""
    result = np.zeros(data[0].shape, np.float64) # init output array
    for dataslice in data: # this for loop isn't such a bad thing cuz the massive add step inside the loop is the limiting factor
        result += dataslice
    result /= len(data)
    return result

def mean_accum2(data, indices):
    """A variant of mean_accum(), where you provide all the data and the indices into it
    to average over. This was Tim Hochberg's version"""
    result = np.zeros(data[0].shape, np.float64)
    for i in indices:
        result += data[i]
    result /= len(indices)
    return result


class neuropyScalarFormatter(mpl.ticker.ScalarFormatter):
    """Overloaded from mpl.ticker.ScalarFormatter for 4 reasons:
    1) turn off stupid offset
    2) increase maximum possible number of sigfigs
    3) increase +ve and -ve order of magnitude thresholds before switching to scientific notation
    4) keep exponents in engineering notation, ie multiples of 3
    """
    def __init__(self, useOffset=False, useMathText=False):
        # useOffset allows plotting small data ranges with large offsets:
        # for example: [1+1e-9,1+2e-9,1+3e-9]
        # useMathText will render the offset an scientific notation in mathtext
        #super(neuropyScalarFormatter, self).__init__(useOffset=useOffset, useMathText=useMathText) # can't use this, cuz derived from an old-style class
        mpl.ticker.ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)

    def _set_orderOfMagnitude(self, range):
        # if scientific notation is to be used, find the appropriate exponent
        # if using an numerical offset, find the exponent after applying the offset
        locs = np.absolute(self.locs)
        if self.offset: oom = math.floor(math.log10(range))
        else:
            if locs[0] > locs[-1]: val = locs[0]
            else: val = locs[-1]
            if val == 0: oom = 0
            else: oom = math.floor(math.log10(val))
        if oom < -3: # decreased -ve threshold for sci notation
            self.orderOfMagnitude = (oom // 3)*3 # stick to engineering notation, multiples of 3
        elif oom > 6: # increased +ve threshold for sci notation
            self.orderOfMagnitude = (oom // 3)*3 # stick to engineering notation, multiples of 3
        else:
            self.orderOfMagnitude = 0

    def _set_format(self):
        # set the format string to format all the ticklabels
        locs = (array(self.locs)-self.offset) / 10**self.orderOfMagnitude+1e-15
        sigfigs = [len(str('%1.10f'% loc).split('.')[1].rstrip('0')) for loc in locs] # '%1.3f' changed to '%1.10f' to increase maximum number of possible sigfigs
        sigfigs.sort()
        self.format = '%1.' + str(sigfigs[-1]) + 'f'
        if self._usetex or self._useMathText: self.format = '$%s$'%self.format


class neuropyAutoLocator(mpl.ticker.MaxNLocator):
    """A tick autolocator that generates more ticks than the standard mpl autolocator"""
    def __init__(self):
        #mpl.ticker.MaxNLocator.__init__(self, nbins=9, steps=[1, 2, 5, 10]) # standard autolocator
        mpl.ticker.MaxNLocator.__init__(self) # use MaxNLocator's defaults instead


def entropy(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob values in p"""
    p = asarray(p)
    assert approx(sum(p), 1.0) # make sure the probs sum to 1
    return -sum(p*np.log2(p))


class Ising(object):
    """Ising maximum entropy model"""
    def __init__(self, means, pairmeans, algorithm='CG'):

        from scipy import maxentropy

        nbits = len(means)
        npairs = len(pairmeans)
        assert npairs == nCr(nbits, 2) # sanity
        self.intsamplespace = range(0, 2**nbits)
        table = getbinarytable(nbits=nbits) # words are in the columns, MSB at bottom row
        self.binsamplespace = ar([ table[::-1, wordi] for wordi in range(0, 2**nbits) ]) # all possible binary words (each is MSB to LSB), as arrays of 0s and 1s
        self.samplespace = self.binsamplespace * 2 - 1 # convert 0s to -1s
        # return the i'th bit (LSB to MSB) of binary word x
        f1s = [ lambda x, i=i: x[-1-i] for i in range(0, nbits) ] # have to do i=i to statically assign value (gets around scope closure problem)
        # return product of the i'th and j'th bit (LSB to MSB) of binary word
        f2s = [ lambda x, i=i, j=j: x[-1-i] * x[-1-j] for i in range(0, nbits) for j in range(i+1, nbits) ]
        f = cat((f1s, f2s))
        self.model = maxentropy.model(f, self.samplespace)
        #self.model.mindual = -10000
        #self.model.log = None # needed to make LBFGSB algorithm work
        # Now set the desired feature expectations
        means = asarray(means)
        pairmeans = asarray(pairmeans)
        #pairmeans /= 2.0 # add the one half in front of each coefficient, NOT TOO SURE IF THIS SHOULD GO HERE! causes convergence problems
        K = cat((means, pairmeans))
        self.model.verbose = False

        # Fit the model
        self.model.fit(K, algorithm=algorithm)

        self.hi = self.model.params[0:nbits]
        self.Jij = self.model.params[nbits:nbits+npairs]
        self.p = self.model.probdist()
        assert (len(self.hi), len(self.Jij), len(self.p)) == (nbits, npairs, 2**nbits) # sanity checks
        #print 'means:', means
        #print 'pairmeans:', pairmeans
        print 'niterations:', self.model.iters
        print 'hi:', self.hi.__repr__()
        print 'Jij:', self.Jij.__repr__()
        '''
        # Output the distribution
        print "\nFitted model parameters are:\n" + str(self.model.params)
        print "\nFitted distribution is:"
        for j in range(len(self.model.samplespace)):
            x = ar(self.model.samplespace[j])
            x = (x+1)/2 # convert from -1s and 1s back to 0s and 1s
            print '\tx:%s, p(x):%s' % (x, p[j])
        '''
        '''
        # Now show how well the constraints are satisfied:
        print
        print "Desired constraints:"
        print "\tp['dans'] + p['en'] = 0.3"
        print ("\tp['dans'] + p['" + a_grave + "']  = 0.5").encode('utf-8')
        print
        print "Actual expectations under the fitted model:"
        print "\tp['dans'] + p['en'] =", p[0] + p[1]
        print ("\tp['dans'] + p['" + a_grave + "']  = " + str(p[0]+p[2])).encode('utf-8')
        # (Or substitute "x.encode('latin-1')" if you have a primitive terminal.)
        '''
