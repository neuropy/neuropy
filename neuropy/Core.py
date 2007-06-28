"""Defines base neuropy Data and Model classes, as well as other miscellaneous functions and classes"""

#print 'importing Core'

import os
import sys
import time
import types
import __main__
import struct
import re
import cStringIO
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
from numpy import arange, array, array as ar, asarray, asarray as aar, log, log2, log10, sqrt, zeros, ones, diff, concatenate, concatenate as cat, mean, median, std
from numpy.core.ma import array as mar
from pylab import figure, plot, loglog, hist, bar, barh, xlabel, ylabel, xlim, ylim, title, gcf, gca, get_current_fig_manager as gcfm, axes, axis, hold, imshow
import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

mpl.use('WXAgg')
mpl.interactive(True)

from dimstim.Core import dictattr # Dictionary with attribute access

# global DEFAULTS

DATAPATH = os.path.join(os.sep, 'data')
MODELPATH = os.path.join(os.sep, 'model')

SPECIES = 'Cat'
CATID = 15
RATID = 0
if SPECIES == 'Cat':
    ANIMALNAME = '%s %d' % (SPECIES, CATID)
elif SPECIES == 'Rat':
    ANIMALNAME = '%s %d' % (SPECIES, RATID)
else:
    raise ValueError, 'unknown species %s' % SPECIES

TRACKID = '7c'
RIPKEYWORDS = ['best'] # a Rip with one of these keywords (listed in decreasing priority) will be loaded as the default Rip for its Recording/Run
MOVIEPATH = os.path.join(os.sep, 'pub', 'Movies')
MOVIENAME = 'mseq32.m'

SYSTEMNAME = 'example model system'

CODEKIND = 'binary'
CODETRES = 20000 # us
CODEWORDLEN = 10 # in bits

TAB = '    ' # 4 spaces


class Data(object):
    """Abstract data class. Data can have multiple Animals in it"""
    def __init__(self, dataPath=DATAPATH):
        self.level = 0 # level in the hierarchy
        self.treebuf = cStringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.name = 'Data'
        self.path = dataPath
        self.a = dictattr() # store Animals in a dictionary with attrib access
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer"""
        self.treebuf.write(string)
        # Data has no parent to write to
    def load(self):

        from Animal import Cat, Rat

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        dirnames = [ dirname for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname)) ]
        for dirname in dirnames:
            if dirname.startswith('Cat '):
                cat = Cat(name=dirname, parent=self) # make an instance using just the dirname
                cat.load() # load the Cat
                self.a[cat.name] = cat # save it, using its (dir)name as the dict key
            elif dirname.startswith('Rat '):
                rat = Rat(name=dirname, parent=self) # make an instance using just the dirname
                rat.load() # load the Rat
                self.a[rat.name] = rat # save it, using its (dir)name as the dict key

_data = Data() # init a default Data object to use as a container for everything that falls under the data object hierarchy


class Model(Data):
    """Abstract model class. Model can have multiple model Systems"""
    def __init__(self, modelPath=MODELPATH):
        self.level = 0 # level in the hierarchy
        self.treebuf = cStringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.name = 'Model'
        self.path = modelPath
        self.s = dictattr() # store model Systems in a dictionary with attrib access
    def load(self):

        from System import System

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        systemNames = [ dirname for dirname in os.listdir(self.path)
                        if os.path.isdir(os.path.join(self.path, dirname)) ]
        for systemName in systemNames:
            system = System(name=systemName, parent=self) # make an instance using the systemName
            system.load() # load the System
            self.s[system.name] = system # save it

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

def innerclass(cls):
    '''Class decorator for making a class behave as a Java (non-static) inner
    class.

    Each instance of the decorated class is associated with an instance of its
    enclosing class. The outer instance is referenced implicitly when an
    attribute lookup fails in the inner object's namespace. It can also be
    referenced explicitly through the property '__outer__' of the inner
    instance.

    Title: Implementing Java inner classes using descriptors
    Submitter: George Sakkis - gsakkis at rutgers.edu
    Last Updated: 2005/07/08
    Version no: 1.1
    Category: OOP
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/409366
    '''
    if hasattr(cls, '__outer__'):
        raise TypeError('Cannot set attribute "__outer__" in inner class')
    class InnerDescriptor(object):
        def __get__(self, outer, outercls):
            if outer is None:
                raise AttributeError('An enclosing instance that contains '
                           '%s.%s is required' % (cls.__name__, cls.__name__))
            clsdict = cls.__dict__.copy()
            # explicit read-only reference to the outer instance
            clsdict['__outer__'] = property(lambda self: outer)
            # implicit lookup in the outer instance
            clsdict['__getattr__'] = lambda self,attr: getattr(outer,attr)
            def __setattr__(this, attr, value):
                # setting an attribute in the inner instance sets the
                # respective attribute in the outer instance if and only if
                # the attribute is already defined in the outer instance
                if hasattr(outer, attr): setattr(outer,attr,value)
                else: super(this.__class__,this).__setattr__(attr,value)
            clsdict['__setattr__'] = __setattr__
            return type(cls.__name__, cls.__bases__, clsdict)
    return InnerDescriptor()


class CanvasFrame(wx.Frame):
    """A minimal wx.Frame containing a matplotlib figure"""
    def __init__(self, title='frame', size=(550, 350)):
        wx.Frame.__init__(self, None, -1, title=title, size=size)
        self.SetBackgroundColour(wx.NamedColor("WHITE"))
        self.figure = mpl.figure.Figure(figsize=(5, 4), dpi=100)
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
    """A wx.Frame for plotting a scrollable 2D grid of receptive fields, with neuron and time labels.
    rfs is a list of (nt, width, height) sized receptive fields of uint8 RGB data, one per neuron"""
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
        """This could be useful..."""
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
def pad0s(val, ndigits=2):
    """Returns a string rep of val, padded with enough leading 0s
    to give you a string rep with ndigits in it"""
    if val.__class__ != int:
        raise ValueError, '%r isn\'t an int' % val
    val = str(val)
    nzerostoadd = ndigits - len(val)
    val = '0'*nzerostoadd + val
    return val

def txtdin2binarydin(fin, fout):
    """Converts a csv text .din file to an int64 binary .din file"""
    fi = file(fin, 'r') # open the din file for reading in text mode
    fo = file(fout, 'wb') # for writing in binary mode
    for line in fi:
        line = line.split(',')
        '''
        # for old NVS display, converts from NVS condition numbers (which increment with repeats) to dimstim sweepis (which don't)
        nruns = 18
        line[1] = int(line[1]) % nruns
        '''
        fo.write( struct.pack('@qq', int(line[0]), int(line[1])) ) # write both values out as a C long longs, using the system's native ('@') byte order
    fi.close()
    fo.close()
    print 'Converted ascii din: %r to binary din: %r' % (fin, fout)

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
        fin = os.path.join(path, dinfname) + '.csv'
        fout = os.path.join(path, dinfname) + '.din'
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
                os.rename(os.path.join(path, fname), os.path.join(path, newfname))

def csv2binary(fin, multiplier=1e6):
    """Exports spike data in a csv file, with cells in the columns and times down the rows,
    into int64 binary files, one for each neuron. Takes csv values and multiplies them by
    multiplier before saving"""
    fin = os.path.normpath(fin)
    fi = file(fin, 'r') # open csv file for reading in text mode
    print 'Exporting %s to:' % fi.name
    firstline = fi.next()
    nneurons = len(firstline.split(','))
    fi.seek(0)
    data = [] # nested list, one entry per neuron
    for ni in range(nneurons):
        data.append([]) # init each neuron's list
    for line in fi:
        line = line.replace('\n', '') # strip the newline character
        line = line.split(',')
        for ni, strval in enumerate(line): # going horizontally across the line
            try:
                data[ni].append(int(round(float(strval)*multiplier)))
            except ValueError: # strval is empty string
                pass
    fi.close()
    #return data
    path = os.path.splitext(fi.name)[0] # extensionless filename
    try:
        os.mkdir(path) # make a dir with that name
    except OSError: # dir already exists
        pass
    tail = os.path.split(path)[-1] # just the extensionless filename
    for ni, neuron in enumerate(data):
        fname = os.path.join(path, tail) + '_t' + pad0s(ni, ndigits=len(str(nneurons))) + '.spk'
        fo = file(fname, 'wb') # for writing in binary mode
        print fo.name
        for spiketime in neuron: # write each spiketime to the file, there should be a more streamlined way to do this
            fo.write( struct.pack('@q', spiketime) ) # write the value out as a C long long, using the system's native ('@') byte order
        fo.close()

def warn(msg, level=2, exit_val=1):
    """Standard warning printer. Gives formatting consistency. Stolen from IPython.genutils"""
    if level > 0:
        header = ['', '', 'WARNING: ', 'ERROR: ', 'FATAL ERROR: ']
        print >> sys.stderr, '%s%s' % (header[level],msg)
        if level == 4:
            print >> sys.stderr,'Exiting.\n'
            sys.exit(exit_val)
'''
def warn(msg):
    import warnings
    warnings.warn(msg, category=RuntimeWarning, stacklevel=2)
'''
def unique(seq):
    """Return unique items from a 1-dimensional sequence. Stolen from numpy.unique().
    Dictionary setting is quite fast"""
    result = {}
    for item in seq:
        result[item] = None
    return result.keys()
'''
def unique(objlist):
    """Returns the input list minus any repeated objects it may have had. Also defined in dimstim"""
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
    """Convert to iterable. If input is iterable, returns it. Otherwise returns it in a list.
    Useful when you want to iterate over an object (like in a for loop),
    and you don't want to have to do type checking or handle exceptions
    when the object isn't a sequence"""
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

def to2d(arr):
    """Converts a 1D array to a 2D array with just a singleton row
    If arr is already 2D, just returns it. If it's anything more than 2D, raises an error"""
    nd = arr.ndim
    assert nd in [1, 2], 'array rank > 2'
    if nd == 1:
        arr = arr.reshape(1, -1)
    return arr

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

def histogram2d(x, y, bins=10, range=None, normed=False):
    """Compute the 2D histogram for a dataset (x,y) given the edges or
    the number of bins.

    Stolen from np.histogram2d() in numpy 1.0
    Modified by mspacek to allow normed='pdf' or normed='pmf' (prob mass function)

    histogram2d(x, y, bins=10, range=None, normed=False) -> H, xedges, yedges

    Compute the 2D histogram from samples x,y.

    Parameters
    ----------
    x,y: 1D data series. Both arrays must have the same length.
    bins: Number of bins -or- [nbin x, nbin y] -or-
         [bin edges] -or- [x bin edges, y bin edges].
    range:  A sequence of lower and upper bin edges (default: [min, max]).
    normed: True or False.

    The histogram array is a count of the number of samples in each
    two dimensional bin.
    Setting normed to 'pdf' returns a density rather than a bin count.
    """
    try:
        N = len(bins)
    except TypeError:
        N = 1
        bins = [bins]
    x = asarray(x)
    y = asarray(y)
    if range is None:
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
    else:
        xmin, xmax = range[0]
        ymin, ymax = range[1]
    if N == 2:
        if np.isscalar(bins[0]):
            xnbin = bins[0]
            xedges = np.linspace(xmin, xmax, xnbin+1)
        else:
            xedges = asarray(bins[0], float)
            xnbin = len(xedges)-1
        if np.isscalar(bins[1]):
            ynbin = bins[1]
            yedges = np.linspace(ymin, ymax, ynbin+1)
        else:
            yedges = asarray(bins[1], float)
            ynbin = len(yedges)-1
    elif N == 1:
        ynbin = xnbin = bins[0]
        xedges = np.linspace(xmin, xmax, xnbin+1)
        yedges = np.linspace(ymin, ymax, ynbin+1)
    else:
        yedges = asarray(bins, float)
        xedges = yedges.copy()
        ynbin = len(yedges)-1
        xnbin = len(xedges)-1

    dxedges = np.diff(xedges)
    dyedges = np.diff(yedges)

    # Flattened histogram matrix (1D)
    hist = np.zeros((xnbin)*(ynbin), int)

    # Count the number of sample in each bin (1D)
    xbin = np.digitize(x, xedges)
    ybin = np.digitize(y, yedges)

    # Values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    xdecimal = int(-np.log10(dxedges.min()))+6
    ydecimal = int(-np.log10(dyedges.min()))+6
    on_edge_x = np.where(np.around(x, xdecimal) == np.around(xedges[-1], xdecimal))[0]
    on_edge_y = np.where(np.around(y, ydecimal) == np.around(yedges[-1], ydecimal))[0]
    xbin[on_edge_x] -= 1
    ybin[on_edge_y] -= 1
    # Remove the true outliers
    outliers = (xbin==0) | (xbin==xnbin+1) | (ybin==0) | (ybin==ynbin+1)
    xbin = xbin[outliers==False] - 1
    ybin = ybin[outliers==False] - 1

    # Compute the sample indices in the flattened histogram matrix.
    if xnbin >= ynbin:
        xy = ybin*(xnbin) + xbin
    else:
        xy = xbin*(ynbin) + ybin

    # Compute the number of repetitions in xy and assign it to the flattened
    # histogram matrix.
    flatcount = np.bincount(xy)
    indices = np.arange(len(flatcount))
    hist[indices] = flatcount

    # Shape into a proper matrix
    shape = np.sort([xnbin, ynbin])
    histmat = hist.reshape(shape)
    if (shape == (ynbin, xnbin)).all():
        histmat = histmat.T

    if normed:
        if normed == 'pdf':
            diff2 = np.outer(dxedges, dyedges)
            histmat = histmat / diff2 / histmat.sum()
        elif normed == 'pmf':
            histmat = histmat / float(histmat.sum())
        else:
            raise ValueError, 'unknown normed value %s' % normed
    return histmat, xedges, yedges
'''
def histogram2dold(x, y, bins, normed=False):
    """Compute the 2D histogram for a dataset (x,y) given the edges or
    the number of bins. Stolen from np.histogram2d() (numpy 1.0b5), modified to allow
    normed='pdf' or normed='pmf' (prob mass function)

    NOTE: THIS HAS A SERIOUS BUG, IT FAILS TO TAKE THE TRANSPOSE AT ONE POINT, OR SOMETHING,
    DON'T USE!!!!!!!!

    Returns histogram, xedges, yedges.
    The histogram array is a count of the number of samples in each bin.
    The array is oriented such that H[i,j] is the number of samples falling
        into binx[j] and biny[i].
    Data falling outside of the edges are not counted.
    """
    try:
        N = len(bins)
    except TypeError:
        N = 1
        bins = [bins]
    if N == 2:
        if np.isscalar(bins[0]):
            xnbin = bins[0]
            xedges = np.linspace(x.min(), x.max(), xnbin+1)
        else:
            xedges = asarray(bins[0], float)
            xnbin = len(xedges)-1
        if np.isscalar(bins[1]):
            ynbin = bins[1]
            yedges = np.linspace(y.min(), y.max(), ynbin+1)
        else:
            yedges = asarray(bins[1], float)
            ynbin = len(yedges)-1
    elif N == 1:
        ynbin = xnbin = bins[0]
        xedges = np.linspace(x.min(), x.max(), xnbin+1)
        yedges = np.linspace(y.max(), y.min(), ynbin+1)
        xedges[-1] *= 1.0001
        yedges[-1] *= 1.0001
    else:
        yedges = asarray(bins, float)
        xedges = yedges.copy()
        ynbin = len(yedges)-1
        xnbin = len(xedges)-1

    # Flattened histogram matrix (1D)
    hist = np.zeros((xnbin)*(ynbin), int)

    # Count the number of sample in each bin (1D)
    xbin = np.digitize(x,xedges)
    ybin = np.digitize(y,yedges)

    # Remove the outliers
    outliers = (xbin==0) | (xbin==xnbin+1) | (ybin==0) | (ybin == ynbin+1)

    xbin = xbin[outliers==False]
    ybin = ybin[outliers==False]

    # Compute the sample indices in the flattened histogram matrix.
    if xnbin >= ynbin:
        xy = ybin*(xnbin) + xbin
        shift = xnbin + 1
    else:
        xy = xbin*(ynbin) + ybin
        shift = ynbin + 1

    # Compute the number of repetitions in xy and assign it to the flattened
    #  histogram matrix.
    flatcount = np.bincount(xy)
    indices = np.arange(len(flatcount)-shift)
    hist[indices] = flatcount[shift:]

    # Shape into a proper matrix
    histmat = hist.reshape(xnbin, ynbin)

    if normed == 'pdf':
        diff2 = np.outer(np.diff(yedges), np.diff(xedges))
        histmat = histmat / diff2 / histmat.sum()
    elif normed == 'pmf':
        histmat = histmat / float(histmat.sum())
    return histmat, xedges, yedges
'''
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
    """Return a string with the binary representation of an integer, or sequence of integers.
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
    """Return a string with the binary representation of an integer. If necessary, will append leading zeros
    if result is less than minbits long. Seems like np.binary_repr() is a somewhat faster alternative.
    First 2 lines stolen from Andrew Gaul <andrew@gaul.org> off the web"""
    l = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    s = ''.join(map(lambda x, l=l: l[int(x, 16)], hex(i)[2:]))
    s = s.lstrip('0') # strip s of leading zeros
    nzerostoadd = minbits - len(s)
    s = '0'*nzerostoadd + s # add enough leading zeros to get requested minbits
    return s

# an alternative would be to use int('10110', base=2) for each column, probably slower though
def binarray2int(bin):
    """Takes a 2D binary array (only 1s and 0s, with rows LSB to MSB from top to bottom)
    and returns the base 10 integer representations of the columns"""
    #assert type(bin) == type(array)
    bin = to2d(bin) # ensure it's 2D. If it's 1D, force it into having a singleton row
    nbits = bin.shape[0] # length of the first dimension, ie the number of rows
    multiplier = []
    for i in range(nbits):
        multiplier.append(2**i)
    multiplier = array(multiplier, ndmin=2).transpose() # convert from list and transpose to a column vector (have to make it 2D to transpose)
    #print multiplier
    x = bin*multiplier
    #print x
    return x.sum(axis=0) # sum over the first dimension (the rows), that way, you're left with only columns in a row vector

def getbinarytable(nbits=8):
    """Generates a 2D binary table containing all possible words for nbits, with bits in the rows and words in the columns (LSB to MSB from top to bottom)"""
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
        row = np.tile(pattern, [1, npatterns])
        x[bit::,::] = row
    return x
    '''
    # this seems to be the fastest method:
    x = []
    for bit in range(nbits): # one row at a time
        pattern = array(0, dtype=np.int8).repeat(2**bit)
        pattern = cat((pattern, array(1, dtype=np.int8).repeat(2**bit)))
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        row = np.tile(pattern, [1, npatterns])
        x.append(row)
    return cat(x)

def enlarge(a, x=2, y=None):
    """Enlarges 2D image array a using simple pixel repetition in both dimensions.
    Enlarges by factor x horizontally and factor y vertically.
    If y is left as None, uses factor x for both dimensions."""
    a = asarray(a)
    assert a.ndim == 2
    if y == None:
        y = x
    for factor in (x, y):
        assert factor.__class__ == int
        assert factor > 0
    return a.repeat(y, axis=0).repeat(x, axis=1)

def charfind(string, char):
    """Finds char in string, returns matching indices. There's gotta be a built-in way to do this somewhere..."""
    assert len(char) == 1
    i = []
    # maybe more efficient to use .find() method on successively smaller slices of string
    for si, s in enumerate(string):
        if s == char:
            i.append(si)
    return i
'''
def shuffle(x):
    """Takes an input list x and returns a shuffled (without replacement) copy. Its only benefit
    over and above random.sample() is that you don't have to pass a second argument len(x)
    every time you use it.
    In NumPy, it's better (and faster) to use np.random.shuffle()"""
    return random.sample(x, len(x))
'''
def shuffle(seq):
    """Takes a sequence and returns a shuffled (without replacement) copy.
    Its only benefit over np.random.shuffle is that it returns a copy instead of shuffling in-place"""
    result = copy(seq)
    np.random.shuffle(result) # shuffles in-place, doesn't convert to an array
    return result
'''
def randomize(x):
    """Takes an input list x and returns a randomized (with replacement) output list of
    the same length, sampled from the input sequence"""
    y = [] # init output list
    for i in range(0, len(x)):
        y.append(random.choice(x))
    return y
'''
def randomize(seq):
    """Returns a randomized (with replacement) output sequence sampled from
    (and of the same length as) the input sequence"""
    n = len(seq)
    i = np.random.randint(n, size=n) # returns random ints from 0 to len(seq)-1
    if seq.__class__ == np.ndarray:
        return np.asarray(seq)[i] # use i as random indices into seq, return as an array
    else:
        return list(np.asarray(seq)[i]) # return as a list

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

def combgen(objects, r=2, i=None, level=0):
    """Generator that yields, without replacement, all length r possible combinations of objects from a length n sequence.
    Eg, if objects=[0,1,2] and r=2, this yields [0,1], [0,2], and [1,2], one at a time.
    A recursive generator is used in order to create the necessary r number of nested for loops.
    This is cool (my first generator!), but deep recursion is slow"""
    objects = asarray(objects)
    assert r <= len(objects)
    try: # recursive case
        if i == None:
            i = [0]*r # stores all the current index values for all r nested for loops
        if level == 0: # handles special case for starting index of top level for loop
            starti = 0
        else:
            starti = i[level-1] + 1 # start this level's loop index at one greater than the previous level's current loop index
        for i[level] in range(starti, len(objects)+1): # not too sure why this is n+1, but it works
            for comb in combgen(objects, r=r, i=i, level=level+1): # iterate over next level's generator
                yield comb # yield whatever the next level (level+1) yields, pass it on up to the previous level (level-1)
    except IndexError: # base case, we're at the deepest recursion level (innermost for loop). IndexError comes from i[level] being out of range
        #if len(i) == 1:
        #    yield objects[i[0]] # no need to yield them in a list
        #else:
            yield objects[i] # use the current index state for all levels to yield a combination of objects

def combs(objects, r=2):
    """Returns all nCr possible combinations of items in objects, in a 1D array of arrays.
    Generates code with the right number of nested for loops, faster than combgen()"""
    objects = asarray(objects)
    dtype = objects.dtype
    n = len(objects)
    assert r <= n
    i = asarray([0]*r)
    combs = np.empty(nCr(n, r), dtype=np.object) # stores all combinations, will be a 1D array of arrays
    combi = -1

    code = ''
    tabs = ''
    code += tabs+'for i[0] in range(0, n):\n' # this is the outermost for loop
    tabs += '\t'
    for level in range(1, r): # here come the inner nested for loops...
        code += tabs+'for i['+str(level)+'] in range(i['+str(level-1)+']+1, n):\n'
        tabs += '\t'

    # here's the innermost part of the nested for loops
    code += tabs + 'combi += 1\n'
    code += tabs + 'combs[combi] = objects[i]\n'
    #print code

    exec(code) # run the generated code
    return combs
    '''
    # example of what the generated code looks like for r==3:
    for i[0] in range(0, n):
        for i[1] in range(i[0]+1, n):
            for i[2] in range(i[1]+1, n):
                combi += 1
                combs[combi] = objects[i]
    '''

def argcombs(objects, r=2):
    """Returns all nCr possible combinations of indices into objects.
    You'd think this would be faster than combs(), but it doesn't seem to be"""
    n = len(objects)
    assert n < 2**8 # this way, we can use uint8's instead of int32's to save memory
    assert r <= n
    i = asarray([0]*r)
    argcombs = np.zeros((nCr(n, r), r), dtype=np.uint8)
    combi = -1

    code = ''
    tabs = ''
    code += tabs+'for i[0] in range(0, n):\n' # this is the outermost for loop
    tabs += '\t'
    for level in range(1, r): # here come the inner nested for loops...
        code += tabs+'for i['+str(level)+'] in range(i['+str(level-1)+']+1, n):\n'
        tabs += '\t'

    # here's the innermost part of the nested for loops
    code += tabs + 'combi += 1\n'
    code += tabs + 'argcombs[combi, :] = i\n'
    #print code

    exec(code) # run the generated code
    return argcombs
    '''
    # example of what the generated code looks like for r==3:
    for i[0] in range(0, n):
        for i[1] in range(i[0]+1, n):
            for i[2] in range(i[1]+1, n):
                combi += 1
                argcombs[combi, :] = i
    '''

def nCrsamples(objects, r, nsamples=None):
    """Returns a list of nsamples unique samples, each of length r, sampled from objects"""
    maxnsamples = nCr(len(objects), r)
    if nsamples == None:
        nsamples = maxnsamples # return all possible combinations
    if nsamples > maxnsamples: # make sure we're not being asked for more than the maximum possible number of unique samples
        raise ValueError, 'requested unique nsamples (%d) is larger than len(objects) choose r (%d C %d == %d)' % (nsamples, len(objects), r, maxnsamples)
    # I've set the criteria for generating a table to be never, cuz generating the table and then sampling it almost always takes longer (at least for maxnsamples as high as 325, say) than just picking combs at random and making sure they're unique
    if maxnsamples < 0: # generate a table of all possible combinations, and then just pick nsamples from it without replacement
        table = combs(objects, r)
        samples = random.sample(table, nsamples)
    elif r == 1: # we're just choosing one item from objects at a time
        samples = random.sample(objects, nsamples)
    else: # the number of possible combs is inconveniently large to completely tabulate, pick some combinations at random and make sure each comb is unique
        samples = []
        samplei = 0
        while samplei < nsamples:
            sample = random.sample(objects, r) # choose r objects at random
            sample.sort() # sort for sake of comparison with other samples, important cuz this removes any differences due to permuatations (as opposed to combs)
            if sample not in samples: # make sure they're not the same set of objects as any previous set in samples
                samples.append(sample) # add it to the list of samples
                samplei += 1
    return samples

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
        self.thousandsSep = '' # default to not using a thousands separator

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

    def pprint_val(self, x):
        xp = (x-self.offset)/10**self.orderOfMagnitude
        if np.absolute(xp) < 1e-8: xp = 0
        s = self.format % xp
        if self.thousandsSep: # add thousands-separating characters
            if s.count('.'): # it's got a decimal in there
                s = re.sub(r'(?<=\d)(?=(\d\d\d)+\.)', self.thousandsSep, s) # use the regexp for floats
            else: # it's an int
                s = re.sub(r'(?<=\d)(?=(\d\d\d)+$)', self.thousandsSep, s) # use the regexp for ints
        return s

class neuropyAutoLocator(mpl.ticker.MaxNLocator):
    """A tick autolocator that generates more ticks than the standard mpl autolocator"""
    def __init__(self):
        #mpl.ticker.MaxNLocator.__init__(self, nbins=9, steps=[1, 2, 5, 10]) # standard autolocator
        mpl.ticker.MaxNLocator.__init__(self) # use MaxNLocator's defaults instead

def normalize(seq):
    """Normalizes a sequence, returning zeros if sum(seq) == 0"""
    a = asarray(seq)
    if a.sum() == 0: # numpy doesn't raise ZeroDivisionErrors for some reason
        return zeros(a.shape) # just return zeros
    else:
        return a / float(a.sum()) # return it normalized

def ensurenormed(p, atol=1e-8):
    """Ensures p is normalized. Returns p unchanged if it's already normalized,
    otherwise, prints a warning and returns it normalized. atol is how close to 1.0
    p.sum() needs to be"""
    p = asarray(p)
    psum = p.sum()
    if not approx(psum, 1.0, atol=atol): # make sure the probs sum to 1
        print 'ps don''t sum to 1, they sum to %f instead, normalizing for you' % psum
        p /= float(psum)
    return p

def logy(x, base=10):
    """Performs log of x with specified base"""
    return log(x) / log(base)

def log_no_sing(x, subval=0.0, base=np.e):
    """Performs log on array x, ignoring any zeros in x to avoid singularities,
    and returning subval in their place in the result"""
    x = asarray(x)
    singi = x==0 # find the singularities
    x[singi] = 1 # replace 'em with 1s, or anything else that's safe to take the log of
    result = logy(x, base=base) # now it's safe to take the log
    result[singi] = subval # substitute the result where the singularities were with the substitution value
    return result

def log10_no_sing(x, subval=0.0):
    """Performs log10 on x, ignoring singularities"""
    return log_no_sing(x, subval=subval, base=10)

def log2_no_sing(x, subval=0.0):
    """Performs log2 on x, ignoring singularities"""
    return log_no_sing(x, subval=subval, base=2)

def entropy(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob values in p"""
    p = ensurenormed(p)
    return -(p * log2(p)).sum()

def entropy_no_sing(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob values in p
    Ignore singularities in p (assumes their contribution to entropy is zero)"""
    p = ensurenormed(p)
    return -(p * log2_no_sing(p, subval=0.0)).sum()


def mutualinfo(XY):
    """Given the joint PDF of two variables, returns the mutual information (in bits) between the two
    I = sum_X sum_Y P(x, y) * log2( P(x, y) / (P(x) * P(y)) )
    where P(x) and P(y) are the marginal distributions taken from the joint"""
    XY = asarray(XY)
    assert XY.ndim == 2
    XY = ensurenormed(XY)
    # calculate the marginal probability distributions for X and Y from the joint
    X = XY.sum(axis=1) # sum over the rows of the joint, get a vector nrows long
    Y = XY.sum(axis=0) # sum over the cols of the joint, get a vector ncols long
    I = 0.0
    for xi, x in enumerate(X):
        for yi, y in enumerate(Y):
            if XY[xi, yi] == 0 or (x * y) == 0: # avoid singularities
                pass # just skip it, assume info contributed is 0 (?????????????????)
            else:
                I += XY[xi, yi] * log2( XY[xi, yi] / (x * y) )
    return I

def DKL(p, q):
    """Kullback-Leibler divergence from true probability distribution p to arbitrary distribution q"""
    assert len(p) == len(q)
    p = ensurenormed(p)
    q = ensurenormed(q)
    return sum([ pi * log2(pi/float(qi)) for pi, qi in zip(p, q) if pi != 0 and qi != 0 ] ) # avoid singularities

def DJS(p, q):
    """Jensen-Shannon divergence, a symmetric measure of divergence between distributions p and q"""
    p = asarray(p) # required for adding p and q
    q = asarray(q)
    m = 1 / 2.0 * (p + q)
    return 1 / 2.0 * ( DKL(p, m) + DKL(q, m) )


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
        print '%d iters,' % self.model.iters,
        #print 'hi:', self.hi.__repr__()
        #print 'Jij:', self.Jij.__repr__()

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
