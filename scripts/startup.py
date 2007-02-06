"""Runs on PyShell/PyCrust startup, as well as IPython startup"""

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
from numpy import arange, array, array as ar, asarray, asarray as aar, log, log2, log10, sqrt, zeros, ones, diff, concatenate, concatenate as cat, mean, median, std
from numpy.core.ma import array as mar
from pylab import figure, plot, loglog, hist, bar, barh, xlabel, ylabel, xlim, ylim, title, gcf, gca, get_current_fig_manager as gcfm, axes, axis, hold, imshow
import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

mpl.use('WXAgg')
mpl.interactive(True)

# set some numpy options
np.set_printoptions(precision=8)
np.set_printoptions(threshold=1000)
np.set_printoptions(edgeitems=5)
np.set_printoptions(linewidth=150)
np.set_printoptions(suppress=True)

#shell.redirectStdout(redirect=True) # does this actually do anything in PyShell?

#import deep_reload # modifies __builtin__.reload() to do a deep reload
#import LazyPython
#sys.excepthook = LazyPython.LazyPython() # this doesn't seem to work

def who():
    """Print all names in namespace"""
    import __main__
    print __main__.__dict__.keys()

def whos():
    """Print all names in namespace, with details"""
    exec 'pprint(locals())' in globals()

def cd(path):
    """Change directories"""
    path = path.replace('~', 'C:/home/mspacek') # make '~' a shortcut to my home
    if path == '..': # go down one directory
        path = '\\'.join(os.getcwd().split('\\')[0:-1])
    try:
        os.chdir(os.getcwd() + path) # path is relative?
    except OSError:
        os.chdir(path) # nope, path is absolute

cd('c:/home/mspacek/Desktop')

def pwd():
    """Returns working directory"""
    return os.getcwd()

def ls():
    """Returns directory contents in a list"""
    print pwd()
    return os.listdir(os.getcwd())

def ll():
    """Returns long-list directory contents in a list"""
    print pwd()
    pprint(os.listdir(os.getcwd()))
    #buf = StringIO.StringIO()
    #pprint(os.listdir(os.getcwd()), stream=buf)
    #return buf.getvalue()

def src(obj):
    """Print object's source code"""
    try:
        import inspect
        source = inspect.getsource(obj)
    except TypeError: # probalby a builtin
        print obj
        return
    except IOError: # probably entered interactively, no source file
        print obj
        return
    print inspect.getfile(obj) + ':'
    print
    print source

try:
    clr = shell.clear
except NameError: # we're running in some environment where shell isn't defined, like ipython
    pass

def refresh(modname):
    """Deletes all modules with 'modname' in their file path (mod.__file__).
    Also, deletes any objects that depended on modname. Then re-imports 'modname' as a module.
    'modname' need not have been previously imported.
     WARNING: Seems to cause problems in IPython"""
    import __main__
    md = __main__.__dict__

    print 'refreshing %s' % modname
    for key, mod in sys.modules.items():
        try:
            if mod.__file__.count(modname):
                #print 'deleting', mod
                del sys.modules[key]
        except AttributeError: # some modules don't have a .__file__ attrib
            pass

    for key in md.keys(): # for all names in the namespace
        if repr(md[key]).count(modname): # if modname shows up in this object's repr, ie if this object depends on modname
            del md[key] # delete the object
            #print 'deleted object:', key

    __import__(modname, globals(), locals(), []) # this is equivalent to "import modname", yet accepts a string for modname
    print '%s refreshed' % modname

def r():
    """Refreshes neuropy.
    WARNING: Seems to cause problems in IPython"""
    refresh('neuropy')

def cf():
    """Closes all figures. This needs to be extended to include all neuropy generated wx.Frames,
    not just matplotlib figures"""
    pl.close('all')
    print 'all figures closed'

def c():
    """Clears all names added to the namespace after the '_original' point.
    Don't try running this in IPython!"""
    import __main__
    md = __main__.__dict__
    for key in md.keys():
        if key not in _original and key != '_original':
            del md[key]
            #print 'deleted md key:', key
    print 'namespace cleared'

_original = __main__.__dict__.keys()

print 'importing neuropy'
import neuropy
from neuropy import *
from neuropy import _data # ensure it's imported, in spite of leading _, useful for user examination of default Data object
