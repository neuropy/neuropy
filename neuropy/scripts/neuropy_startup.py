"""Runs on PyShell/PyCrust startup, as well as IPython startup"""

try:
    shell.run('from __future__ import division')
except: # we're running in some environment where shell isn't defined, like ipython
    pass


from minimal_startup import * # I keep this in ~/scripts
import __main__ # has _, needs to be imported explicitly

import time
import types
import struct
import re
import StringIO
import random
import math

import matplotlib as mpl
mpl.use('WXAgg')
mpl.interactive(True)

import scipy as sp
import scipy.signal as sig
import scipy.weave as weave

from numpy.core.ma import array as mar
import wx


#shell.redirectStdout(redirect=True) # does this actually do anything in PyShell?

#import deep_reload # modifies __builtin__.reload() to do a deep reload
#import LazyPython
#sys.excepthook = LazyPython.LazyPython() # this doesn't seem to work

cd('c:/home/mspacek/Desktop')

try:
    clr = shell.clear
except NameError: # we're running in some environment where shell isn't defined, like ipython
    pass

def r():
    """Refreshes neuropy.
    WARNING: Seems to cause problems in IPython"""
    refresh('neuropy')

def cf():
    """Closes all figures. This needs to be extended to include all neuropy generated wx.Frames,
    not just matplotlib figures"""
    pl.close('all')
    print 'all figures closed'

_original = __main__.__dict__.keys()

print 'importing neuropy'
import neuropy
from neuropy import *
from neuropy import _data # ensure it's imported, in spite of leading _, useful for user examination of default Data object

import pylab as pl
from pylab import figure, plot, loglog, hist, bar, barh, xlabel, ylabel, xlim, ylim, title, gcf, gca, get_current_fig_manager as gcfm, axes, axis, hold, imshow
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
