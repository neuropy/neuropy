"""Code executed in user namespace on startup.

ip = get_ipython() # get current interactive shell
ip.ex() # execute strings in user namespace
ip.ev() # evaluates string
ip.write() # write a string to the shell
ip.push() # push a variable in a dict to user namespace

"""

from __future__ import division
import os

import pylab as pl

from animal import Animal
from track import Track
from recording import Recording


def cf():
    import pylab
    pylab.close('all')

def fontsize(pts=None):
    """Get/set font size in points of various typical plot features. Needs work"""
    rc = get_ipython().user_ns['rcParams']
    if pts == None:
        return rc['font.size']
    rc['axes.titlesize'] = pts
    rc['axes.labelsize'] = pts
    rc['xtick.labelsize'] = pts
    rc['ytick.labelsize'] = pts
    rc['legend.fontsize'] = pts
    rc['font.size'] = pts
    print('font size set to %r points' % pts)

ip = get_ipython()
ip.pdb = True # not sure what the difference is between .pdb and .call_pdb
ip.call_pdb = True
ip.magic("pylab")
ip.Completer.greedy = True
ip.Completer.omit__names = 1 # all 'magic' names (``__foo__``) will be excluded
rc = get_ipython().user_ns['rcParams']
rc['mathtext.default'] = 'regular' # use same font for math mode as regular text mode
