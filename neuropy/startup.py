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
    pl.close('all')
    
ip = get_ipython()
ip.call_pdb = True
ip.magic("pylab")
ip.Completer.greedy = True
ip.Completer.omit__names = 1 # all 'magic' names (``__foo__``) will be excluded
