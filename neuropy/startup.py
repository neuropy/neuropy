"""Code executed in user namespace on startup.

ip = get_ipython()
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
