r"""Tools for analyzing experimental and model neuronal and stimulus data in python

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

__author__ = "Martin Spacek"

'''
TODO:

- Nah!: Rips should really have ids to make them easier to reference to: r[83].rip[0] instead of r[83].rip['conservative spikes'] - this means adding id prefixes to rip folder names (or maybe suffixes: 'conservative spikes.0.rip', 'liberal spikes.1.rip', etc...). Prefixes would be better cuz they'd force sorting by id in explorer (which uses alphabetical order) - ids should be 0-based of course
- worry about conversion of ids to strings: some may be only 1 digit and may have a leading zero!
- maybe make two load() f'ns for Experiment and Neuron: one from files, and a future one from a database
- make a save() f'n that pickles the object (including any of its results, like its STA, tuning curve points, etc)? - just use IPython's %store
'''

print 'importing neuropy'

from Core import *
from Recording import Recording
from Run import Run
from Experiment import Experiment
from Rip import Rip
from Movie import Movie, MSEQ32, MSEQ16
from Neuron import Neuron

from Test import test

# init and load some neuropy objects:
'''
print 'Initing and loading Track(\'7c\'):'
t = Track('7c')
t.load
'''
print 'Initing and loading Recording(71):'
r71 = Recording(71)
r71.load()
print 'Initing and loading Recording(92):'
r92 = Recording(92)
r92.load()
