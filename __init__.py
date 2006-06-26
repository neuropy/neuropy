r"""Tools for analyzing neuronal and stimulus data in python

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
print 'importing neuropy'

from Core import *
from Recording import Recording
from Run import Run
from Experiment import Experiment
from Rip import Rip
from Movie import Movie
from Neuron import Neuron

from Test import test

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



'''
# set neuropy's __all__ attrib for designating which names will be imported when
# you go 'from neuropy import *'
__importall__ = [Core] # list of mods imported above as "from mod import *"
__all__ = ['Neuron', 'Experiment', 'test']
import types
for mod in __importall__:
    for key, value in mod.__dict__.iteritems():
        # don't add modules or special names designated with a leading underscore
        if type(value) is not types.ModuleType and not key.startswith('_'):
            #print key
            __all__.append(key)
'''
