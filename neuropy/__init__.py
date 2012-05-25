r"""Experimental and model neuronal and stimulus data analysis in Python

          object hierarchy:

                                level:

                Animal            1
                  |
                Track             2
                  |
              Recording           3
             /         \
       Experiment      Sort       4
            |           |
      (Cat15Movie)    Neuron      5

As for variable scope, this is how things work:

Main Qt event loop creates a
    NeuropyWindow which runs an
        IPython session which imports
            neuropy modules
    
So, something like get_ipython() only works from within those neuropy modules that have been
imported during the IPython session, but not in any of NeuropyWindow's methods.
"""

__authors__ = ["Martin Spacek"]
__version__ = 0.2
'''
from animal import Animal
from track import Track
from recording import Recording
from experiment import Experiment
from sort import Sort
from neuron import Neuron
'''
#from test import test

# init and load some neuropy objects:
'''
print 'Initing and loading Track(\'7c\'):'
t = Track('7c')
t.load()'''
'''
print 'Initing and loading Recording(71):'
r71 = Recording(71)
r71.load()
''''''
print 'Initing and loading Recording(75):'
r75 = Recording(75)
r75.load()
''''''
print 'Initing and loading Recording(76):'
r76 = Recording(76)
r76.load()
''''''
print 'Initing and loading Recording(92):'
r92 = Recording(92)
r92.load()
'''
