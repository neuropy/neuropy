r"""Experimental and model neuronal and stimulus data analysis in Python

          object hierarchy:

                                level:

                Animal            1
                  |
                Track             2             Track
                  |                               |
              Recording           3           TrackSort
             /         \                          |
       Experiment      Sort       4          TrackNeuron
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
__version__ = '0.3'
