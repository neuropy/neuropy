[neuropy](http://neuropy.github.io) is a Python application for interactive analysis
of sorted neuronal spike data and raw LFP waveforms.

The user-facing data object hierarchy looks like this:
```
                                level:

                Animal            1
                  |
                Track             2             Track
                  |                               |
              Recording           3           TrackSort
             /         \                          |
       Experiment      Sort       4          TrackNeuron
                        |
                      Neuron      5
```
Once data are loaded, they can be accessed in various ways:
```
animal.track.recording
animal.tr[trackid].r[recordingid]
recording.n[neuronid]
```
For example, to plot a population raster plot of all active cells in
recording 5 of track 1 of animal example1, you can type:
```
example1.tr1.r05.praster()
```
Different object types have different analysis methods associated with them.
Thanks to [IPython](http://ipython.org), tab completion should help you discover
analysis methods, and typing `?` after the method name or simply typing the opening
parenthesis should give you a docstring.

As for internal variable scope, this is how things work:
```
Main Qt event loop creates a
    NeuropyWindow which runs an
        IPython session which imports
            neuropy modules
```
So, something like `get_ipython()` only works from within those neuropy modules
that have been imported during the IPython session, but not in any of
`NeuropyWindow`'s methods.

Dependencies:

neuropy requires recent versions of the following to be installed:

* [Python](http://python.org) (2.7.x, 3.x hasn't been tested)
* [IPython](http://ipython.org) 1.0.dev
* [numpy](http://numpy.org)
* [scipy](http://scipy.org)
* [matplotlib](http://matplotlib.org)
* [PyQt4](http://www.riverbankcomputing.co.uk/software/pyqt)
  ([PySide](http://pyside.org) will probably work too, but is untested)
* [Cython](http://cython.org)

Optional:

* [scikit-learn](http://scikit-learn.org)

neuropy is developed in Xubuntu 12.10. It should work in other Linux distributions.
In principle, it should also work in Windows and OSX.

An older version is described in the paper [Python for large-scale electrophysiology]
(http://www.frontiersin.org/Neuroinformatics/10.3389/neuro.11.009.2008/abstract).

To run:
```
$ python main.py # in the neuropy folder
```
To install for use as a library:
```
$ python setup.py install
```
