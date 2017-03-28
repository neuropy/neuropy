[neuropy](http://neuropy.github.io) is a Python application for interactive analysis
of neuronal spike, LFP, and stimulus data.

The user-facing data object hierarchy looks like this:
```
                                level:

                Animal            1
                  |
                Track             2             Track
                  |                               |
              Recording           3           TrackSort
             /    |    \                          |
          Sort   LFP  Experiment  4          TrackNeuron
           |
         Neuron                   5
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

Some more details are available in Appendix C.3 of my
[thesis](http://mspacek.github.io/mspacek_thesis.pdf).

Dependencies:

neuropy requires recent versions of the following to be installed:

* [Python](http://python.org) (3.x, 2.x no longer supported)
* [IPython](http://ipython.org) 5.3.0
* [numpy](http://numpy.org)
* [scipy](http://scipy.org)
* [matplotlib](http://matplotlib.org)
* [PyQt4](http://www.riverbankcomputing.co.uk/software/pyqt)
  ([PySide](http://pyside.org) will probably work too, but is untested)
* [Cython](http://cython.org)

Optional:

* [scikit-learn](http://scikit-learn.org)

To ensure that the [Qt v2
API](http://ipython.org/ipython-doc/dev/interactive/reference.html#pyqt-and-pyside) is used by
IPython, run `export QT_API=pyqt` at the command line before running neuropy, or add it to
your `.bashrc` or `.bash_aliases` file.

Some analyses assume the use of custom matplotlib default settings in the included
[matplotlibrc](matplotlibrc) file.

neuropy is developed in Xubuntu 14.04. It should work in other Linux distributions. In
principle, it should also work in Windows and OSX.

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
