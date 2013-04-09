"""Global variables that can be modified by the user at the IPython command line.
Access programatically using:

get_ipython().user_ns['VARNAME']
"""

import os

MOVIEPATH = os.path.expanduser('~/data/mov')

# for each recording, load all Sorts, or just the most recent one?
LOADALLSORTS = False

"""Mean spike rate that delineates normal vs "quiet" neurons. 0.1 Hz seems reasonable if you
plot mean spike rate distributions for all the neurons in a given track. But, if you want a
reasonable looking DJS histogram withouth a lot of missing netstates, you need to exclude
more low firing rate cells, 0.5 works better"""
MINRATE = 0.05 # Hz
"""Calculate a TrackNeuron's meanrate according to its trange (period between its first and
last spike), or according to its track's entire duration. Need to reload the track or call
Track.calc_meanrates() after changing this on the fly"""
TRACKNEURONPERIOD = 'track' # 'trange' or 'track'
# ditto for recordings:
RECNEURONPERIOD = 'recording' # 'trange' or 'recording'

"""NeuronCode (Ising matrix) and network state parameters"""
CODEKIND = 'binary'
# values to use for CODEKIND codes, doesn't seem to make any difference to correlation
# calcs, unless set to really extreme values like [-100s, 100s], which is probably due to
# int8 overflow
CODEVALS = [0, 1]
CODETRES = 20000 # us
CODEPHASE = 0 # deg
CODEWORDLEN = 10 # in bits

"""LFP synchrony index time range windows"""
SIWIDTH = 32.768 # sec (2**15 ms)
SITRES = 1 # sec
SILOWBAND = 0.5, 7 # Hz
SIHIGHBAND = 7, 100 # Hz

"""Code correlation time range windows"""
CCWIDTH = 16.384 # sec
CCTRES = 1 # sec

"""Multiunit activity time range windows"""
MUAWIDTH = CCWIDTH # sec
MUATRES = CCTRES # sec
# time resolution for alternative method of convolving with smoothing window
MUASMOOTHRES = 1 # sec

"""Distance ranges down vertical axis of probe designating superficial and deep cells"""
SUPRANGE = 0, 500 # um
DEEPRANGE = 700, 2000 # um
