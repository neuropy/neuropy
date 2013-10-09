"""Global variables that can be modified by the user at the IPython command line.
Access programatically using:

get_ipython().user_ns['VARNAME']
"""
# __future__ imports don't seem to work when executing this file in IPython, see main.py:
#from __future__ import division
#from __future__ import print_function

import os
from core import mergeuniquedictvals

DATAPATH = os.path.expanduser('~/data')
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

"""Spike correlation time range windows"""
SCWIDTH = 10 # sec
SCTRES = 1 # sec
SCLIMITS = -0.01, 0.13

"""Multiunit activity time range windows"""
MUAWIDTH = 0.25 # sec
MUATRES = 0.1 # sec

"""MUA synchrony index time range windows"""
MUASIWIDTH = 10 # sec
MUASITRES = 1 # sec

"""LFP synchrony index time range windows"""
LFPSIWIDTH = 32.768 # sec (2**15 ms)
LFPSITRES = 1 # sec
LFPSILOWBAND = 0.5, 7 # Hz
LFPSIHIGHBAND = 7, 100 # Hz

"""List of sorted track IDs"""
TRACKS = ['ptc15.tr7c', 'ptc22.tr1', 'ptc22.tr2']

"""Track-specific superficial, middle and deep layer ranges (um), inferred from
track.pospdf and sc.pos"""
LAYERS = {'ptc15.tr7c': [(0, 900), (900, 1100), (1100, 2000)],
          'ptc22.tr1':  [(0, 500), (500,  700), ( 700, 2000)],
          'ptc22.tr2':  [(0, 550), (550,  700), ( 700, 2000)]}
          
"""IDs of blankscreen recordings"""
BLANKRIDS = {'ptc15.tr7c': ['87'],
             'ptc22.tr1':  ['07', '09', '11', '21'],
             'ptc22.tr2':  ['27', '32', '36']}
             
"""IDs of msequence recordings"""
MSEQRIDS = {'ptc15.tr7c': ['70', '81', '91', '92', '94'],
            'ptc22.tr1':  ['04', '17'],
            'ptc22.tr2':  ['26', '34']}
            
"""IDs of movie recordings"""
MOVRIDS = {'ptc15.tr7c': ['76', '96'],
           'ptc22.tr1':  ['05', '06', '08', '10', '19', '20'],
           'ptc22.tr2':  ['28', '33']}
           
"""IDs of drift bar recordings"""
DRIFTRIDS = {'ptc15.tr7c': ['71', '82'],
             'ptc22.tr1':  ['03', '18'],
             'ptc22.tr2':  ['25', '31']}

"""IDs of full field flash recordings"""
FFFRIDS = {'ptc15.tr7c': ['69', '78', '79', '88', '93'],
           'ptc22.tr1':  ['01', '02', '12', '16', '22'],
           'ptc22.tr2':  ['23', '24', '35']}

"""IDs of flash grating recordings"""
FGRIDS = {'ptc15.tr7c': ['73'],
          'ptc22.tr1':  ['13'],
          'ptc22.tr2':  ['30']}

"""Per-track list of relevant recordings, merged from those above"""
RIDS = mergeuniquedictvals([BLANKRIDS, MSEQRIDS, MOVRIDS, DRIFTRIDS])
"""Per-track list of blankscreen and msequence recordings"""
BLANKMSEQRIDS = mergeuniquedictvals([BLANKRIDS, MSEQRIDS])
"""Per-track list of movie and driftbar recordings"""
MOVDRIFTRIDS = mergeuniquedictvals([MOVRIDS, DRIFTRIDS])

NULLDIN = 65535 # integer value in stimulus .din files used as NULL (stimulus off)
