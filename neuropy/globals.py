"""Global variables that can be modified by the user at the IPython command line.
Access programatically using:

get_ipython().user_ns['VARNAME']
"""
# __future__ imports don't seem to work when executing this file in IPython, see main.py:
#from __future__ import division
#from __future__ import print_function

import os

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

"""LFP synchrony index time range windows"""
SIWIDTH = 32.768 # sec (2**15 ms)
SITRES = 1 # sec
SILOWBAND = 0.5, 7 # Hz
SIHIGHBAND = 7, 100 # Hz

"""Code correlation time range windows"""
SCWIDTH = 16.384 # sec
SCTRES = 1 # sec

"""Multiunit activity time range windows"""
MUAWIDTH = SCWIDTH # sec
MUATRES = SCTRES # sec

"""Distance ranges down vertical axis of probe designating superficial, middle and
deep cells"""
SUPRANGE = 0, 500 # um
MIDRANGE = 500, 700 # um
DEEPRANGE = 700, 2000 # um

"""Track-specific superficial, middle and deep layer ranges (um), inferred from
track.pospdf and sc.pos"""
LAYERS = {'ptc15.tr7c': [(0, 900), (900, 1100), (1100, 2000)],
          'ptc22.tr1':  [(0, 500), (500,  700), ( 700, 2000)],
          'ptc22.tr2':  [(0, 550), (550,  700), ( 700, 2000)],
         }
"""IDs of blankscreen recordings"""
BLANKRIDS = {'ptc15.tr7c': ['87'],
             'ptc22.tr1': ['07', '09', '11', '21'],
             'ptc22.tr2': ['27', '32', '36'],
            }
"""IDs of msequence recordings"""
MSEQRIDS = {'ptc15.tr7c': ['70', '81', '91', '92', '94'],
            'ptc22.tr1': ['04', '17'],
            'ptc22.tr2': ['26', '34'],
           }
"""IDs of movie recordings"""
MOVRIDS = {'ptc15.tr7c': ['76', '96'],
           'ptc22.tr1': ['05', '06', '08', '10', '19', '20'],
           'ptc22.tr2': ['28', '33'],
          }
"""IDs of drift bar recordings"""
DRIFTRIDS = {'ptc15.tr7c': ['71', '82'],
             'ptc22.tr1': ['03', '18'],
             'ptc22.tr2': ['25', '31'],
            }

NULLDIN = 65535 # integer value in stimulus .din files used as NULL (stimulus off)
