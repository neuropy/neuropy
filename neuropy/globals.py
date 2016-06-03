"""Global variables that can be modified by the user at the IPython command line.
Access programatically using:

get_ipython().user_ns['VARNAME']
"""
# __future__ imports don't seem to work when executing this file in IPython, see main.py:
#from __future__ import division
#from __future__ import print_function

import os
from core import mergeuniquedictvals, dictattr

DATAROOTPATH = os.path.expanduser('~/data')
LABPATHNAME = 'slab'
#LABPATHNAME = 'blab'
DATAPATH = os.path.join(DATAROOTPATH, LABPATHNAME)
MOVIEPATH = os.path.join(DATAPATH, 'mov')
MOVIES = dictattr()

# for each recording, load all Sorts, or just the most recent one?
LOADALLSORTS = False

"""Mean spike rate that delineates normal vs "quiet" neurons. 0.1 Hz seems reasonable if you
plot mean spike rate distributions for all the neurons in a given track. But, if you want a
reasonable looking DJS histogram without a lot of missing netstates, you need to exclude
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
TMUAWIDTH = 0.02 # sec
TMUATRES = 0.005 # sec

"""MUA synchrony index time range windows"""
MUASIWIDTH = 10 # sec
MUASITRES = 1 # sec
MUASIKIND = 'nstdmed'

"""LFP power time range windows"""
LFPWIDTH = 2 # sec
LFPTRES = 0.5 # sec

"""LFP synchrony index time range windows"""
LFPSIWIDTH = 10 # sec
LFPSITRES = 2 # sec
LFPSILOWBAND = 0.5, 7 # Hz
LFPSIHIGHBAND = 15, 100 # Hz
LFPSIKIND = 'L/(L+H)'#'n3stdmed'

"""List of sorted track IDs"""
TRACKS = ['ptc15.tr7c', 'ptc22.tr1', 'ptc22.tr2']

"""Track-specific superficial, middle and deep layer ranges (um), inferred from
track.pospdf and sc.pos"""
LAYERS = {'ptc15.tr7c': [(0, 900), (900, 1100), (1100, 2000)],
          'ptc22.tr1':  [(0, 500), (500,  700), ( 700, 2000)],
          'ptc22.tr2':  [(0, 550), (550,  700), ( 700, 2000)]}

"""Polytrode type to shank width (um) mapping, from NeuroNexus 2008 catalog"""
PTSHANKWIDTHS = {'uMap54_1a':207, 'uMap54_1b':210, 'uMap54_1c':208,
                 'uMap54_2a':200, 'uMap54_2b':207}

"""Polytrode type to shank tip length (um) mapping, from NeuroNexus 2008 catalog"""
PTTIPLENGTHS = {'uMap54_1a':325, 'uMap54_1b':324, 'uMap54_1c':324,
                'uMap54_2a':324, 'uMap54_2b':324}

"""IDs of blankscreen recordings"""
BSRIDS = {'ptc15.tr7c': ['87'],
          'ptc22.tr1':  ['07', '09', '11', '21'],
          'ptc22.tr2':  ['27', '32', '36']}

"""IDs of msequence recordings"""
MSRIDS = {'ptc15.tr7c': ['70', '81', '91', '92', '94'],
          'ptc22.tr1':  ['04', '17'],
          'ptc22.tr2':  ['26', '34']}

"""IDs of natural scene movie recordings"""
NSRIDS = {'ptc15.tr7c': ['76', '96'],
          'ptc22.tr1':  ['05', '06', '08', '10', '19', '20'],
          'ptc22.tr2':  ['28', '33']}

"""IDs of drift bar recordings"""
DBRIDS = {'ptc15.tr7c': ['71', '82'],
          'ptc22.tr1':  ['03', '18'],
          'ptc22.tr2':  ['25', '31']}

"""IDs of drift grating recordings"""
DGRIDS = {'ptc15.tr7c': ['85'],
          'ptc22.tr1':  ['14']}

"""IDs of flash grating recordings"""
FGRIDS = {'ptc15.tr7c': ['73'],
          'ptc22.tr1':  ['13'],
          'ptc22.tr2':  ['30']}

"""IDs of full field flash recordings"""
FFRIDS = {'ptc15.tr7c': ['69', '78', '79', '88', '93'],
          'ptc22.tr1':  ['01', '02', '12', '16', '22'],
          'ptc22.tr2':  ['23', '24', '35']}

"""Per-track list of blankscreen, mseq, natscene, and driftbar recordings"""
BSMSNSDBRIDS = mergeuniquedictvals([BSRIDS, MSRIDS, NSRIDS, DBRIDS])
"""Per-track list of mseq, driftbar, driftgrating, and flashed grating recordings"""
MSDBDGFGRIDS = mergeuniquedictvals([MSRIDS, DBRIDS, DGRIDS, FGRIDS])
"""Per-track list of blankscreen, natscene, mseq, driftbar, driftgrating, and flashed grating
recordings"""
BSNSMSDBDGFGRIDS = mergeuniquedictvals([BSRIDS, NSRIDS, MSRIDS, DBRIDS, DGRIDS, FGRIDS])


NULLDIN = 65535 # integer value in stimulus .din files used as NULL (stimulus off)
ALPHA = 0.05 # threshold p value for rejecting null hypothesis

# mapping of recording absname to list of desynched and synched tranges, in that order.
# tranges are in us relative to start of ADC clock:
REC2STATETRANGES = {'ptc17.tr2b.r58': [(5.7e6, 700e6), # desynched trange, 66 Hz refresh rate
                                       (800e6, 1117.1e6)], # synched trange, 66 Hz refresh rate
                    'ptc18.tr1.r38':  [(43.8e6, 425e6), # desynched trange, ends ~ trial 76
                                       (550e6, 2243.8e6)], # synched trange, starts ~ trial 98
                    'ptc18.tr2c.r58': [(49e6, 750e6), # desynched trange
                                       (1000e6, 2248.9e6)], # synched trange
                    'ptc22.tr1.r08':  [(11e6, 1500e6), # desynched trange
                                       (1550e6, 2329.9e6)], # synched trange
                    'ptc22.tr1.r10':  [(1480e6, 2330.9e6), # desynched trange
                                       (12.1e6, 1400e6)], # synched trange
                    'ptc22.tr4b.r49': [(12.7e6, 1475e6), # desynched trange
                                       (1500e6, 2331.6e6)], # synched trange
                   }
