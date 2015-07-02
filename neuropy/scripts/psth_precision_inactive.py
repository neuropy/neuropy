"""Detect response events within PSTHs of specific inactive neurons, and plot them. Run from
within neuropy using `run -i scripts/psth_precision_inactive.py`"""

from __future__ import division, print_function

from core import intround

from psth_funcs import get_psth_peaks_simple, plot_psth

rec = ptc22.tr1.r08
strange = 1550e6, np.inf # r08 synched, us, end is ~ 2300s
nids = [5, 23, 24] # 3 example inactive yet responsive nids in ptc22.tr1.r08

EPS = np.spacing(1) # epsilon, smallest representable non-zero number

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
# 2.5 Hz thresh is 1 spike in the same 20 ms wide bin every 20 trials, assuming 0 baseline:
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
FWFRACTION = 0.5 # full width fraction of max
WIDTHMAX = 200 # maximum width, ms
WIDTHMAXPOINTS = intround(WIDTHMAX / 1000 / TRES) # maximum width, number of PSTH timepoints

# plotting params:
PLOTPSTH = True
FIGSIZE = 3.14, 2
YMAX = 6 # Hz
YTICKS = 0, 3, 6
MS = 5

t, psths, spikets = rec.psth(nids=nids, natexps=False, strange=strange, plot=False,
                             binw=BINW, tres=TRES, norm='ntrials')

psthparams = {} # params returned for each PSTH by get_psth_peaks
for nid, psth in zip(nids, psths):
    psthparams[nid] = get_psth_peaks_simple(t, psth, nid, WIDTHMAXPOINTS, minthresh=MINTHRESH,
                                            medianx=MEDIANX, fwfraction=FWFRACTION)
    if PLOTPSTH:
        plot_psth(psthparams, nid, fmt='k-', ms=MS, ymax=YMAX, yticks=YTICKS, figsize=FIGSIZE)
    #t, psth, thresh, baseline, peakis, lis, ris = psthparams[nid] # unpack

pl.show()
