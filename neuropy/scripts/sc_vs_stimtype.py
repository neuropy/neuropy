"""Scatter plot some summary statistic of spike correlations of each recording vs what
stimulus group each recording falls into.

Run by calling `%run -i scripts/sc_vs_stimtype.py track_absname` within neuropy"""

## TODO: weight corrs summary statistic by the number of spikes and/or cell pairs in
## each recording. Or, split recordings up into pieces of x min long, and then plot
## corr of all pieces of one stim type vs those of the other

## TODO: for each pair of recordings, find common subset of active neurons and calculate
## pairwise corrs for each recording in that pair using just those neurons

## TODO: designate recording type by colour, plot median corrs on x axis vs something else,
## like median SI or MUA on y axis

## TODO: maybe limit to visually responsive cells

#import os
#import numpy as np
#from animal import Animal
#from globals import DATAPATH

import argparse
from pylab import get_current_fig_manager as gcfm

method = 'weighted mean' # 'median', 'weighted median', 'weighted mean'
width = None #60 # sec
tres = None #60 # sec
figsize = (7.5, 6.5)

parser = argparse.ArgumentParser()
parser.add_argument('track')
track_absname = parser.parse_args().track

if track_absname == 'ptc22.tr1':
    tr = ptc22.tr1
    blank_mseq_rids = ['04', '07', '09', '11', '17', '21']
    mov_drift_rids = ['03', '05', '06', '08', '10', '18', '19', '20']
elif track_absname == 'ptc22.tr2':
    tr = ptc22.tr2
    blank_mseq_rids = ['26', '27', '32', '34', '36']
    mov_drift_rids = ['25', '28', '31', '33']
else:
    raise ValueError("can't handle track %r" % track_absname)

corrs = {}
for rid in (blank_mseq_rids + mov_drift_rids):
    r = tr.r[rid]
    sc = r.sc(width=width, tres=tres)
    sc.calc()
    totalcounts = sc.counts.sum(axis=0) # len(ntranges)
    if method == 'median':
        corrs[rid] = np.median(sc.corrs, axis=0)
    elif method == 'weighted median': # not entirely sure this is right:
        corrs[rid] = np.median(sc.corrs * sc.counts / totalcounts, axis=0) * sc.npairs
    elif method == 'weighted mean':
        corrs[rid] = (sc.corrs * sc.counts / totalcounts).sum(axis=0)
    else:
        raise ValueError("unknown method %r" % method)

data = []
rpairs = []
for rid0 in blank_mseq_rids:
    for rid1 in mov_drift_rids:
        data.append((corrs[rid0], corrs[rid1]))
        rpairs.append((rid0, rid1))
data = np.asarray(data)
rpairs = np.asarray(rpairs)
print(data)
print(rpairs)

f = pl.figure(figsize=figsize)
a = f.add_subplot(111)
lim = min(data.min(), 0), data.max()
a.plot(lim, lim, c='e', ls='--', marker=None) # y=x line
a.plot(data[:, 0], data[:, 1], 'k.')
#a.set_xlim(lim)
#a.set_ylim(lim)
a.set_xlabel('%s spike correlations: blankscreen and mseq' % method)
a.set_ylabel('%s spike correlations: movie and drift bar' % method)
winstr = 'sc_vs_stimtype.py_%s' % tr.absname
gcfm().window.setWindowTitle(winstr)
titlestr = 'sc_vs_stimtype.py: %s' % tr.absname
a.set_title(titlestr)
f.tight_layout(pad=0.3) # crop figure to contents
f.show()
