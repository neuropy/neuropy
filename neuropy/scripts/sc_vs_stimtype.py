"""Scatter plot some summary statistic of spike correlations of each recording vs what
stimulus group each recording falls into.

Run by calling `%run -i scripts/sc_vs_stimtype.py track_absname` within neuropy"""

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
width = 60 # sec, or None
tres = 60 # sec, or None
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

blank_mseq_corrs = []
mov_drift_corrs = []
for rid in (blank_mseq_rids + mov_drift_rids):
    r = tr.r[rid]
    sc = r.sc(width=width, tres=tres)
    sc.calc()
    totalcounts = sc.counts.sum(axis=0) # len(ntranges)
    if method == 'median':
        corrs = np.median(sc.corrs, axis=0)
    elif method == 'weighted median': # not entirely sure this is right:
        corrs = np.median(sc.corrs * sc.counts / totalcounts, axis=0) * sc.npairs
    elif method == 'weighted mean':
        corrs = (sc.corrs * sc.counts / totalcounts).sum(axis=0)
    else:
        raise ValueError("unknown method %r" % method)
    if rid in blank_mseq_rids:
        blank_mseq_corrs.append(corrs)
    else:
        mov_drift_corrs.append(corrs)
    #print('rid: %s, ntranges: %d, len(corrs): %d' % (rid, len(totalcounts), len(corrs)))

blank_mseq_corrs = np.hstack(blank_mseq_corrs)
mov_drift_corrs = np.hstack(mov_drift_corrs)
# repeat each element in blank_mseq_corrs len(mov_drift_corrs) times:
x = np.repeat(blank_mseq_corrs, len(mov_drift_corrs))
# tile mov_drift_corrs len(blank_mseq_corrs) times:
y = np.tile(mov_drift_corrs, len(blank_mseq_corrs))

f = pl.figure(figsize=figsize)
a = f.add_subplot(111)
lim = min([x.min(), y.min(), 0]), max([x.max(), y.max()])
a.plot(lim, lim, c='e', ls='--', marker=None) # y=x line
a.plot(x, y, 'k.')
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
