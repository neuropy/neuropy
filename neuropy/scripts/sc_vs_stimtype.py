"""Scatter plot some summary statistic of spike correlations of each recording vs what
stimulus group each recording falls into.

Run by calling `%run -i scripts/sc_vs_stimtype.py` within neuropy"""

## TODO: weight corrs summary statistic by the number of spikes and/or cell pairs in
## each recording

## TODO: for each pair of recordings, find common subset of active neurons and calculate
## pairwise corrs for each recording in that pair using just those neurons

## TODO: designate recording type by colour, plot median corrs on x axis vs something else,
## like median SI or MUA on y axis 

#import os
#import numpy as np

#from animal import Animal
#from globals import DATAPATH
from pylab import get_current_fig_manager as gcfm
figsize = (7.5, 6.5)
'''
try:
    ptc22.tr1;
except NameError:
    ptc22 = Animal(os.path.join(DATAPATH, 'ptc22'))
    ptc22.load('tr1')
'''

tr = ptc22.tr1
blank_mseq_rids = ['04', '07', '09', '11', '17', '21']
mov_drift_rids = ['03', '05', '06', '08', '10', '18', '19', '20']
'''
tr = ptc22.tr2
blank_mseq_rids = ['26', '27', '32', '34', '36']
mov_drift_rids = ['25', '28', '31', '33']
'''

corrs = {}
for rid in (blank_mseq_rids + mov_drift_rids):
    r = tr.r[rid]
    sc = r.sc()
    sc.calc()
    corrs[rid] = np.median(sc.corrs)

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
a.set_xlabel('median spike correlations: blankscreen and mseq')
a.set_ylabel('median spike correlations: movie and drift bar')
winstr = 'sc_vs_stimtype.py_%s' % tr.absname
gcfm().window.setWindowTitle(winstr)
titlestr = 'sc_vs_stimtype.py: %s' % tr.absname
a.set_title(titlestr)
f.tight_layout(pad=0.3) # crop figure to contents
f.show()
