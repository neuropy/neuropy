"""Analyses pooled across neuropy objects, typically across Recordings"""

import numpy as np
import scipy
from scipy.stats import linregress
import pylab as pl
from pylab import get_current_fig_manager as gcfm

from core import lastcmd
from track import Track


def sc_si(source, method='mean', kind='ncv', layers=False, ms=1, sirange=(-1, 1),
          figsize=(7.5, 6.5)):
    """Pool recording.sc().si() results across recordings specified by source,
    plot the result"""
    uns = get_ipython().user_ns
    if layers == False:
        layers = ['all']
    elif layers == True:
        layers = ['sup', 'deep']
    LAYER2I = {'all':0, 'sup':1, 'mid':2, 'deep':3, 'other':4}
    layeris = [ LAYER2I[layer] for layer in layers ]

    # decide which recordings to pool across, collect recordings in a list:
    if type(source) == Track:
        rids = uns['RIDS'][source.absname]
        recs = [ source.r[rid] for rid in rids ]
    elif type(source) == list: # assume it's a list of recordings
        recs = source
    elif type(source) == dict: # assume it's a dict of {animal.tr: [rids]} key-value pairs
        recs = []
        for animal_track, rids in source.items():
            animalname, trackname = animal_track.split('.')
            tr = uns[animalname].__getattribute__(trackname)
            recs.extend([ tr.r[rid] for rid in rids ])

    # calculate
    corrss, sis = [], []
    for rec in recs:
        print(rec.absname)
        corrs, si, ylabel = rec.sc().si(method=method, kind=kind, plot=False)
        corrss.append(corrs)
        sis.append(si)
    corrs = np.hstack(corrss)
    si = np.hstack(sis)

    # plot
    if kind in ('L/(L+H)', 'L/H'):
        sisource = 'lfp'
    else:
        sisource = 'mua'

    f = pl.figure(figsize=figsize)
    a = f.add_subplot(111)

    #ylim = corrs[layeris].min(), corrs[layeris].max()
    #yrange = ylim[1] - ylim[0]
    #extra = yrange*0.03 # 3 %
    #ylim = ylim[0]-extra, ylim[1]+extra
    ylim = uns['SCLIMITS']

    # keep only those points whose synchrony index falls within sirange:
    if sirange == None:
        finitesi = si[np.isfinite(si)]
        sirange = finitesi.min(), finitesi.max()
    sirange = np.asarray(sirange)
    keepis = (sirange[0] <= si[0]) * (si[0] <= sirange[1]) # boolean index array
    si = si[:, keepis]
    corrs = corrs[:, keepis]
    # plot linear regressions of corrs vs si[0]:
    if 'all' in layers:
        m0, b0, r0, p0, stderr0 = linregress(si[0], corrs[0])
        a.plot(sirange, m0*sirange+b0, 'e--')
    if 'sup' in layers:
        m1, b1, r1, p1, stderr1 = linregress(si[0], corrs[1])
        a.plot(sirange, m1*sirange+b1, 'r--')
    if 'mid' in layers:
        m2, b2, r2, p2, stderr2 = linregress(si[0], corrs[2])
        a.plot(sirange, m2*sirange+b2, 'g--')
    if 'deep' in layers:
        m3, b3, r3, p3, stderr3 = linregress(si[0], corrs[3])
        a.plot(sirange, m3*sirange+b3, 'b--')
    if 'other' in layers:
        m4, b4, r4, p4, stderr4 = linregress(si[0], corrs[4])
        a.plot(sirange, m4*sirange+b4, 'y--', zorder=0)

    # scatter plot corrs vs si, one colour per laminarity:
    if 'all' in layers:
        a.plot(si[0], corrs[0], 'e.', ms=ms, label='all, m=%.3f, r=%.3f'
                                                   % (m0, r0))
    if 'sup' in layers:
        a.plot(si[0], corrs[1], 'r.', ms=ms, label='superficial, m=%.3f, r=%.3f'
                                                   % (m1, r1))
    if 'mid' in layers:
        a.plot(si[0], corrs[2], 'g.', ms=ms, label='middle, m=%.3f, r=%.3f'
                                                   % (m2, r2))
    if 'deep' in layers:
        a.plot(si[0], corrs[3], 'b.', ms=ms, label='deep, m=%.3f, r=%.3f'
                                                   % (m3, r3))
    if 'other' in layers:
        a.plot(si[0], corrs[4], 'y.', ms=ms, label='other, m=%.3f, r=%.3f'
                                                   % (m4, r4), zorder=0)
    #a.set_xlim(sirange)
    if sisource == 'lfp':
        a.set_xlim(0, 1)
    elif sisource == 'mua' and kind[0] == 'n':
        a.set_xlim(-1, 1)
    a.set_ylim(ylim)
    #a.autoscale(enable=True, axis='y', tight=True)
    a.set_xlabel(kind)
    a.set_ylabel(ylabel)
    titlestr = lastcmd()
    gcfm().window.setWindowTitle(titlestr)
    a.set_title(titlestr)
    a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
    f.tight_layout(pad=0.3) # crop figure to contents
