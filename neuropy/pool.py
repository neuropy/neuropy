"""Analyses pooled across neuropy objects, typically across Recordings"""

import numpy as np
import scipy
import pylab as pl
from pylab import get_current_fig_manager as gcfm

from core import lastcmd
from track import Track


def sc_si(source, method='weighted mean', kind='ncv', layers=False, sirange=(-1, 1),
          figsize=(7.5, 6.5)):
    """Pool recording.sc().si() results across recordings specified by source,
    plot the result"""
    uns = get_ipython().user_ns

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
        corrs, si, ylabel = rec.sc().si(method=method, plot=False)
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
    if layers:
        ylim = corrs[:5].min(), corrs[:5].max()
    else:
        ylim = corrs[0].min(), corrs[0].max()
    yrange = ylim[1] - ylim[0]
    extra = yrange*0.03 # 3 %
    ylim = ylim[0]-extra, ylim[1]+extra

    # keep only those points whose synchrony index falls within sirange:
    if sirange == None:
        finitesi = si[np.isfinite(si)]
        sirange = finitesi.min(), finitesi.max()
    sirange = np.asarray(sirange)
    keepis = (sirange[0] <= si[0]) * (si[0] <= sirange[1]) # boolean index array
    si = si[:, keepis]
    corrs = corrs[:, keepis]
    # plot linear regressions of corrs vs si[0]:
    m0, b0, r0, p0, stderr0 = scipy.stats.linregress(si[0], corrs[0])
    a.plot(sirange, m0*sirange+b0, 'e--')
    if layers:
        m1, b1, r1, p1, stderr1 = scipy.stats.linregress(si[0], corrs[1])
        a.plot(sirange, m1*sirange+b1, 'r--')
        m2, b2, r2, p2, stderr2 = scipy.stats.linregress(si[0], corrs[2])
        a.plot(sirange, m2*sirange+b2, 'g--')
        m3, b3, r3, p3, stderr3 = scipy.stats.linregress(si[0], corrs[3])
        a.plot(sirange, m3*sirange+b3, 'b--')
        #m4, b4, r4, p4, stderr4 = scipy.stats.linregress(si[0], corrs[4])
        #a.plot(sirange, m4*sirange+b4, 'y--', zorder=0)

    # scatter plot corrs vs si, one colour per laminarity:
    a.plot(si[0], corrs[0], 'e.', label='all, m=%.3f, r=%.3f' % (m0, r0))
    if layers:
        a.plot(si[0], corrs[1], 'r.', label='superficial, m=%.3f, r=%.3f'
                                         % (m1, r1))
        a.plot(si[0], corrs[2], 'g.', label='middle, m=%.3f, r=%.3f'
                                         % (m2, r2))
        a.plot(si[0], corrs[3], 'b.', label='deep, m=%.3f, r=%.3f'
                                         % (m3, r3))
        #a.plot(si[0], corrs[4], 'y.', label='other, m=%.3f, r=%.3f'
        #                                 % (m4, r4), zorder=0)
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
