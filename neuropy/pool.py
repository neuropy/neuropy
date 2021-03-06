"""Analyses pooled across neuropy objects, typically across Recordings"""

import numpy as np
import scipy
from scipy.stats import linregress
import pylab as pl
from pylab import get_current_fig_manager as gcfm

from core import lastcmd, parse_source


def sc_si(source, method='mean', sisource='lfp', kind=None, chani=-1, sirange=None,
          layers=False, ms=1, figsize=(7.5, 6.5)):
    """Pool recording.sc().si() results across recordings specified by source,
    plot the result"""
    uns = get_ipython().user_ns
    if layers == False:
        layers = ['all']
    elif layers == True:
        layers = ['sup', 'deep']
    LAYER2I = {'all':0, 'sup':1, 'mid':2, 'deep':3, 'other':4}
    layeris = [ LAYER2I[layer] for layer in layers ]

    recs, tracks = parse_source(source)

    if sisource not in ['lfp', 'mua']:
        raise ValueError('unknown sisource %r' % sisource)

    if kind == None:
        if sisource == 'lfp':
            kind = uns['LFPSIKIND']
        else:
            kind = uns['MUASIKIND']

    # calculate
    corrss, sis = [], []
    for rec in recs:
        print(rec.absname)
        corrs, si, ylabel = rec.sc().si(method=method, sisource=sisource, kind=kind,
                                        chani=chani, sirange=sirange, plot=False)
        corrss.append(corrs)
        sis.append(si)
    corrs = np.hstack(corrss)
    si = np.hstack(sis)

    # plot
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
    if kind[0] == 'n':
        a.set_xlim(-1, 1)
    a.set_ylim(ylim)
    #a.autoscale(enable=True, axis='y', tight=True)
    a.set_xlabel('%s SI (%s)' % (sisource.upper(), kind))
    a.set_ylabel(ylabel)
    titlestr = lastcmd()
    gcfm().window.setWindowTitle(titlestr)
    a.set_title(titlestr)
    a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
    f.tight_layout(pad=0.3) # crop figure to contents


def mua_si_lfp_si(source, layers=False, ms=1, figsize=(7.5, 6.5)):
    """Pool recording.mua_si_lfp_si() results across recordings specified by source,
    plot the result"""
    uns = get_ipython().user_ns
    recs, tracks = parse_source(source)
    lfpsis, muasis = [], []
    for rec in recs:
        print(rec.absname)
        lfpsi, muasi, t = rec.mua_si_lfp_si(ms=ms, layers=layers, plot=False, plotseries=False,
                                            figsize=figsize)
        lfpsis.append(lfpsi)
        muasis.append(muasi)
    lfpsi = np.hstack(lfpsis)
    muasi = np.hstack(muasis)
    # plot:
    f = pl.figure(figsize=figsize)
    a = f.add_subplot(111)
    a.plot([-1, 1], [-1, 1], 'e--') # underplot y=x line
    a.plot(lfpsi, muasi[0], 'e.', ms=ms)
    if layers:
        a.plot(lfpsi, muasi[1], 'r.', ms=ms)
        a.plot(lfpsi, muasi[2], 'g.', ms=ms)
        a.plot(lfpsi, muasi[3], 'b.', ms=ms)
    a.set_xlabel('LFP SI (%s)' % uns['LFPSIKIND'])
    a.set_ylabel('MUA SI (%s)' % uns['MUASIKIND'])
    a.set_xlim(-1, 1)
    a.set_ylim(-1, 1)
    titlestr = lastcmd()
    gcfm().window.setWindowTitle(titlestr)
    a.set_title(titlestr)
    f.tight_layout(pad=0.3) # crop figure to contents
    #return lfpsi, muasi, t


def sc_ising_vs_cch(source, ms=5, figsize=(7.5, 6.5)):
    """Scatter plot spike corrs calculated from Ising matrix against those calculated
    from CCH. INCOMPLETE.

    - find tracks in common, get allnids from each track

    - how to deal with big time gaps between experiments in a single recording? I constrain
    to the set of tranges of each experiment in rec.codes()

    - maybe i can convert the core.SpikeCorr object to take a source argument instead of
    recording/experiment objects
        - do all the spikecorr analyses make sense for multiple recordings, or for recordings
        from different tracks?

    - for each track absname
    """

    isingscs = {}
    cchscs = {}
    # init a dict
    # for each rec, find out which track it's from

    recs, tracks = parse_source(source)
    isingscs, cchscs = [], []
    for rec in recs:
        print(rec.absname)
        sc = rec.sc()
        sc.calc()
        isingscs.append(sc.corrs)
        cchscs.append(rec.sc_cch())
    # for something...

    isingsc = np.hstack(isingscs)
    cchsc = np.hstack(cchscs)

    # plot:
    f = pl.figure(figsize=figsize)
    a = f.add_subplot(111)
    a.plot(isingsc, cchsc, 'e.', ms=ms)
    a.set_xlabel('Ising spike corrs')
    a.set_ylabel('CCH spike corrs')
    a.set_xlim(-0.05, 0.2)
    a.set_ylim(-0.5, 1)
    titlestr = lastcmd()
    gcfm().window.setWindowTitle(titlestr)
    a.set_title(titlestr)
    f.tight_layout(pad=0.3) # crop figure to contents
