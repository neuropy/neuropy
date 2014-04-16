"""Plot distributions of synchrony index (SI) for various tracks. Run from within neuropy
using `run -i scripts/sihist.py`"""

from __future__ import division

figsize = (2.5, 2.5)
lfpwidth = 30 # sec
lfptres = 5 # sec
binw = 0.02 # SI
simin, simax = 0, 1
kind = 'L/(L+H)'

try:
    SICACHE
except NameError:
    SICACHE = {} # init, caches sicalc results

def calcsi(recs):
    sis = []
    for rec in recs:
        si, t = rec.lfp.si(kind=kind, lfpwidth=lfpwidth, lfptres=lfptres, plot=False)
        sis.append(si)
    return np.hstack(sis)

def sihist(recs=None, basetitle='', sis=None):
    """Plot SI histogram across all recordings in recs"""
    # collect SI signals from all recs:
    if sis == None:
        try:
            sis = SICACHE[basetitle]
        except KeyError:
            sis = calcsi(recs)
            SICACHE[basetitle] = sis # cache it
    else:
        pass # use sis provided by caller

    # plot SI histogram:
    figure(figsize=figsize)
    bins = np.arange(simin, simax+binw, binw) # left edges + rightmost edge
    n = np.histogram(sis, bins=bins)[0]
    bar(left=bins[:-1], height=n, width=binw, color='k', ec='k')
    xlim(xmin=simin, xmax=simax)
    ylim(ymax=n.max()) # max out at max bin height
    ticks = [0, 0.2, 0.4, 0.6, 0.8, 1] # minimal decimal places
    xticks(ticks, [str(t) for t in ticks])
    yticks([]) # turn off y ticks to save space
    title(basetitle)
    #yticks([0, n.max()])
    #xlabel('SI (%s)' % kind)
    gcfm().window.setWindowTitle(basetitle + '_sihist')
    tight_layout(pad=0.3)

# quick tests:
#sihist([ptc15.tr7c.r74], 'ptc15.tr7c.r74')
#sihist([ptc22.tr1.r08], 'ptc15.tr7c.r08')
#sihist([ptc22.tr1.r10], 'ptc15.tr7c.r10')

tracks = [ptc15.tr7c, ptc17.tr1, ptc17.tr2b, ptc18.tr1, ptc18.tr2c, ptc20.tr1, ptc20.tr2,
          ptc20.tr3, ptc21.tr2, ptc21.tr5c, ptc21.tr6b, ptc22.tr1, ptc22.tr2, ptc22.tr4b,
          ptc22.tr5b]
for track in tracks:
    # sorted() calls each rec's __cmp__ for sorting, not that recording order matters:
    sihist(sorted(track.r.values()), track.absname)

# plot histogram across all tracks:
sis = np.hstack([ SICACHE[track.absname] for track in tracks ])
sihist(basetitle='all tracks', sis=sis)
dthourall = np.sum([ track.dthour for track in tracks ])
print('all tracks: %.3f h' % dthourall)

show()
