"""Scatter plot MUA vs LFP SI for all specified recordings.
Run from within neuropy using `run -i scripts/muasi.py`"""

from __future__ import division

import scipy

# if choosing a cell spiketype or rftype, they will be picked active cells:
neurons = 'complex' # None (active), 'all', 'quiet', 'fast', 'slow', 'fastasym', 'slowasym',
               # 'simple', 'complex', 'LGN', 'unknown'
layers = True
rectype = 'bs'
kind = 'L/(L+H)'
axeslabels = False
ticklabels = False
figsize = 2.4, 2.4
pad = 0.3
ms = 3
ls = '--'
lw = 2
alpha = 0.6

rectype2rids = {'bs':BSRIDS, 'ns':NSRIDS, 'db':DBRIDS,
                'dg':DGRIDS, 'ms':MSRIDS, 'art':MSDBDGFGRIDS}

# choose desired tracks, need to be loaded ahead of time:
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
#tracks = [ptc22.tr1, ptc22.tr2]

recs = []
if rectype in rectype2rids:
    # collect specific subset of recordings, by recording type
    ridtype = rectype2rids[rectype]
    for track in tracks:
        rids = ridtype[track.absname]
        for rid in rids:
            recs.append(track.r[rid])
    recs = np.hstack(recs)
else:
    # collect all recordings:
    recs = []
    for track in tracks:
        # may as well do this in rid order
        rids = sorted(track.r)
        recs.append([ track.r[rid] for rid in rids ])
    recs = np.hstack(recs)

# override recordings if desired:
#recs = [ptc22.tr1.r['08'], ptc22.tr1.r['09'], ptc22.tr1.r['10']]
muas, sis = [], []
for rec in recs:
    mua, si = rec.mua_lfpsi(neurons=neurons, muawidth=30, muatres=10, lfpwidth=30, lfptres=10,
                            lfpsiwidth=30, lfpsitres=10, plot=False)
    muas.append(mua)
    sis.append(si)
mua = np.hstack(muas)
si = np.hstack(sis)

figure(figsize=figsize)

# calculate linear regression:
if layers:
    m1, b1, r1, p1, stderr1 = scipy.stats.linregress(si, mua[1]) # sup
    #m2, b2, r2, p2, stderr2 = scipy.stats.linregress(si, mua[2])
    m3, b3, r3, p3, stderr3 = scipy.stats.linregress(si, mua[3]) # deep
else:
    m0, b0, r0, p0, stderr0 = scipy.stats.linregress(si, mua[0]) # all

# scatter plot MUA vs SI:
if layers:
    plot(si, mua[1], 'r.', ms=ms, mew=0, alpha=alpha, label='sup, m=%.3f, r=%.3f' % (m1, r1))
    #plot(si, mua[2], 'g.', ms=ms, alpha=alpha, label='mid, m=%.3f, r=%.3f' % (m2, r2))
    plot(si, mua[3], 'b.', ms=ms, mew=0, alpha=alpha, label='deep, m=%.3f, r=%.3f' % (m3, r3))
else:
    plot(si, mua[0], 'e.', ms=ms, mew=0, alpha=alpha, label='all, m=%.3f, r=%.3f' % (m0, r0))

# plot linear regressions
sirange = np.array([0, 1])
if layers:
    plot(sirange, m1*sirange+b1, 'r', ls=ls, lw=lw, alpha=0.9) # sup
    #plot(sirange, m2*sirange+b2, 'g', ls=ls, lw=lw)
    plot(sirange, m3*sirange+b3, 'b', ls=ls, lw=lw, alpha=0.9) # deep
else:
    plot(sirange, m0*sirange+b0, 'k', ls=ls, lw=lw, alpha=0.9) # all

xmin, xmax = 0, 1
ymin, ymax = 0, 5
xlim(xmin, xmax)
ylim(ymin, ymax)
xt = np.arange(xmin, xmax+0.25, 0.25)
yt = np.arange(ymin, ymax+1, 1)
gca().set_xticks(xt)
gca().set_yticks(yt)
if not ticklabels:
    #tick_params(labelbottom=False)
    #tick_params(labelleft=False)
    gca().set_xticklabels([])
    gca().set_yticklabels([])
if axeslabels:
    xlabel('SI (%s)' % kind)
    ylabel('MUA (Hz/neuron)')

# add info text in upper left corner:
if layers:
    text(0.01, 0.99, 'super: r=%.2f, p=%.2g' % (r1, p1), color='r',
         transform=gca().transAxes, horizontalalignment='left', verticalalignment='top')
    #text(0.01, 0.99, 'mid: r=%.2f, p=%.2g' % (r1, p1), color='r',
    #     transform=gca().transAxes, horizontalalignment='left', verticalalignment='top')
    text(0.01, 0.91, 'deep: r=%.2f, p=%.2g' % (r3, p3), color='b',
         transform=gca().transAxes, horizontalalignment='left', verticalalignment='top')
else:
    text(0.01, 0.99, 'all: r=%.2f, p=%.2g' % (r0, p0), color='k',
         transform=gca().transAxes, horizontalalignment='left', verticalalignment='top')

tracknames = ', '.join([track.absname for track in tracks])
gcfm().window.setWindowTitle('muasi rectype=%s, neurons=%s, layers=%s, tracks=[%s]'
                             % (rectype, neurons, layers, tracknames))
tight_layout(pad=pad)

show()
