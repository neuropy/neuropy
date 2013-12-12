"""Plot spike paramaters vs. time of all spikes from each neuron in specified tracks. Params
are Vpp, sx, x0, y0. Run with 'run -i scripts/sparams.py' or copy and paste into neuropy
console"""

from __future__ import division
from __future__ import print_function

from colour import CCBLACKDICT0, CCWHITEDICT0 # for plotting on black or white

# style: 'points' or 'lines'. Binned lines have less detail but better visibility:
style = 'lines'
if style == 'points':
    rc['savefig.format'] = 'png' # too many points for pdf
elif style == 'lines':
    rc['savefig.format'] = 'pdf' # restore to default pdf
    bw = 10/60 # bin width for lines style, hours
    tres = 10/60 # bin tres for lines style, hours

# define stuff associated with each desired track:
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
# .spike files saved by spyke, these are what have all the required spike-level data:
spikefnames = ['/home/mspacek/data/ptc15/tr7c/track7c.track_2012-08-07_08.59.10.spike',
               '/home/mspacek/data/ptc22/tr1/track1.track_2012-03-02_15.56.59.spike',
               '/home/mspacek/data/ptc22/tr2/track2.track_2012-10-04_11.00.00.spike']
yposylims = (0, 1800), (0, 1250), (0, 1250) # um
#yposylims = (0, 1800), (0, 1800), (0, 1800) # um

fontsize(18)
vppylim = 0, 600
vppyticks = 0, 300, 600 # uV
sxylim = 0, 150 # um
sxyticks = 50, 100, 150 # um
xposylim = -100, 100 # um
xposyticks = -50, 0, 50 # um
bg = 'black'
CCDICT = CCBLACKDICT0
#bg = 'white'
#CCDICT = CCWHITEDICT0


# Vpp, sx, posx and posy axes left, bottom, width, height and spacing (inches):
vs = 0.15 # vertical spacing
hs = 0.75 # horizontal spacing
ths = 0.15 # tight horizontal spacing, with no ylabel
l  = 1 # left position (for first track), also left margin
rm = 0.1 # right margin
vppb = 0.85 # Vpp bottom, also bottom margin
vpph = 1.25 # Vpp height
sxb = vppb+vpph+vs # sx bottom
sxh = 1.25 # sx height
xb = sxb+sxh+vs # posx bottom
xh = 1.25 # posx height
yb = xb+xh+vs # posy bottom
yh = 8 # posy height
tm = 0.35 # top margin

vpplbwh, sxlbwh, xlbwh, ylbwh  = [], [], [], []
spikess, maxts, xtickss = [], [], []
for tracki, spikefname in enumerate(spikefnames):
    spikes = np.load(spikefname) # record array
    spikess.append(spikes)
    maxt = spikes['t'].max() / 1e6 / 3600 # convert from us to hours
    xticks = range(0, intround(maxt), 2) # steps of 2 hours
    maxts.append(maxt)
    xtickss.append(xticks)
    w = maxt * 1/3 # inches
    vpplbwh.append((l, vppb, w, vpph))
    sxlbwh.append((l, sxb, w, sxh))
    xlbwh.append((l, xb, w, xh))
    ylbwh.append((l, yb, w, yh))
    # inc l for next track (next column of axes)
    if tracki == 1:
        l += w + ths
    else:
        l += w + hs

fw = l - hs + rm # inches
fh = yb+yh+tm # inches
f = figure(figsize=(fw, fh))

# convert to fractional figure units, annoying:
vpplbwh = np.float64(vpplbwh)
sxlbwh = np.float64(sxlbwh)
xlbwh = np.float64(xlbwh)
ylbwh = np.float64(ylbwh)
vpplbwh[:, [0,2]] /= fw; vpplbwh[:, [1,3]] /= fh
sxlbwh[:, [0,2]] /= fw; sxlbwh[:, [1,3]] /= fh
xlbwh[:, [0,2]] /= fw; xlbwh[:, [1,3]] /= fh
ylbwh[:, [0,2]] /= fw; ylbwh[:, [1,3]] /= fh

vppas, sxas, xas, yas = [], [], [], []
for tracki, track in enumerate(tracks):
    spikes = spikess[tracki]
    maxt = maxts[tracki]
    xticks = xtickss[tracki]
    yposylim = yposylims[tracki]
    # manually position each axes, one track per column:
    vppa = f.add_axes(vpplbwh[tracki], axisbg=bg)
    sxa = f.add_axes(sxlbwh[tracki], axisbg=bg)
    xa = f.add_axes(xlbwh[tracki], axisbg=bg)
    ya = f.add_axes(ylbwh[tracki], axisbg=bg)
    nids, ts = sorted(track.alln), spikes['t']
    vpps, sxs, x0s, y0s = spikes['Vpp'], spikes['sx'], spikes['x0'], spikes['y0']
    if style == 'lines': # generate time bins of bw and tres, in hours:
        t0 = ts[0] / 1e6 / 3600 # convert from us to hours
        t1 = ts[-1] / 1e6 / 3600 # convert from us to hours
        tranges = core.split_tranges([(t0, t1)], bw, tres)
    # plot data for each nid, one at a time:
    for nidi, nid in enumerate(nids):
        sids, = np.where(spikes['nid'] == nid)
        t, vpp, sx, x, y = ts[sids], vpps[sids], sxs[sids], x0s[sids], y0s[sids]
        t = t / 1e6 / 3600 # convert from us to hours
        c = CCDICT[nidi] # use nidi to maximize colour alternation
        if style == 'points':
            vppa.plot(t, vpp, '.', ms=1, c=c)
            sxa.plot(t, sx, '.', ms=1, c=c)
            xa.plot(t, x, '.', ms=1, c=c)
            ya.plot(t, y, '.', ms=1, c=c)
        elif style == 'lines':
            tiranges = t.searchsorted(tranges)
            # split each series into values that fall within each bin (defined by its trange),
            # then take mean of each bin, using np.nan for empty bins:
            binvpp = [ vpp[tirange[0]:tirange[1]] for tirange in tiranges ]
            binvpp = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in binvpp ])
            binsx = [ sx[tirange[0]:tirange[1]] for tirange in tiranges ]
            binsx = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in binsx ])
            biny = [ y[tirange[0]:tirange[1]] for tirange in tiranges ]
            biny = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in biny ])
            binx = [ x[tirange[0]:tirange[1]] for tirange in tiranges ]
            binx = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in binx ])
            # take left edges of tranges as bin time:
            bint = tranges[:, 0]
            vppa.plot(bint, binvpp, '-', marker=None, lw=1, c=c)
            sxa.plot(bint, binsx, '-', marker=None, lw=1, c=c)
            xa.plot(bint, binx, '-', marker=None, lw=1, c=c)
            ya.plot(bint, biny, '-', marker=None, lw=1, c=c)
            # an alternative is to use smaller bins and plot disconnected points, but then
            # it's again more difficult to track a single neuron, as it is with
            # style='points':
            #vppa.plot(bint, binvpp, '.', ms=1, lw=1, c=c)
            #sxa.plot(bint, binsx, '.', ms=1, lw=1, c=c)
            #xa.plot(bint, binx, '.', ms=1, lw=1, c=c)
            #ya.plot(bint, biny, '.', ms=1, lw=1, c=c)
            
    vppa.set_xlim(0, maxt)
    sxa.set_xlim(0, maxt)
    xa.set_xlim(0, maxt)
    ya.set_xlim(0, maxt)
    vppa.set_ylim(vppylim)
    sxa.set_ylim(sxylim)
    xa.set_ylim(xposylim)
    ya.set_ylim(yposylim)
    ya.invert_yaxis()
    vppa.set_xlabel('time (hours)')
    sxa.set_xticks(xticks)
    if tracki == 0:
        vppa.set_ylabel('$V_{pp}$ ($\mu$V)')
        sxa.set_ylabel('$\sigma$ ($\mu$m)')
        xa.set_ylabel('x ($\mu$m)')
        ya.set_ylabel('y ($\mu$m)')
    vppa.set_yticks(vppyticks)
    sxa.set_yticks(sxyticks)
    xa.set_yticks(xposyticks)
    sxa.set_xticklabels([])
    xa.set_xticklabels([])
    ya.set_xticklabels([])
    if tracki == 0: # leftmost column of axes, give them the same y label x position
        vppa.yaxis.set_label_coords(-0.11, 0.5) # fractional axes coords
        sxa.yaxis.set_label_coords(-0.11, 0.5)
        xa.yaxis.set_label_coords(-0.11, 0.5)
        ya.yaxis.set_label_coords(-0.11, 0.5)
    elif tracki == 1: # keep ylabels
        pass
    else: # remove y labels
        vppa.set_yticklabels([])
        sxa.set_yticklabels([])
        xa.set_yticklabels([])
        ya.set_yticklabels([])
        
    ya.set_title(track.absname, y=1.005) # move it up slightly
    vppas.append(vppa)
    sxas.append(xa)
    xas.append(xa)
    yas.append(ya)

if style == 'points':
    titlestr = 'sparams'
elif style == 'lines':
    titlestr = 'sparams_lines'
f.canvas.manager.set_window_title(titlestr)
show()
