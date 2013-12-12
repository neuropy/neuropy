"""Plot spike paramaters vs. time of all spikes from each neuron in specified tracks. Params
are sx, x0, y0. Copy and paste into neuropy console"""

from colour import CCBLACKDICT0, CCWHITEDICT0 # for plotting on black or white

# style: 'points' or 'lines'. Binned lines have less detail but better visibility:
style = 'lines'
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
rc['savefig.format'] = 'png' # too many points for pdf
sxylim = 0, 160 # um
sxyticks = range(0, sxylim[1]+40, 40)
xposylim = -100, 100 # um
xposyticks = -50, 0, 50 # um
bg = 'black'
CCDICT = CCBLACKDICT0
#bg = 'white'
#CCDICT = CCWHITEDICT0


# posx and posy axes left, bottom, width, height and spacing (inches):
vs = 0.15 # vertical spacing
hs = 0.75 # horizontal spacing
ths = 0.15 # tight horizontal spacing, with no ylabel
l  = 1 # posx and posy left (for first track), also left margin
rm = 0.1 # right margin
sxb = 0.85 # sx bottom, also bottom margin
sxh = 2 # sx height
xb = sxb+sxh+vs # posx bottom, also bottom margin
xh = 2 # posx height
yb = xb+xh+vs # posy bottom
yh = 8 # posy height
tm = 0.35 # top margin

xlbwh, ylbwh, sxlbwh = [], [], []
spikess, maxts, xtickss = [], [], []
for tracki, spikefname in enumerate(spikefnames):
    spikes = np.load(spikefname) # record array
    spikess.append(spikes)
    maxt = spikes['t'].max() / 1e6 / 3600 # convert from us to hours
    xticks = range(0, intround(maxt), 2) # steps of 2 hours
    maxts.append(maxt)
    xtickss.append(xticks)
    w = maxt * 1/3 # inches
    xlbwh.append((l, xb, w, xh))
    ylbwh.append((l, yb, w, yh))
    sxlbwh.append((l, sxb, w, sxh))
    # inc l for next track (next column of axes)
    if tracki == 1:
        l += w + ths
    else:
        l += w + hs

fw = l - hs + rm # inches
fh = yb+yh+tm # inches
f = figure(figsize=(fw, fh))

# convert to fractional figure units, annoying:
sxlbwh = np.float64(sxlbwh)
xlbwh = np.float64(xlbwh)
ylbwh = np.float64(ylbwh)
sxlbwh[:, [0,2]] /= fw; sxlbwh[:, [1,3]] /= fh
xlbwh[:, [0,2]] /= fw; xlbwh[:, [1,3]] /= fh
ylbwh[:, [0,2]] /= fw; ylbwh[:, [1,3]] /= fh

sxas, xas, yas = [], [], []
for tracki, track in enumerate(tracks):
    spikes = spikess[tracki]
    maxt = maxts[tracki]
    xticks = xtickss[tracki]
    yposylim = yposylims[tracki]
    # manually position each axes, one track per column:
    sxa = f.add_axes(sxlbwh[tracki], axisbg=bg)
    xa = f.add_axes(xlbwh[tracki], axisbg=bg)
    ya = f.add_axes(ylbwh[tracki], axisbg=bg)
    nids, ts = sorted(track.alln), spikes['t']
    sxs, x0s, y0s = spikes['sx'], spikes['x0'], spikes['y0']
    if style == 'lines': # generate time bins of bw and tres, in hours:
        t0 = ts[0] / 1e6 / 3600 # convert from us to hours
        t1 = ts[-1] / 1e6 / 3600 # convert from us to hours
        tranges = core.split_tranges([(t0, t1)], bw, tres)
    # plot data for each nid, one at a time:
    for nidi, nid in enumerate(nids):
        sids, = np.where(spikes['nid'] == nid)
        t, sx, x, y = ts[sids], sxs[sids], x0s[sids], y0s[sids]
        t = t / 1e6 / 3600 # convert from us to hours
        c = CCDICT[nidi] # use nidi to maximize colour alternation
        if style == 'points':
            sxa.plot(t, sx, '.', ms=1, c=c)
            xa.plot(t, x, '.', ms=1, c=c)
            ya.plot(t, y, '.', ms=1, c=c)
        elif style == 'lines':
            tiranges = t.searchsorted(tranges)
            # split each series into values that fall within each bin (defined by its trange),
            # then take mean of each bin, using np.nan for empty bins:
            binsx = [ sx[tirange[0]:tirange[1]] for tirange in tiranges ]
            binsx = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in binsx ])
            biny = [ y[tirange[0]:tirange[1]] for tirange in tiranges ]
            biny = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in biny ])
            binx = [ x[tirange[0]:tirange[1]] for tirange in tiranges ]
            binx = np.array([ vals.mean() if len(vals) > 0 else np.nan for vals in binx ])
            # take left edges of tranges as bin time:
            bint = tranges[:, 0]
            sxa.plot(bint, binsx, '-', ms=None, lw=1, c=c)
            xa.plot(bint, binx, '-', ms=None, lw=1, c=c)
            ya.plot(bint, biny, '-', ms=None, lw=1, c=c)
            
    sxa.set_xlim(0, maxt)
    xa.set_xlim(0, maxt)
    ya.set_xlim(0, maxt)
    sxa.set_ylim(sxylim)
    xa.set_ylim(xposylim)
    ya.set_ylim(yposylim)
    ya.invert_yaxis()
    sxa.set_xticks(xticks)
    sxa.set_xlabel('time (hours)')
    if tracki == 0:
        sxa.set_ylabel('spread ($\mu$m)')
        xa.set_ylabel('x position ($\mu$m)')
        ya.set_ylabel('y position ($\mu$m)')
    sxa.set_yticks(sxyticks)
    xa.set_yticks(xposyticks)
    xa.set_xticklabels([])
    ya.set_xticklabels([])
    if tracki == 0: # leftmost column of axes, give them the same y label x position
        sxa.yaxis.set_label_coords(-0.11, 0.5) # fractional axes coords
        xa.yaxis.set_label_coords(-0.11, 0.5) # fractional axes coords
        ya.yaxis.set_label_coords(-0.11, 0.5)
    elif tracki == 1: # keep ylabels
        pass
    else: # remove y labels
        sxa.set_yticklabels([])
        xa.set_yticklabels([])
        ya.set_yticklabels([])
        
    ya.set_title(track.absname, y=1.005) # move it up slightly
    sxas.append(xa)
    xas.append(xa)
    yas.append(ya)

if style == 'points':
    titlestr = 'sparams'
elif style == 'lines':
    titlestr = 'sparams_lines'
f.canvas.manager.set_window_title(titlestr)
