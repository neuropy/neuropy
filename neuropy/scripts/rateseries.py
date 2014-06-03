"""Plot each neuron's firing rate (binned at some fairly coarse resolution) over the span of
an entire track, repeated for each track, so 3 full width panels. Run from within neuropy
using `run -i scripts/rateseries.py`"""

from __future__ import division
from __future__ import print_function

from core import split_tranges
from colour import CCBLACKDICT0, CCWHITEDICT0 # for plotting on black or white

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
figsize = (10, 3)
alpha = 1
maxrate = 50 # max rate to display, Hz
widthsec = 5*60 # bin width, in sec
tressec = 60 # tres, in sec
width = widthsec * 1e6 # bin width, in us
tres = tressec * 1e6 # time resolution of potentially overlapping bins, in us
# any MPL cmap like jet_r or gist_rainbow or hsv or spectral_r. hsv is problematic because it's
# red at both extremes. None cycles through colours in CCDICT:
cmap = None #cm.gist_rainbow
cmapmax = 1.0 # use only first cmapmax of cmap, useful for hsv to prevent red at both extremes
CCDICT = CCBLACKDICT0
bg = 'k'
lw = 0.5
trtrange = [0, 5] # track trange to plot, in hours. Set to None to plot entire tracks

rates = {} # one entry per track
for track in tracks:
    trackrates = []
    # to prevent AssertionError in split_tranges(), filter out track tranges < width wide,
    # also to prevent unnecessary plotting, exclude track tranges that fall completely
    # outside trtrange:
    rectranges = []
    for trange in track.tranges:
        if trange[1] - trange[0] >= width: # trange is sufficiently wide
            if trtrange and trange[0] < trtrange[1]*3600*1e6: # trange starts within trtrange
                rectranges.append(trange)
    rectranges = np.array(rectranges)
    nids = np.sort(track.alln.keys())
    nn = len(nids)
    figure(figsize=figsize)
    axes(axisbg=bg) # set background color
    # plot rates for each rectrange separately, so lines aren't drawn across time gaps:
    for rectrange in rectranges: # one per recording
        recrates = [] # one row per neuron for this recording
        tranges = split_tranges([rectrange], width, tres) # possibly overlapping bins, in us
        midtranges = tranges.mean(axis=1) / 1e6 / 3600 # midpoints of tranges, in hours
        for nidi, nid in enumerate(nids): # nid order is (or should be) also depth order
            n = track.alln[nid]
            spikeis = n.spikes.searchsorted(tranges) # slice indices into spikes
            nspikes = spikeis[:, 1] - spikeis[:, 0] # number of spikes in each trange
            rate = nspikes / (width / 1e6) # spike rate per bin, in Hz
            rate[rate == 0.0] = np.nan # replace 0s with nans so they're ignored by plot()
            if cmap:
                cmapi = nidi/nn # from 0 to just under 1, cmaps wrap at 1
                c = cmap(cmapmax*cmapi)
            else:
                c = CCDICT[nidi] # use nidi to maximize colour alternation
            plot(midtranges, rate, '-', lw=lw, c=c, alpha=alpha)
            recrates.append(rate)
        recrates = np.vstack(recrates)
        trackrates.append(recrates)
        logmeanrate = np.nansum(log10(recrates), axis=0) / len(nids) # mean across neurons
        # underplot fat logmeanrate line in grey, take 10^ cuz yscale is set to log below:
        plot(midtranges, 10**logmeanrate, '-', lw=lw*10, c='w', alpha=0.5)
    trackrates = np.hstack(trackrates) # one row per neuron
    minrate = np.nanmin(trackrates) # lowest (non-nan) rate
    yscale('log')
    xlabel('time (hours)')
    ylabel('firing rate (Hz)')
    if trtrange:
        xlim(trtrange)
    ylim(minrate, maxrate)
    if cmap:
        cmapstr = cmap.name
    else:
        cmapstr = cmap # presumable None
    titlestr = ("%s_width=%d_tres=%d_trange=%s_cmap=%s_cmapmax=%.1f"
                % (track.absname, widthsec, tressec, trtrange, cmapstr, cmapmax))
    gcfm().window.setWindowTitle(titlestr)
    tight_layout(pad=0.3)
    rates[track.absname] = trackrates

show()
