"""Plot each neuron's firing rate (binned at some fairly coarse resolution) over the span of
an entire track, repeated for each track, so 3 full width panels. Run from within neuropy
using `run -i scripts/rateseries.py`"""

from __future__ import division
from __future__ import print_function

from core import split_tranges
from colour import CCBLACKDICT0, CCWHITEDICT0 # for plotting on black or white


tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
figsize = (8, 3)
alpha = 1
widthsec = 2*60 # bin width, in sec
tressec = 30 # tres, in sec
width = widthsec * 1e6 # bin width, in us
tres = tressec * 1e6 # time resolution of potentially overlapping bins, in us
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
    figure(figsize=figsize)
    axes(axisbg=bg) # set background color
    for nidi, nid in enumerate(nids):
        neuronrates = []
        n = track.alln[nid]
        # plot rates for each rectrange separately, so lines aren't drawn across time gaps:
        for rectrange in rectranges: # one per recording
            tranges = split_tranges([rectrange], width, tres) # possibly overlapping bins, in us
            midtranges = tranges.mean(axis=1) / 1e6 / 3600 # midpoints of tranges, in hours
            spikeis = n.spikes.searchsorted(tranges) # slice indices into spikes
            nspikes = spikeis[:, 1] - spikeis[:, 0] # number of spikes in each trange
            rate = nspikes / (width / 1e6) # spike rate per bin, in Hz
            rate[rate == 0.0] = np.nan # replace 0's with nans so they're ignored by plot()
            c = CCDICT[nidi] # use nidi to maximize colour alternation
            plot(midtranges, rate, '-', lw=lw, c=c, alpha=alpha)
            neuronrates.append(rate)
        trackrates.append(np.hstack(neuronrates))
    yscale('log')
    xlabel('time (h)')
    ylabel('firing rate (Hz)')
    trackrates = np.vstack(trackrates) # one row per neuron
    minrate = np.unique(trackrates)[1] # lowest rate, excluding 0 Hz
    if trtrange:
        xlim(trtrange)
    ylim(ymin=minrate)
    titlestr = "%s_width=%d_tres=%d_trange=%s" % (track.absname, widthsec, tressec, trtrange)
    gcfm().window.setWindowTitle(titlestr)
    tight_layout(pad=0.3)
    rates[track.absname] = trackrates

show()
