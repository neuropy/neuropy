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

rates = {}
for track in tracks:
    trackrates = []
    # filter out track tranges < width wide (like short CSD recordings):
    tranges = []
    for trange in track.tranges:
        if trange[1] - trange[0] >= width:
            tranges.append(trange)
    tranges = np.array(tranges)
    # split tranges, so that gaps between recordings are excluded from rate calcs
    tranges = split_tranges(tranges, width, tres) # possibly overlapping, in us
    midtranges = tranges.mean(axis=1) / 1e6 / 3600 # midpoints of tranges, in hours
    nids = np.sort(track.alln.keys())
    figure(figsize=figsize)
    axes(axisbg=bg) # set background color
    for nidi, nid in enumerate(nids):
        n = track.alln[nid]
        spikeis = n.spikes.searchsorted(tranges) # slice indices into spikes
        nspikes = spikeis[:, 1] - spikeis[:, 0] # number of spikes in each trange
        rate = nspikes / (width / 1e6) # spike rate per bin, in Hz
        #rate[rate == 0.0] = np.nan # replace 0's with nans so they're ignored by plot()
        c = CCDICT[nidi] # use nidi to maximize colour alternation
        plot(midtranges, rate, '-', lw=lw, c=c, alpha=alpha)
        trackrates.append(rate)
    yscale('log')
    xlabel('time (h)')
    ylabel('firing rate (Hz)')
    gcfm().window.setWindowTitle(track.absname)
    tight_layout(pad=0.3)
    trackrates = np.vstack(trackrates)
    minrate = np.unique(trackrates)[1] # lowest rate, excluding 0 Hz
    xlim(trtrange)
    ylim(ymin=minrate)
    rates[track.absname] = trackrates

show()
