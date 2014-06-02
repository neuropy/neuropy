"""Plot each neuron's firing rate (binned at some fairly coarse resolution) over the span of
an entire track, repeated for each track, so 3 full width panels. Run from within neuropy
using `run -i scripts/rateseries.py`"""

from __future__ import division
from __future__ import print_function

from colour import CCBLACKDICT0, CCWHITEDICT0 # for plotting on black or white

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
figsize = (8, 3)
alpha = 1
binwmin = 5 # bin width, in min
binwhour = binwmin / 60
CCDICT = CCBLACKDICT0
bg = 'k'
lw = 0.5

rates = {}
for track in tracks:
    trackrates = []
    dthour = track.dthour
    bins = np.arange(0, dthour+binwhour, binwhour) # nonoverlapping bins, end inclusive
    nids = np.sort(track.alln.keys())
    figure(figsize=figsize)
    axes(axisbg=bg) # set background color
    for nidi, nid in enumerate(nids):
        n = track.alln[nid]
        spikes = n.spikes / 1e6 / 3600 # spike times, floats, in hours
        nspikes, edges = np.histogram(spikes, bins) # spike counts per bin
        rate = nspikes / binwhour / 3600 # spike rate per bin, in Hz
        rate[rate == 0.0] = np.nan # replace 0's with nans so that they're ignored by plot()
        c = CCDICT[nidi] # use nidi to maximize colour alternation
        plot(bins[:-1], rate, '-', lw=lw, c=c, alpha=alpha)
        trackrates.append(rate)
    yscale('log')
    xlim(0, dthour)
    xlabel('time (h)')
    ylabel('firing rate (Hz)')
    gcfm().window.setWindowTitle(track.absname + '_binwmin=%d' % binwmin)
    tight_layout(pad=0.3)
    trackrates = np.vstack(trackrates)
    minrate = np.unique(trackrates)[1] # lowest rate, excluding 0 Hz
    xlim(0, 6) # show only first 6 hours
    ylim(ymin=minrate)
    rates[track.absname] = trackrates

show()
