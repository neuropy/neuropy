"""Plot each neuron's firing rate (binned at some fairly coarse resolution) over the span of
an entire track, repeated for each track, so 3 full width panels. Run from within neuropy
using `run -i scripts/rateseries.py`"""

from __future__ import division
from __future__ import print_function

from core import split_tranges
from colour import CCBLACKDICT0, CCWHITEDICT0 # for plotting on black or white
from scipy.stats import nanmean

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
figsize = (10, 3)
alpha = 1
maxrate = 50 # max rate to display, Hz
widthsec = 5*60 # bin width, in sec
tressec = 60 # time resolution, in sec
width = widthsec * 1e6 # bin width, in us
tres = tressec * 1e6 # time resolution of potentially overlapping bins, in us
# any MPL cmap like jet_r or gist_rainbow or hsv or spectral_r. hsv is problematic because it's
# red at both extremes. None cycles through colours in CCDICT:
cmap = cm.gist_rainbow
cmapmax = 1.0 # use only first cmapmax of cmap, useful for hsv to prevent red at both extremes
CCDICT = CCBLACKDICT0
bg = 'k'
lw = 0.5
trtrange = None #[0, 5] # track trange to plot, in hours. Set to None to plot entire tracks
plotlogmeanhists = True

rates, logmeans = {}, {} # one entry per track
fnanss, meanratess = [], [] # accumulate over all tracks
for track in tracks:
    trackrates = []
    meanrates = []
    logmeanrates = []
    # to prevent AssertionError in split_tranges(), filter out track tranges < width wide,
    # also to prevent unnecessary plotting, exclude track tranges that fall completely
    # outside trtrange:
    rectranges = []
    for trange in track.tranges:
        if trange[1] - trange[0] >= width: # trange is sufficiently wide
            if trtrange:
                if trange[0] < trtrange[1]*3600*1e6: # trange starts within trtrange
                    rectranges.append(trange)
            else:
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
        #nmuspikes = np.zeros(len(tranges), dtype=np.int64) # count multiunit spikes per bin
        for nidi, nid in enumerate(nids): # nid order is (or should be) also depth order
            n = track.alln[nid]
            spikeis = n.spikes.searchsorted(tranges) # slice indices into spikes
            nspikes = spikeis[:, 1] - spikeis[:, 0] # number of spikes in each trange
            #nmuspikes += nspikes # accumulate
            rate = nspikes / (width / 1e6) # spike rate per bin, in Hz
            rate[rate == 0.0] = np.nan # replace 0s with nans so they're ignored by plot()
            if cmap:
                cmapi = nidi / nn # from 0 to just under 1, cmaps wrap at 1
                c = cmap(cmapmax*cmapi)
            else:
                c = CCDICT[nidi] # use nidi to maximize colour alternation
            plot(midtranges, rate, '-', lw=lw, c=c, alpha=alpha)
            recrates.append(rate)
        recrates = np.vstack(recrates)
        trackrates.append(recrates)
        ## NOTE: np.nansum replaces nans with 0. scipy.stats.nanmean ignores nans completely
        # plot (arithmetic) meanrate line in transparent red:
        meanrate = np.nansum(recrates, axis=0) / nn # mean across neurons, replace nan with 0
        meanrates.append(meanrate)
        #plot(midtranges, meanrate, '-', lw=lw*10, c='r', alpha=0.5)
        # plot (geometric) logmeanrate line in transparent grey, take 10^ due to log yscale:
        logmeanrate = nanmean(log10(recrates), axis=0) # logmean across neurons, exclude nans
        logmeanrates.append(logmeanrate)
        plot(midtranges, 10**logmeanrate, '-', lw=lw*10, c='w', alpha=0.5)
        # plot multiunit rate in transparent blue, this is identical to the meanrate:
        #murate = nmuspikes / (width / 1e6) / nn # multiunit rate per bin per neuron, in Hz
        #plot(midtranges, murate, '-', lw=lw*10, c='b', alpha=0.5)
    trackrates = np.hstack(trackrates) # one row per neuron
    meanrates = np.hstack(meanrates) # concatenate meanrates from all recs in this track
    logmeanrates = np.hstack(logmeanrates) # concatenate logmeans from all recs in this track
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
    logmeans[track.absname] = logmeanrates
    if plotlogmeanhists:
        figure(figsize=(3, 3))
        nbins = 20 #intround(np.sqrt(len(logmeanrates)))
        logmin, logmax = -2, 0
        edges = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
        y = hist(10**logmeanrates, bins=edges, color='k')[0]
        ymax = max(y)
        xscale('log')
        ylim(0, ymax)
        yticks((0, ymax))
        xlabel('log-average firing rate (Hz)')
        ylabel('time bin count')
        titlestr = ("%s_logmeanhist_width=%d_tres=%d_trange=%s"
                    % (track.absname, widthsec, tressec, trtrange))
        gcfm().window.setWindowTitle(titlestr)
        tight_layout(pad=0.3)
        '''
        # plot meanrates vs fraction of cells at 0 Hz:
        figure(figsize=(3, 3))
        # fraction of cells with no rate, as a f'n of time:
        fnans = np.isnan(trackrates).sum(axis=0) / nn
        plot(fnans, meanrates/meanrates.mean(), 'k.', alpha=0.3)
        #yscale('log')
        fnanss.append(fnans)
        #normmeanrates.append(meanrates)
        meanratess.append(meanrates)
        '''
'''
fnanss = np.hstack(fnanss)
meanratess = np.hstack(meanratess)
figure(figsize=(3, 3))
plot(fnanss, meanratess/meanratess.mean(), 'k.', alpha=0.3)
#yscale('log')
xlabel('fraction of cells not firing')
ylabel('meanrate relative to meanrate.mean()')
'''
show()
