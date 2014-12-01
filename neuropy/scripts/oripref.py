"""Plot orientation preference vs tuning strength and cell position of all oriented cells from
a number of orientated stimulus recordings. Run with 'run -i scripts/oripref.py' or copy and
paste into neuropy console."""

from __future__ import division
from __future__ import print_function

import scipy

from polar_demo import fractional_polar_axes

# drift bar, drift grating and flashed grating recordings:
trackrecs = {ptc15.tr7c: [ptc15.tr7c.r71, ptc15.tr7c.r73, ptc15.tr7c.r85],
             ptc22.tr1: [ptc22.tr1.r03, ptc22.tr1.r13, ptc22.tr1.r14, ptc22.tr1.r18],
             ptc22.tr2: [ptc22.tr2.r25, ptc22.tr2.r30, ptc22.tr2.r31]}
#trackrecs = {ptc15.tr7c: [ptc15.tr7c.r71, ptc15.tr7c.r73, ptc15.tr7c.r85]}
#trackrecs = {ptc15.tr7c: [ptc15.tr7c.r71]}

# keep tracks in a fixed order:
tracknames = sorted([ track.absname for track in trackrecs ])
tracks = [ eval(trackname) for trackname in tracknames ]

ALPHA = 0.01 # p value threshold for significance
ONLYACTIVE = True # only consider active neurons (mean rate >= MINRATE global value)?
if ONLYACTIVE:
    print('ALPHA=%g, ONLYACTIVE=%r, MINRATE=%r' % (ALPHA, ONLYACTIVE, MINRATE))
else:
    print('ALPHA=%g, ONLYACTIVE=%r' % (ALPHA, ONLYACTIVE))
ec = 'grey'
allnids, totals = {}, [] # all nids, by track, >= MINRATE, or all if not ONLYACTIVE
allsnids, allthetas, allrs, alldepths, allps, allrates, allbestrec = {}, {}, {}, {}, {}, {}, {}
fs = fontsize() # save original font size
for track in tracks:
    tracknids = [] # nids for this track >= MINRATE, or all if not ONLYACTIVE
    # theta in deg, r in fraction of total spikes, depth in um:
    thetas, rs, depths, ps, rates, bestrec = {}, {}, {}, {}, {}, {}
    for rec in trackrecs[track]:
        if ONLYACTIVE:
            n = rec.n
        else:
            n = rec.alln
        # sorted list of relevant neurons for this recording:
        neurons = [ n[nid] for nid in np.array(sorted(n)) ]
        nids = [ neuron.id for neuron in neurons ] # corresponding nids
        tracknids.append(nids)
        snids = [] # significantly tuned nids for this rec
        for nid in nids:
            neuron = rec.alln[nid]
            tune = neuron.tune()
            theta, r, z, p = tune.pref(var='ori')
            if not p < ALPHA:
                continue # insignificant tuning, skip to next nid
            if nid in rs and r <= rs[nid]:
                continue # skip nid if it's less tuned than same nid from earlier rec
            snids.append(nid)
            thetas[nid] = theta
            rs[nid] = r
            depths[nid] = neuron.pos[1]
            ps[nid] = p
            rates[nid] = neuron.meanrate
            bestrec[nid] = rec.name
        print('%s: %d/%d neurons tuned' % (rec.absname, len(snids), len(nids)))
    nids = np.unique(np.hstack(tracknids)) # nids that were relevant in at least one rec
    snids = sorted(thetas.keys()) # significantly tuned nids across all recs
    thetas = np.asarray([ thetas[nid] for nid in snids ])
    rs = np.asarray([ rs[nid] for nid in snids ])
    depths = np.asarray([ depths[nid] for nid in snids ])
    ps = np.asarray([ ps[nid] for nid in snids ])
    rates = np.asarray([ rates[nid] for nid in snids ])
    cellmaxdepth = depths.max()
    chanmaxdepth = track.sort.chanpos[:, 1].max()
    maxdepth = max(cellmaxdepth, chanmaxdepth)
    depths /= maxdepth # normalize by cell or channel of greatest depth
    allnids[track.absname] = nids
    allsnids[track.absname] = snids
    allthetas[track.absname] = thetas
    allrs[track.absname] = rs
    alldepths[track.absname] = depths
    allps[track.absname] = ps
    allrates[track.absname] = rates
    allbestrec[track.absname] = bestrec

    # plot tuning strength vs theta, in half polar plot, colour by depth, with darker colours
    # indicating greater depth:
    f = figure(figsize=(7.5, 6.5))
    a = fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1.02),
                              thlabel=None, rlabel=None, ticklabels=False)
    a.scatter(thetas, rs, marker='o', edgecolors=ec, linewidth=1, s=175,
              cmap=cm.hot_r, c=depths, vmin=0, vmax=1, zorder=-0.1)
    # plot overlapping edges to help distinguish points:
    a.scatter(thetas, rs, marker='o', edgecolors=ec, facecolors='none', linewidth=1, s=175)
    #colorbar(sc)
    f.tight_layout(pad=0.3)
    f.canvas.manager.set_window_title(track.absname)

    # plot depth vs theta, and colour by r, with lighter colours indicating higher r:
    f = figure(figsize=(4.88, 3))
    fontsize(20)
    scatter(thetas, depths, marker='o', edgecolors=ec, linewidth=0.5, s=50,
            cmap=cm.hot, c=rs, vmin=0, vmax=1, zorder=-0.1)
    # plot overlapping edges to help distinguish points:
    scatter(thetas, depths, marker='o', edgecolors=ec, facecolors='none', linewidth=0.5, s=50)
    xlim(0, 180)
    xticks([0, 90, 180])
    ylim(1, 0)
    yticks([1, 0.5, 0])
    #xlabel('orientation preference ($^{\circ}$)')
    #ylabel('normalized depth')
    f.tight_layout(pad=0.3)
    f.canvas.manager.set_window_title(track.absname+'_depth')
    fontsize(fs) # restore
    if ONLYACTIVE:
        tracktotal = len(nids)
    else: # count all cells in each track, even if they didn't fire at all during oriented recs
        tracktotal = track.nallneurons
    totals.append(tracktotal)
    nsn, tn = len(snids), tracktotal
    print('%s: %d/%d (%g%%) neurons tuned' % (track.absname, nsn, tn, nsn/tn*100))

# create flattened versions:
snids = np.hstack([ allsnids[track.absname] for track in tracks ])
thetas = np.hstack([ allthetas[track.absname] for track in tracks ])
rs = np.hstack([ allrs[track.absname] for track in tracks ])
depths = np.hstack([ alldepths[track.absname] for track in tracks ])
ps = np.hstack([ allps[track.absname] for track in tracks ])
rates = np.hstack([ allrates[track.absname] for track in tracks ])

nsn, tn = len(snids), sum(totals)
print('total: %d/%d (%g%%) neurons tuned' % (nsn, tn, nsn/tn*100))

# plot tuning strength vs theta, in half polar plot, colour by depth, with darker colours
# indicating greater depth:
f = figure(figsize=(7.5, 6.5))
a = fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1.02),
                          thlabel='orientation preference', rlabel='tuning strength')
a.scatter(thetas, rs, marker='o', edgecolors=ec, linewidth=0.5, s=25,
          cmap=cm.hot_r, c=depths, vmin=0, vmax=1, zorder=-0.1)
# plot overlapping edges to help distinguish points:
a.scatter(thetas, rs, marker='o', edgecolors=ec, facecolors='none',
          linewidth=0.5, s=25)
#colorbar(sc)
f.tight_layout(pad=0.3)
f.canvas.manager.set_window_title('all_tracks')

# plot depth vs theta, and colour by r, with lighter colours indicating higher r:
f = figure(figsize=(7.5, 3.25))
scatter(thetas, depths, marker='o', edgecolors=ec, linewidth=0.5, s=25,
        cmap=cm.hot, c=rs, vmin=0, vmax=1, zorder=-0.1)
# plot overlapping edges to help distinguish points:
scatter(thetas, depths, marker='o', edgecolors=ec, facecolors='none', linewidth=0.5, s=25)
xlim(0, 180)
xticks([0, 90, 180])
ylim(1, 0)
yticks([1, 0.5, 0])
xlabel('orientation preference ($^{\circ}$)')
ylabel('normalized depth')
f.tight_layout(pad=0.3)
f.canvas.manager.set_window_title('all_tracks_depth')

# plot tuning strength vs mean firing rates during best recording:
f = figure(figsize=(3.5, 3.5))
scatter(rates, rs, marker='o', s=30, c='k', edgecolor='none', alpha=0.5)
xlim(1e-3, 1e2)
ylim(0, 1.02)
xscale('log')
xlabel('mean firing rate (Hz)')
ylabel('tuning strength')
# plot linear-log regression:
m, b, r, p, stderr = scipy.stats.linregress(np.log10(rates), rs)
rr = np.asarray(xlim()) # rates range
plot(rr, m*np.log10(rr)+b, 'r', ls='--', lw=4, alpha=0.7)
# add text in upper right corner:
text(0.75, 0.99, 'r=%.2f\np=%.1g' % (r, p),
     color='k', alpha=1, transform=gca().transAxes,
     horizontalalignment='left', verticalalignment='top')
f.tight_layout(pad=0.3)
f.canvas.manager.set_window_title('tuning_strength_vs_meanrates')

# plot distribution of mean firing rates of significantly tuned cells during best recording
figure(figsize=(3.5, 3.5))
logrates = np.log10(rates)
logmin, logmax = min(logrates), max(logrates)
logstart, logend = np.floor(logmin), np.ceil(logmax)
nbins = 15
edges = np.logspace(logstart, logend, nbins+1) # nbins+1 points in log space
hist(rates, bins=edges, color='k')
xscale('log')
xlabel('mean firing rate (Hz)')
ylabel('tuned neuron count')
ymin, ymax, lw = 18, 20, 2
# mark the distribution mean:
logmean = mean(log10(rates))
# arrow doesn't display correctly on log axis, use annotate instead:
#vlines(10**logmean, ymin, ymax, colors='k', lw=lw)
annotate('', xy=(10**logmean, 17), xycoords='data',
             xytext=(10**logmean, 20), textcoords='data',
             arrowprops=dict(fc='k', ec='none', width=1, headwidth=10, frac=0.5))
text(0.43, 0.97, '$\mu$=%.2f Hz' % 10**logmean,
                 horizontalalignment='right',
                 verticalalignment='top', linespacing=1.15,
                 transform=gca().transAxes,
                 color='k')
tight_layout(pad=0.3) # crop figure to contents
gcfm().window.setWindowTitle('tuned_meanrates_hist')

## TODO: why are some p-values negative??

show()


## TODO: add scatter plots of ori prefs of cells from driftbar vs flashed grating, one per
## track, or overplot scatter of all three tracks, giving each a different colour. Also might
## scatter plot tuning strengths of same cells. Mention delta t in hours among each pair of
## recordings

## TODO: plot tuning pref over time, coloured by strength? Would need to decide whether to
## include those cells that lost or gained significance?
