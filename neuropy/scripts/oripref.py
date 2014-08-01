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

alpha = 0.01 # p value threshold for significance
ec = 'gray'
allnids, allthetas, allrs, alldepths, allps, allrates = {}, {}, {}, {}, {}, {}
fs = fontsize() # save original font size
for track in tracks:
    # theta in deg, r in fraction of total spikes, depth in um:
    thetas, rs, depths, ps, rates = {}, {}, {}, {}, {}
    for rec in trackrecs[track]:
        nids = np.array(sorted(rec.alln))
        neurons = [ rec.alln[nid] for nid in nids ]
        for nid in nids:
            neuron = rec.alln[nid]
            tune = neuron.tune()
            theta, r, z, p = tune.pref(var='ori')
            if not p < alpha:
                continue # insignificant tuning, skip to next nid
            if nid in rs and r <= rs[nid]:
                continue # skip nid if it's less tuned than same nid from earlier rec
            thetas[nid] = theta
            rs[nid] = r
            depths[nid] = neuron.pos[1]
            ps[nid] = p
            rates[nid] = neuron.meanrate
        print('%s: %d of %d neurons tuned' % (rec.absname, len(thetas), rec.nallneurons))
    nids = sorted(thetas.keys()) # just the significant ones that were kept
    thetas = np.asarray([ thetas[nid] for nid in nids ])
    rs = np.asarray([ rs[nid] for nid in nids ])
    depths = np.asarray([ depths[nid] for nid in nids ])
    ps = np.asarray([ ps[nid] for nid in nids ])
    rates = np.asarray([ rates[nid] for nid in nids ])
    cellmaxdepth = depths.max()
    chanmaxdepth = track.sort.chanpos[:, 1].max()
    maxdepth = max(cellmaxdepth, chanmaxdepth)
    depths /= maxdepth # normalize by cell or channel of greatest depth
    allnids[track.absname] = nids
    allthetas[track.absname] = thetas
    allrs[track.absname] = rs
    alldepths[track.absname] = depths
    allps[track.absname] = ps
    allrates[track.absname] = rates

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
    print('%s: %d of %d neurons tuned' % (track.absname, len(nids), track.nallneurons))

# create flattened versions:
nids = np.hstack([ allnids[track.absname] for track in tracks ])
thetas = np.hstack([ allthetas[track.absname] for track in tracks ])
rs = np.hstack([ allrs[track.absname] for track in tracks ])
depths = np.hstack([ alldepths[track.absname] for track in tracks ])
ps = np.hstack([ allps[track.absname] for track in tracks ])
rates = np.hstack([ allrates[track.absname] for track in tracks ])

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
# plot linear log regression:
m, b, r, p, stderr = scipy.stats.linregress(np.log10(rates), rs)
rr = np.asarray(xlim()) # log rates range
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
