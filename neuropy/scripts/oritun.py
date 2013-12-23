"""Plot orientation preference vs tuning strength and cell position of all oriented cells from
a number of orientated stimulus recordings. Run with 'run -i scripts/oripop.py' or copy and
paste into neuropy console

TODO: for each unique nid in each track, go through all driftbar and drift/flash grating
recordings and use the ori tuning from the earliest recording, or perhaps from the recording
that provides the strongest tuning. Need to do as much as possible to keep n as close to total
n in the track as possible

"""

from __future__ import division
from __future__ import print_function

from scripts.polar_demo import fractional_polar_axes

# maybe include drifting and flashed grating experiments as well, if necessary?
#recs = [ptc15.tr7c.r71, ptc22.tr1.r03, ptc22.tr1.r18, ptc22.tr2.r25, ptc22.tr2.r31]
#recs = [ptc15.tr7c.r71, ptc22.tr1.r03, ptc22.tr2.r25]
#recs = [ptc15.tr7c.r71]
#recs = [ptc22.tr1.r03, ptc22.tr1.r18]
trackrecs = {ptc15.tr7c: [ptc15.tr7c.r71],
             ptc22.tr1: [ptc22.tr1.r03, ptc22.tr1.r18],
             ptc22.tr2: [ptc22.tr2.r25, ptc22.tr2.r31]}
alpha = 0.01 # p value threshold for significance
ec = 'gray'
allthetas, allrs, alldepths = [], [], []
fs = fontsize() # save original font size
for track, recs in trackrecs.items():
    thetas, rs, depths = {}, {}, {} # theta in deg, r in fraction of total spikes, depth in um
    for rec in recs:
        nids = np.array(sorted(rec.alln))
        for nid in nids:
            neuron = rec.alln[nid]
            tune = neuron.tune()
            theta, r, z, p = tune.pref(var='ori')
            if not p < alpha:
                continue # insignificant tuning, skip to next nid
            if nid in rs and r <= rs[nid]:
                continue # skip nid if it's less tuned than same nid from earlier rec
            thetas[nid] = theta # off by 90 deg for some reason
            rs[nid] = r
            depths[nid] = neuron.pos[1]
    nids = sorted(thetas.keys()) # just the significant ones that were kept
    thetas = np.asarray([ thetas[nid] for nid in nids ])
    rs = np.asarray([ rs[nid] for nid in nids ])
    depths = np.asarray([ depths[nid] for nid in nids ])
    depths /= depths.max() # normalize by cell of greatest depth
    #chanmaxdepth = track.sort.chanpos[:, 1].max()
    #depths /= chanmaxdepth # normalize by channel of maximum depth
    alldepths.append(depths)
    allthetas.append(thetas)
    allrs.append(rs)

    # plot tuning strength vs theta, in half polar plot, colour by depth, with darker colours
    # indicating greater depth:
    f = figure()
    a = fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1.02),
                              thlabel=None, rlabel=None, ticklabels=False)
    a.scatter(thetas, rs, marker='o', edgecolors=ec, linewidth=1, s=175,
              cmap=cm.hot_r, c=depths, vmin=0, vmax=1, zorder=-0.1)
    # plot overlapping edges to help distinguish points:
    a.scatter(thetas, rs, marker='o', edgecolors=ec, facecolors='none', linewidth=1, s=175)
    #colorbar(sc)
    f.tight_layout(pad=0.3)
    f.canvas.manager.set_window_title(track.absname)
    f.show()

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
    f.show()
    fontsize(fs) # restore
    print()
    print('%s: %d of %d neurons' % (track.absname, len(nids), track.nallneurons))

allthetas = np.hstack(allthetas)
allrs = np.hstack(allrs)
alldepths = np.hstack(alldepths)

af = figure()
aa = fractional_polar_axes(af, thlim=(0, 180), rlim=(0, 1.02),
                           thlabel='orientation preference', rlabel='tuning strength')
aa.scatter(allthetas, allrs, marker='o', edgecolors=ec, linewidth=0.5, s=25,
           cmap=cm.hot_r, c=alldepths, vmin=0, vmax=1, zorder=-0.1)
# plot overlapping edges to help distinguish points:
aa.scatter(allthetas, allrs, marker='o', edgecolors=ec, facecolors='none',
           linewidth=0.5, s=25)

#colorbar(sc)
af.tight_layout(pad=0.3)
af.canvas.manager.set_window_title('all tracks')
af.show()
