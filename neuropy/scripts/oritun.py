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

# maybe include grating experiments as well, if necessary?
#recs = [ptc15.tr7c.r71, ptc22.tr1.r03, ptc22.tr1.r18, ptc22.tr2.r25, ptc22.tr2.r31]
recs = [ptc15.tr7c.r71, ptc22.tr1.r03, ptc22.tr2.r25]
ec = 'gray'
allthetas, allrs, alldepths = [], [], []
fs = fontsize() # save original font size
for rec in recs:
    thetas, rs, depths = [], [], [] # theta in deg, r in fraction of total spikes, depth in um
    nids = np.array(sorted(rec.alln))
    for nid in nids:
        neuron = rec.alln[nid]
        tune = neuron.tune()
        tune.calc()
        tune.plot(var='ori', plot=False)
        # don't allow oris > 180 deg, otherwise completely direction independent responses
        # will cancel out, resulting in no apparent tuning. Do this by angle doubling:
        oris, counts = tune.x % 180, tune.y # deg
        orisrad = 2 * oris * np.pi/180 # convert from deg to rad
        x = (counts*np.cos(orisrad)).sum()
        y = (counts*np.sin(orisrad)).sum()
        # arctan2 takes sign of x and y into account, then undo the angle doubling:
        theta = np.arctan2(y, x) / 2 # rad
        theta = theta * 180 / np.pi + 90 # rad to deg, off by 90 deg for some reason
        r = np.sqrt(x**2+y**2) / counts.sum() # fraction of total spikes
        thetas.append(theta)
        rs.append(r)
        depths.append(neuron.pos[1])
    thetas = np.asarray(thetas)
    rs = np.asarray(rs)
    chanmaxdepth = rec.sort.chanpos[:, 1].max()
    depths = np.asarray(depths)
    #depths = (chanmaxdepth - depths) / chanmaxdepth # invert and normalize
    depths /= chanmaxdepth # normalize
    alldepths.append(depths)
    allthetas.append(thetas)
    allrs.append(rs)

    # plot tuning strength vs theta, in half polar plot:
    f = figure()
    a = fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1.02),
                              thlabel=None, rlabel=None, ticklabels=False)
    a.scatter(thetas, rs, marker='o', edgecolors=ec, linewidth=1, s=175,
                   cmap=cm.hot_r, c=depths, zorder=-0.1)
    # plot overlapping edges to help distinguish points:
    a.scatter(thetas, rs, marker='o', edgecolors=ec, facecolors='none', linewidth=1, s=175)
    #colorbar(sc)
    f.tight_layout(pad=0.3)
    f.canvas.manager.set_window_title(rec.absname)
    f.show()

    # plot depth vs theta:
    f = figure(figsize=(4.88, 3))
    fontsize(20)
    scatter(thetas, depths, marker='o', edgecolors=ec, linewidth=0.5, s=50,
            cmap=cm.hot_r, c=depths, zorder=-0.1)
    # plot overlapping edges to help distinguish points:
    scatter(thetas, depths, marker='o', edgecolors=ec, facecolors='none', linewidth=0.5, s=50)
    xlim(0, 180)
    xticks([0, 90, 180])
    ylim(1, 0)
    yticks([1, 0.5, 0])
    #xlabel('orientation preference ($^{\circ}$)')
    #ylabel('normalized depth')
    f.tight_layout(pad=0.3)
    f.canvas.manager.set_window_title(rec.absname+'_depth')
    f.show()
    fontsize(fs) # restore

allthetas = np.hstack(allthetas)
allrs = np.hstack(allrs)
alldepths = np.hstack(alldepths)

af = figure()
aa = fractional_polar_axes(af, thlim=(0, 180), rlim=(0, 1.02),
                           thlabel='orientation preference', rlabel='tuning strength')
aa.scatter(allthetas, allrs, marker='o', edgecolors=ec, linewidth=0.5, s=25,
           cmap=cm.hot_r, c=alldepths, zorder=-0.1)
# plot overlapping edges to help distinguish points:
aa.scatter(allthetas, allrs, marker='o', edgecolors=ec, facecolors='none',
           linewidth=0.5, s=25)

#colorbar(sc)
af.tight_layout(pad=0.3)
af.canvas.manager.set_window_title('all tracks')
af.show()
