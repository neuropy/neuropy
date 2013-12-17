"""Plot orientation preference vs cell position of all oriented cells from a number of drift
bar recordings. Run with 'run -i scripts/oripop.py' or copy and paste into neuropy console

TODO: use a colourmap, normalized from 0 to nn, indexed into by nidi, or from 0 to 1 and
indexed into by nidi/nn. This would show ori and tuning strength as a function of depth

TODO: for each unique nid in each track, go through all driftbar and grating recordings and
use the ori tuning from the earliest recording, or perhaps from the recording that provides
the strongest tuning. Need to do as much as possible to keep n as close to total n in the
track as possible

"""

from __future__ import division
from __future__ import print_function

from scripts.polar_demo import fractional_polar_axes

# maybe include grating experiments as well, if necessary?
#recs = [ptc15.tr7c.r71, ptc22.tr1.r03, ptc22.tr1.r18, ptc22.tr2.r25, ptc22.tr2.r31]
recs = [ptc15.tr7c.r71, ptc22.tr1.r03, ptc22.tr2.r25]
#recs = [ptc15.tr7c.r71]
colours = ['r', 'g', 'b']
allthetas, allrs = [], []

for rec, c in zip(recs, colours):
    thetas, rs = [], [] # theta in deg, r in fraction of total spikes
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
    thetas = np.asarray(thetas)
    rs = np.asarray(rs)
    allthetas.append(thetas)
    allrs.append(rs)

    f = figure()
    a = fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1.02),
                              thlabel='preferred orientation', rlabel='tuning strength')
    a.plot(thetas, rs, marker='.', ls='None', ms=7, c=c)
    f.tight_layout(pad=0.3)
    f.canvas.manager.set_window_title(rec.absname)
    f.show()

allthetas = np.hstack(allthetas)
allrs = np.hstack(allrs)

af = figure()
aa = fractional_polar_axes(af, thlim=(0, 180), rlim=(0, 1.02),
                           thlabel='preferred orientation', rlabel='tuning strength')
aa.plot(allthetas, allrs, marker='.', ls='None', ms=7, c='k')
af.tight_layout(pad=0.3)
af.canvas.manager.set_window_title('all tracks')
af.show()
