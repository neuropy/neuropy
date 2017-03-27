"""Calculate average cross-correlation histogram for all possible pairs of spike trains in a
given recording"""

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

nids = np.sort(list(rec.n))
"""
Shuffle the nids, which are sorted by depth - this allows examination of whether assymmetry in
the resulting mean CCH is a result of some kind of causal regularity as a function of cell
depth. To do this properly, when shuffling, should really generate lots of mean CCHs with
shuffled nids and average them. That would guarantee no assymmetry in the super averaged CCH:
"""
np.random.shuffle(nids)

nn = len(nids)
binw = 2000 # us
trange = np.array([-100000, 100000]) # us
bins = np.arange(trange[0], trange[1]+binw, binw)
hists = []
for nii0 in range(nn):
    for nii1 in range(nii0+1, nn):
        spikes0 = rec.n[nids[nii0]].spikes
        spikes1 = rec.n[nids[nii1]].spikes
        dts = util.xcorr(spikes0, spikes1, trange) # spike time differences in us
        hist = np.histogram(dts, bins=bins)[0]
        # if we don't normalize, we treat our confidence in the CCHs of cell pairs
        # proportional to the firing rates of cell pairs, which may be the optimal thing to do:
        #hist = hist / hist.sum() # pmf: normalize so that sum of each hist is 1
        hists.append(hist)
#hists = np.vstack(hists).sum(axis=0)
hists = np.vstack(hists).mean(axis=0)

figure()
bar(bins[:-1]/1000, hists, width=binw/1000)
xlim(trange/1000)
xlabel('time (ms)')
ylabel('mean cross-correlation histogram across all pairs')
title(rec.name)
tight_layout(pad=0.3)
