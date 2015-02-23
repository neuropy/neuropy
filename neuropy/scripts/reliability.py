"""Measure trial raster plot reliability, during natural scene movie stimulation as a function
of cortical state, by taking all pairwise correlations between all trials for a given neuron,
and then averaging them to get a single number between 0 and 1 representing reliability of
that trial raster plot. See Goard & Dan, 2009"""

from __future__ import division, print_function

from core import get_ssnids

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

recs = ptc22tr1r08s + ptc22tr1r10s
stranges = strangesr08s + strangesr10s

MINTRIALRATE = 0.2 # Hz, that's at least 1 spike per trial for 4.5 s trials
WEIGHTED = True
BINW, TRES = 0.02, 0.005
BLANK = False
NIDSKIND = 'all' # 'active' or 'all'
figsize = 3.5, 3.5

# get active or all neuron ids for each section of both r08 and r10:
ssnids, recsecnids = get_ssnids(recs, stranges, kind=NIDSKIND)

recsecrel = []
for rec, nids, strange in zip(recs, recsecnids, stranges):
    n2rel = {} # nid:reliability mapping for this recsec
    n2count, n2totcount, ts = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                                             binw=BINW, tres=TRES)
    for nid in nids:
        cs = n2count[nid] # trial counts
        totcs = n2totcount[nid] # total spike counts per trial
        trialwidth = ts[-1, 1] - ts[0, 0] # end of last bin minus start of first bin, in sec
        totalntrials = len(cs)
        # keep only those trials with sufficiently high spike rate:
        cs = cs[totcs/trialwidth >= MINTRIALRATE]
        ntrials = len(cs)
        if ntrials < 2:
            print("nid%d has insufficient trials that pass rate thresh" % nid)
            continue # not enough trials to correlate
        print("nid%d: %d --> %d trials after applying %g Hz trial rate thresh"
              % (nid, totalntrials, ntrials, MINTRIALRATE))
        rhos, weights = core.pairwisecorr(cs, weighted=WEIGHTED)
        #n2rel[nid] = np.mean(rhos)
        n2rel[nid] = np.median(rhos)
        #n2rel[nid] = (rhos*weights).sum()
    recsecrel.append(n2rel)

# 1st transition: desynched to synched
nids0 = np.intersect1d(sorted(recsecrel[0]), sorted(recsecrel[1]))
desynchrel0 = [ recsecrel[0][nid] for nid in nids0 ]
synchrel0 = [ recsecrel[1][nid] for nid in nids0 ]
# 2nd transition: synched to desynched
nids1 = np.intersect1d(sorted(recsecrel[2]), sorted(recsecrel[3]))
synchrel1 = [ recsecrel[2][nid] for nid in nids1 ]
desynchrel1 = [ recsecrel[3][nid] for nid in nids1 ]
# plot all transitions in one plot:
synchrel = synchrel0 + synchrel1
desynchrel = desynchrel0 + desynchrel1
figure(figsize=figsize)
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
# scatter plot reliability in the two states
plot(synchrel, desynchrel, 'k.', alpha=0.5, mec='None')
xlabel('synchronized reliability')
ylabel('desynchronized reliability')
xlim(min(-0.1, min(synchrel)), 1)
ylim(min(-0.1, min(desynchrel)), 1)
x0, x1, y0, y1 = axis()
gcfm().window.setWindowTitle('reliability')
tight_layout(pad=0.3)

show()

'''
# test single neuron in single recording section:
nid = 17
n2count, n2totcount, ts = ptc22.tr1.r08.bintraster(nids=[nid], blank=BLANK,
                                                   strange=strangesr08s[0],
                                                   binw=BINW, tres=TRES)
cs = n2count[nid] # binned trial counts, ntrials x nbins
totcs = n2totcount[nid] # total spike counts per trial
trialwidth = ts[-1, 1] - ts[0, 0] # end of last bin minus start of first bin, in sec
# filter out trials with too few spikes
totalntrials = len(cs)
cs = cs[totcs/trialwidth >= MINTRIALRATE]
ntrials = len(cs)
print("ntrials: %d --> %d after applying %g Hz trial rate thresh"
      % (totalntrials, ntrials, MINTRIALRATE))
rhos, weights = core.pairwisecorr(cs, weighted=WEIGHTED)

print("mean: %g, weighted mean: %g, median: %g"
      % (np.mean(rhos), (rhos*weights).sum(), np.median(rhos)))
figure()
hist(rhos, bins=100)
show()
'''
