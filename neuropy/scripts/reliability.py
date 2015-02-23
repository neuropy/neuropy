"""Measure trial raster plot reliability, during natural scene movie stimulation as a function
of cortical state, by taking all pairwise correlations between all trials for a given neuron,
and then averaging them to get a single number between 0 and 1 representing reliability of
that trial raster plot. See Goard & Dan, 2009"""

from __future__ import division, print_function

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

MINTRIALRATE = 0.2 # Hz, that's at least 1 spike per trial for 4.5 s trials
WEIGHTED = True
BINW, TRES = 0.02, 0.005
BLANK = False
NIDSKIND = 'all' # 'active' or 'all'

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
