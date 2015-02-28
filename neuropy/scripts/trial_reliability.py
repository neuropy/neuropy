"""Measure trial raster plot reliability, during natural scene movie stimulation as a function
of cortical state, by taking all pairwise correlations between all trials for a given neuron,
and then averaging them to get a single number between 0 and 1 representing reliability of
that trial raster plot. See Goard & Dan, 2009"""

from __future__ import division, print_function

from scipy.stats import mannwhitneyu

from core import get_ssnids, ceilsigfig

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

recs = ptc22tr1r08s + ptc22tr1r10s
stranges = strangesr08s + strangesr10s

#MINTRIALRATE = 0.5 # Hz, 0.2 is at least 1 spike per trial for 4.5 s trials
#MINTRIALFRACTION = 0.25
WEIGHT = False
BINW, TRES = 0.02, 0.005
BLANK = False
NIDSKIND = 'all' # 'active' or 'all'
RELSTEP = 0.05
scatterfigsize = 3, 3
histfigsize = 3, 3

# get active or all neuron ids for each section of both r08 and r10:
ssnids, recsecnids = get_ssnids(recs, stranges, kind=NIDSKIND)

rels = []
for rec, nids, strange in zip(recs, recsecnids, stranges):
    n2rel = {} # nid:reliability mapping for this recsec
    n2count, n2totcount, ts = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                                             binw=BINW, tres=TRES)
    for nid in nids:
        cs = n2count[nid] # trial counts
        '''
        totcs = n2totcount[nid] # total spike counts per trial
        trialwidth = ts[-1, 1] - ts[0, 0] # end of last bin minus start of first bin, in sec
        totalntrials = len(cs)
        # keep only those trials with sufficiently high spike rate:
        cs = cs[totcs/trialwidth >= MINTRIALRATE]
        ntrials = len(cs)
        print("nid%d: %d --> %d trials after applying %g Hz trial rate thresh"
              % (nid, totalntrials, ntrials, MINTRIALRATE))
        # if MINTRIALRATE cuts down the number of trials by too much, skip this nid:
        trialfraction = ntrials / totalntrials
        if trialfraction < MINTRIALFRACTION:
            print("nid%d has insufficient trials that pass rate thresh" % nid)
            continue # not enough trials to correlate
        '''
        rhos, weights = core.pairwisecorr(cs, weight=WEIGHT, invalid='ignore')
        # replace any nans with 0s to represent corrs of pairs in which one or both trials
        # have no spikes:
        nanis = np.isnan(rhos)
        rhos[nanis] = 0.0
        n2rel[nid] = np.mean(rhos)
        #n2rel[nid] = np.median(rhos)
        #n2rel[nid] = (rhos*weights).sum()
    rels.append(n2rel)

# nids during at least one state on either side of 1st transition: desynched to synched
nids0 = np.union1d(sorted(rels[0]), sorted(rels[1]))
nids1 = np.union1d(sorted(rels[2]), sorted(rels[3]))

def filterdict(d, key, fallback=0):
    try:
        return d[key]
    except KeyError:
        return fallback

# 1st transition: desynched to synched
desynchrel0 = [ filterdict(rels[0], nid) for nid in nids0 ]
synchrel0 = [ filterdict(rels[1], nid) for nid in nids0 ]
# 2nd transition: synched to desynched
synchrel1 = [ filterdict(rels[2], nid) for nid in nids1 ]
desynchrel1 = [ filterdict(rels[3], nid) for nid in nids1 ]
# combine all transitions into one plot:
synchrel = synchrel0 + synchrel1
desynchrel = desynchrel0 + desynchrel1

# scatter plot reliability in the two states, note that some cells are discarded because
# they don't have enough trials with enough spikes in only one of the states (usually
# desynched):
figure(figsize=scatterfigsize)
plot([-1, 1], [-1, 1], 'e--') # plot y=x line
plot(synchrel, desynchrel, 'o', mec='k', mfc='None')
xlabel('synchronized trial reliability')
ylabel('desynchronized trial reliability')
xmin, xmax = min(-0.1, min(synchrel)), max(0.8, max(synchrel))
ymin, ymax = min(-0.1, min(desynchrel)), max(0.8, max(desynchrel))
xlim(xmin, xmax)
ylim(ymin, ymax)
xticks(np.arange(0, xlim()[1], 0.2))
yticks(np.arange(0, ylim()[1], 0.2))
titlestr = 'trial reliability ptc22.tr1.r08 ptc22.tr1.r10'
gcfm().window.setWindowTitle(titlestr)
#gcfm().window.setWindowTitle(titlestr+', MINTRIALRATE=%g, MINTRIALFRACTION=%g'
#                             % (MINTRIALRATE, MINTRIALFRACTION))
tight_layout(pad=0.3)

# plot distributions of trial reliability in the two states, gives higher N:
figure(figsize=histfigsize)
bins = np.arange(0, 1+RELSTEP, RELSTEP)
synchrel = rels[1].values() + rels[2].values()
desynchrel = rels[0].values() + rels[3].values()
n1 = hist(synchrel, bins=bins, color='r')[0] # synched
n0 = hist(desynchrel, bins=bins, color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=xmax)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(np.arange(0, xmax, 0.2))
xlabel('trial reliability')
ylabel('cell count')
#t, p = ttest_ind(desynchrel, synchrel, equal_var=False) # Welch's T-test
u, p = mannwhitneyu(desynchrel, synchrel) # 1-sided
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.2f' % mean(synchrel), # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.2f' % mean(desynchrel), # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'trial reliability hist ptc22.tr1.r08 ptc22.tr1.r10'
gcfm().window.setWindowTitle(titlestr)
#gcfm().window.setWindowTitle(titlestr+', MINTRIALRATE=%g, MINTRIALFRACTION=%g'
#                             % (MINTRIALRATE, MINTRIALFRACTION))
tight_layout(pad=0.3)

show()
