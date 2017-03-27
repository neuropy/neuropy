"""Plot SI vs sparseness for each trial for each unit in all state change natural movie
recordings. Should show single-trial positive correlation between the two, but the problem
with this is that any measure, whether precision, reliability or sparseness, is very noisy at
the single trial level"""

from core import sparseness, intround

from psth_funcs import get_psth_peaks_gac


LFPWIDTH, LFPTRES = 30, 5

TRASTERBINW, TRASTERTRES = 0.02, 0.001 # trial raster bins, sec
WEIGHT = False # weight trials by spike count for reliability measure?
BLANK = False
GAUSS = True
#BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
BINW, TRES = 0.02, 0.001 # PSTH time bins, sec
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
figsize = 5, 5 # inches
MINSPIKES = 5 # per trial, for calculating sparseness

TRIALWINWIDTH = 20 # number of trials to average over for each window of trials
TRIALWINTRES = 5 # trial offset between consecutive trial windows

def get_responsive_nids(rec):
    """Get responsive nids (those with at least 1 detected response event in their PSTH) by
    finding the responsive nids in each state separately, then taking their superset"""
    rnids = []
    allnids = sorted(rec.alln)
    for statei in [0, 1]: # 0 is desynched, 1 is synched
        tranges = REC2STATETRANGES[rec.absname]
        trange = tranges[statei]
        t, psths, spikets = rec.psth(nids=allnids, natexps=False, blank=BLANK, strange=trange,
                                     plot=False, binw=BINW, tres=TRES, gauss=GAUSS,
                                     norm='ntrials')
        for nid, psth, ts in zip(allnids, psths, spikets):
            # run PSTH peak detection:
            baseline = MEDIANX * np.median(psth)
            thresh = baseline + MINTHRESH # peak detection threshold
            print("n%d" % nid, end='')
            peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh)
            npeaks = len(peakis)
            if npeaks > 0:
                rnids.append(nid)
        print() # newline
    rnids = np.unique(rnids)
    return rnids

# sort state change recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted
urecnames = ' '.join([rec.absname for rec in urecs])

sis, rels, spars = [], [], []
for rec in urecs:
    print(rec.absname)
    si, sit = rec.lfp.si(kind='L/(L+H)', lfpwidth=LFPWIDTH, lfptres=LFPTRES, states=False,
                         relative2t0=False, lim2stim=False, plot=False)
    sit = intround(sit * 1e6) # convert from s to us for later comparison with strange
    #nids = sorted(rec.alln)
    nids = get_responsive_nids(rec)
    print('nids:', nids)
    nn = len(nids)

    # get reliability and sparseness of PSTHs calculated from sliding window of TRIALWINWIDTH
    # trials at a time, at tres of TRIALWINTRES trials. Represent SI for each range of trials
    # by the mean over those trials:
    ttranges = rec.trialtranges()[0]
    ntrials = len(ttranges)
    # build up trial indices of sliding window of trials:
    trialiranges = core.split_tranges([(0, ntrials-1)], TRIALWINWIDTH-1, TRIALWINTRES)
    for (trial0i, trial1i) in trialiranges: # iterate over all trial windows
        # get spike time range to consider for this window, from start of first trial to end
        # of last trial in this trial range:
        strange = ttranges[trial0i, 0], ttranges[trial1i, 1]
        # slice out SI based on strange, save mean SI for this trial range:
        sit0i, sit1i = sit.searchsorted(strange)
        sitrials = si[sit0i:sit1i]
        sis.append(np.tile(sitrials.mean(), nn)) # SI is the same for all units
        # calculate a new set of spike counts for all nids for this strange:
        n2count = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                                 binw=TRASTERBINW, tres=TRASTERTRES, gauss=GAUSS)[0]
        # get PSTH for all nids over this strange:
        pstht, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK,
                                         strange=strange, plot=False, binw=BINW, tres=TRES,
                                         gauss=GAUSS, norm='ntrials')
        # iterate over units:
        for nid, psth in zip(nids, psths):
            # calculate reliability:
            cs = n2count[nid] # 2D array of spike counts over trial time, one row per trial
            rhos, weights = core.pairwisecorr(cs, weight=WEIGHT, invalid='ignore')
            # set rho to 0 for trial pairs with undefined rho (one or both trials with 0
            # spikes):
            nanis = np.isnan(rhos)
            rhos[nanis] = 0.0
            rels.append(np.mean(rhos))
            # calculate sparseness:
            spars.append(sparseness(psth))
            

sis = np.concatenate(sis)
rels = np.asarray(rels)
spars = np.asarray(spars)

# plot reliability vs SI:
figure(figsize=figsize)
plot(sis, rels, 'k.', ms=0.5)
xlabel('trial range SI')
ylabel('trial range reliability')
titlestr = ('SI_reliability_trials_trialwinwidth=%d_trialwintres=%d'
            % (TRIALWINWIDTH, TRIALWINTRES))
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# plot sparseness vs SI:
figure(figsize=figsize)
plot(sis, spars, 'k.', ms=0.5)
xlabel('trial range SI')
ylabel('trial range sparseness')
titlestr = ('SI_sparseness_trials_trialwinwidth=%d_trialwintres=%d'
            % (TRIALWINWIDTH, TRIALWINTRES))
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

show()
