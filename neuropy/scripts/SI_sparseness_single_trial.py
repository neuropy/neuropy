"""Plot SI vs sparseness for each trial for each unit in all state change natural movie
recordings. Should show single-trial positive correlation between the two, but the problem
with this is that any measure, whether precision, reliability or sparseness, is very noisy at
the single trial level"""

from __future__ import division, print_function

from core import sparseness, intround

from psth_funcs import get_psth_peaks_gac


LFPWIDTH, LFPTRES = 30, 5

TRASTERBINW, TRASTERTRES = 0.02, 0.001 # trial raster bins, sec
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
    allnids = sorted(rec.alln.keys())
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

spars = []
sis = []
for rec in urecs:
    print(rec.absname)
    si, sit = rec.lfp.si(kind='L/(L+H)', lfpwidth=LFPWIDTH, lfptres=LFPTRES, states=False,
                         relative2t0=False, lim2stim=False, plot=False)
    sit = intround(sit * 1e6) # convert from s to us for later comparison with strange
    #nids = sorted(rec.alln.keys())
    nids = get_responsive_nids(rec)
    print('nids:', nids)
    nn = len(nids)

    # calculate sparseness of PSTH calculated from sliding window of TRIALWINWIDTH
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
        # get PSTH for all nids over this strange:
        pstht, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK,
                                         strange=strange, plot=False, binw=BINW, tres=TRES,
                                         gauss=GAUSS, norm='ntrials')
        # slice out SI based on strange
        sit0i, sit1i = sit.searchsorted(strange)
        sitrials = si[sit0i:sit1i]
        sis.append(np.tile(sitrials.mean(), nn)) # SI is the same for each unit
        for psth in psths: # iterate over units
            spars.append(sparseness(psth))

sis = np.concatenate(sis)
spars = np.asarray(spars)

figure(figsize=figsize)
plot(sis, spars, 'k.', ms=0.5)
xlabel('trial range SI')
ylabel('trial range sparseness')
titlestr = 'SI_sparseness_single_trial'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

show()
