"""Plot SI vs sparseness for each trial for each unit in all state change natural movie
recordings. Should show single-trial positive correlation between the two, but the problem
with this is that any measure, whether precision, reliability or sparseness, is very noisy at
the single trial level"""

from __future__ import division, print_function

from core import sparseness

from psth_funcs import get_psth_peaks_gac


LFPWIDTH, LFPTRES = 30, 1

TRASTERBINW, TRASTERTRES = 0.02, 0.001 # trial raster bins, sec
BLANK = False
GAUSS = True
BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
figsize = 5, 5 # inches
MINSPIKES = 5 # per trial, for calculating sparseness

def get_responsive_nids(rec):
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
    rnids = np.unique(rnids)
    return rnids

# sort state change recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted
urecnames = ' '.join([rec.absname for rec in urecs])

spars = []
sis = []
for rec in urecs:
    si, t = rec.lfp.si(kind='L/(L+H)', lfpwidth=LFPWIDTH, lfptres=LFPTRES, states=False,
                       relative2t0=False, lim2stim=False, plot=False)
    #nids = sorted(rec.alln.keys())
    nids = get_responsive_nids(rec)
    nn = len(nids)
    n2count, n2totcount, bins, ttranges = rec.bintraster(nids=nids, blank=BLANK, strange=None,
        binw=TRASTERBINW, tres=TRASTERTRES, gauss=GAUSS)
    ttranges = ttranges / 1e6 # convert to s
    ntrials = len(ttranges)
    print('%d trials' % ntrials)

    # t are SI time bin midpoints, find all values of SI that fall within each trial trange,
    # average those SI values to get mean SI for each trial:
    trialis, trialsis = [], []
    for triali, ttrange in enumerate(ttranges):
        si0, si1 = t.searchsorted(ttrange)
        if si0 == si1:
            print(si0, si1)
            continue # SI undefined for this trial, skip to next trial
        trialis.append(triali)
        trialsis.append(si[si0:si1].mean())
    # SI as a function of trials for which it is defined, repeated for every neuron:
    sis.append(np.tile(trialsis, nn))

    # calculate sparseness for every single trial for every neuron:
    for nid in nids:
        count = n2count[nid]
        totcount = n2totcount[nid]
        for triali in trialis: # only iterate over those trials that have an SI value
            if totcount[triali] < MINSPIKES:
                # can't calculate sparseness reliably for single trials with too few spikes:
                spars.append(0)
            else:
                spars.append(sparseness(count[triali]))

sis = np.concatenate(sis)
spars = np.asarray(spars)

figure(figsize=figsize)
plot(sis, spars, 'k.', ms=0.5)
ylim(ymin=0.4)
xlabel('trial SI')
ylabel('trial sparseness')
titlestr = 'SI_sparseness_single_trial'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

show()
