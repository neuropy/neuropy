"""Examine trial-averaged LFP and MUA during synched and desynched state of the specified
natural scene movies. Run from within neuropy using `run -i scripts/lfp_reliability.py`"""

from __future__ import division, print_function

import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy.stats import mannwhitneyu

from core import sparseness, ceilsigfig, corrcoef

FIGSIZE = (6, 3)
YLABELX = -0.06

# sort recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted
#urecnames = ' '.join([rec.absname for rec in urecs])
TLFP, TMUA = {}, {} # dicts storing results from rec.tlfp() and rec.tmua()
MAXMUA = {}
fmts = ('b-', 'r-') # desynched, then synched
LFPCORRS = [[], []] # desynched, then synched
MUACORRS = [[], []] # desynched, then synched
CORRBINW = 0.05
aw = 0.04 # arrow width

# calculate LFP and MUA time series, as well as correlations, one per trial:
for rec in urecs:
    print(rec.absname)
    TLFP[rec.absname] = [] # one entry per state
    TMUA[rec.absname] = []
    MAXMUA[rec.absname] = 0 # init, used for setting y limits during plotting
    stranges = REC2STATETRANGES[rec.absname]
    for statei, strange in enumerate(stranges): # desynched, then synched
        lfpt, lfptrials = rec.tlfp(trange=strange, plot=False)
        muat, muatrials = rec.tmua(trange=strange, plot=False) # Hz/unit
        TLFP[rec.absname].append((lfpt, lfptrials))
        TMUA[rec.absname].append((muat, muatrials))
        MAXMUA[rec.absname] = max(MAXMUA[rec.absname], muatrials.max())
        ntrials = len(lfptrials)
        assert ntrials == len(muatrials)
        for triali in range(ntrials):
            # measure reliability as correlation of each trial with mean of all others
            lfptrial, muatrial = lfptrials[triali], muatrials[triali]
            otheris = np.ones(ntrials, dtype=bool)
            otheris[triali] = False # exclude current trial
            LFPCORRS[statei].append(corrcoef(lfptrial, lfptrials[otheris].mean(axis=0)))
            MUACORRS[statei].append(corrcoef(muatrial, muatrials[otheris].mean(axis=0)))

for statei in range(2): # desynched, then synched
    LFPCORRS[statei] = np.asarray(LFPCORRS[statei]) # convert from list to array
    MUACORRS[statei] = np.asarray(MUACORRS[statei])

# plot LFP and MUA time series, plus mean and stdevs, and SNR time series:
for rec in urecs:
    print(rec.absname)
    # subplotting trickery from
    # http://stackoverflow.com/questions/22511550/gridspec-with-shared-axes-in-python
    lfpf = plt.figure(figsize=FIGSIZE)
    muaf = plt.figure(figsize=FIGSIZE)
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 0.5])
    LFPa1 = lfpf.add_subplot(gs[0])
    LFPa2 = lfpf.add_subplot(gs[1], sharex=LFPa1)
    LFPas = [LFPa1, LFPa2]
    plt.setp(LFPa1.get_xticklabels(), visible=False)
    plt.setp(LFPa2.get_xticklabels(), visible=True)
    MUAa1 = muaf.add_subplot(gs[0])
    MUAa2 = muaf.add_subplot(gs[1], sharex=MUAa1)
    MUAas = [MUAa1, MUAa2]
    plt.setp(MUAa1.get_xticklabels(), visible=False)
    plt.setp(MUAa2.get_xticklabels(), visible=True)
    # desynched, then synched:
    for statei, (fmt, LFPa, MUAa) in enumerate(zip(fmts, LFPas, MUAas)):
        lfpt, lfptrials = TLFP[rec.absname][statei]
        muat, muatrials = TMUA[rec.absname][statei]
        # normalize MUA across both states to get arbitrary units (0-1):
        nmuatrials = muatrials / MAXMUA[rec.absname] # normalized
        lfpmean, lfpstd = lfptrials.mean(axis=0), lfptrials.std(axis=0)
        nmuamean, nmuastd = nmuatrials.mean(axis=0), nmuatrials.std(axis=0)
        # to make saturation represent reliability, scale transparency inversely
        # with ntrials:
        alpha = 10 / ntrials
        # plot all trials for this rec, in mV
        # LFP:
        LFPa.plot(lfpt, lfptrials.T/1e3, fmt, alpha=alpha)
        LFPa.plot(lfpt, lfpmean/1e3, 'w-', alpha=1)
        LFPa.plot(lfpt, (lfpmean+lfpstd)/1e3, 'k-', alpha=1)
        LFPa.plot(lfpt, (lfpmean-lfpstd)/1e3, 'k-', alpha=1)
        LFPa.set_xlim(xmax=5.5)
        LFPa.set_ylim(-0.5, 0.5) # mV
        LFPa.set_yticks([-0.5, 0, 0.5])
        LFPa.set_yticklabels(['-0.5', '0', '0.5'])
        LFPa.set_ylabel("LFP (mV)")
        LFPa.yaxis.set_label_coords(YLABELX, 0.5)
        # MUA:
        MUAa.plot(muat, nmuatrials.T, fmt, alpha=alpha)
        MUAa.plot(muat, nmuamean, 'w-', alpha=1)
        MUAa.plot(muat, (nmuamean+nmuastd), 'k-', alpha=1)
        MUAa.plot(muat, (nmuamean-nmuastd), 'k-', alpha=1)
        MUAa.set_xlim(xmax=5.5)
        MUAa.set_ylim(0, 1) # AU
        MUAa.set_yticks([0, 1])
        MUAa.set_ylabel("MUA (AU)")
        MUAa.yaxis.set_label_coords(YLABELX, 0.5)
        print('sparseness(abs(mean(LFP))):', sparseness(np.abs(lfpmean)))
        print('sparseness(mean(MUA)):', sparseness(nmuamean))

    LFPa.set_xlabel("time (s)")
    lfpf.canvas.manager.set_window_title("LFP trials %s" % rec.absname)
    lfpf.tight_layout(pad=0.3) # crop figure to contents
    MUAa.set_xlabel("time (s)")
    muaf.canvas.manager.set_window_title("MUA trials %s" % rec.absname)
    muaf.tight_layout(pad=0.3) # crop figure to contents
    pl.show() # ensure figures pop up in order

# calculate significance
u, lfpp = mannwhitneyu(LFPCORRS[0], LFPCORRS[1]) # 1-sided
print('lfpp = %.2g' % lfpp)
u, muap = mannwhitneyu(MUACORRS[0], MUACORRS[1]) # 1-sided
print('muap = %.2g' % muap)

# plot LFPCORRS PDFs:
corrbins = np.arange(0, 1+CORRBINW, CORRBINW) # left edges + rightmost edge
nd = np.histogram(LFPCORRS[0], bins=corrbins, density=False)[0]
ns = np.histogram(LFPCORRS[1], bins=corrbins, density=False)[0]
dmax, smax = nd.max(), ns.max()
nmax = max(np.hstack([nd, ns]))
ah = nmax / 7 # arrow height
nd = np.hstack([nd, nd[-1]]) # repeat last point for right edge
ns = np.hstack([ns, ns[-1]])
figure(figsize=(3, 3))
plot(corrbins, nd, 'b-', lw=2)
plot(corrbins, ns, 'r-', lw=2)
# draw arrows at means:
dmean, smean = LFPCORRS[0].mean(), LFPCORRS[1].mean()
arrow(dmean, dmax+ah*1.1, 0, -ah, head_width=aw, head_length=ah/2,
      length_includes_head=True, color='b')
arrow(smean, smax+ah*1.1, 0, -ah, head_width=aw, head_length=ah/2,
      length_includes_head=True, color='r')
xlabel('LFP trial correlation')
ylabel('trial count')
xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '0.2', '0.4', '0.6', '0.8', '1'])
ylim(0, nmax+ah*1.1)
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.2f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.2f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(lfpp, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('LFP_trial_correlation_hist')
tight_layout(pad=0.3)
pl.show() # ensure figures pop up in order

# plot MUACORRS PDFs:
nd = np.histogram(MUACORRS[0], bins=corrbins, density=False)[0]
ns = np.histogram(MUACORRS[1], bins=corrbins, density=False)[0]
dmax, smax = nd.max(), ns.max()
nmax = max(np.hstack([nd, ns]))
ah = nmax / 7 # arrow height
nd = np.hstack([nd, nd[-1]]) # repeat last point for right edge
ns = np.hstack([ns, ns[-1]])
figure(figsize=(3, 3))
plot(corrbins, nd, 'b-', lw=2)
plot(corrbins, ns, 'r-', lw=2)
# draw arrows at means:
dmean, smean = MUACORRS[0].mean(), MUACORRS[1].mean()
arrow(dmean, dmax+ah*1.1, 0, -ah, head_width=aw, head_length=ah/2,
      length_includes_head=True, color='b')
arrow(smean, smax+ah*1.1, 0, -ah, head_width=aw, head_length=ah/2,
      length_includes_head=True, color='r')
xlabel('MUA trial correlation')
ylabel('trial count')
xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '0.2', '0.4', '0.6', '0.8', '1'])
ylim(0, nmax+ah*1.1)
# display means and p value:
text(0.02, 0.98, '$\mu$ = %.2f' % smean, # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.02, 0.90, '$\mu$ = %.2f' % dmean, # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.02, 0.82, 'p < %.1g' % ceilsigfig(muap, 1),
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('MUA_trial_correlation_hist')
tight_layout(pad=0.3)
pl.show() # ensure figures pop up in order

pl.show()
