"""Examine trial-averaged LFP and MUA during synched and desynched state of the specified
natural scene movies. Run from within neuropy using `run -i scripts/lfp_reliability.py`"""

from __future__ import division, print_function

import matplotlib.pyplot as plt
from matplotlib import gridspec

from scipy.stats import mannwhitneyu

from core import sparseness, ceilsigfig

FIGSIZE = (8, 4)
YLABELX = -0.06

# sort recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted
#urecnames = ' '.join([rec.absname for rec in urecs])
fmts = ('b-', 'r-') # desynched and synched
LFPSNRs = [[], []] # desynched and synched
MUASNRs = [[], []]

TLFP, TMUA = {}, {} # dicts storing results from rec.tlfp() and rec.tmua()
MAXMUA = {}

# calculate LFP and MUA time series, one per trial:
for rec in urecs:
    print(rec.absname)
    TLFP[rec.absname] = [] # one entry per state
    TMUA[rec.absname] = []
    MAXMUA[rec.absname] = 0 # init, used for setting y limits during plotting
    stranges = REC2STATETRANGES[rec.absname]
    for strange in stranges: # desynched, then synched
        lfpt, lfptrials = rec.tlfp(trange=strange, plot=False)
        muat, muatrials = rec.tmua(trange=strange, plot=False) # Hz/unit
        TLFP[rec.absname].append((lfpt, lfptrials))
        TMUA[rec.absname].append((muat, muatrials))
        MAXMUA[rec.absname] = max(MAXMUA[rec.absname], muatrials.max())

# plot LFP and MUA time series, plus mean and stdevs, and SNR time series:
for rec in urecs:
    print(rec.absname)
    # subplotting trickery from
    # http://stackoverflow.com/questions/22511550/gridspec-with-shared-axes-in-python
    lfpf = plt.figure(figsize=FIGSIZE)
    muaf = plt.figure(figsize=FIGSIZE)
    gs = gridspec.GridSpec(3, 1, height_ratios=[0.3, 0.3, 0.4])
    LFPa1 = lfpf.add_subplot(gs[0])
    LFPa2 = lfpf.add_subplot(gs[1], sharex=LFPa1)
    LFPas = [LFPa1, LFPa2]
    LFPSNRa = lfpf.add_subplot(gs[2], sharex=LFPa1)
    plt.setp(LFPa1.get_xticklabels(), visible=False)
    plt.setp(LFPa2.get_xticklabels(), visible=False)
    MUAa1 = muaf.add_subplot(gs[0])
    MUAa2 = muaf.add_subplot(gs[1], sharex=MUAa1)
    MUAas = [MUAa1, MUAa2]
    MUASNRa = muaf.add_subplot(gs[2], sharex=MUAa1)
    plt.setp(MUAa1.get_xticklabels(), visible=False)
    plt.setp(MUAa2.get_xticklabels(), visible=False)
    # desynched, then synched:
    for statei, (fmt, LFPa, MUAa) in enumerate(zip(fmts, LFPas, MUAas)):
        lfpt, lfptrials = TLFP[rec.absname][statei]
        muat, muatrials = TMUA[rec.absname][statei]
        lfpmean, lfpstd = lfptrials.mean(axis=0), lfptrials.std(axis=0)
        muamean, muastd = muatrials.mean(axis=0), muatrials.std(axis=0)
        ntrials = len(lfptrials)
        assert ntrials == len(muatrials)
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
        LFPa.set_ylabel("LFP (mV)")
        LFPa.yaxis.set_label_coords(YLABELX, 0.5)
        # Fano-factor and CV don't work very well when mean approaches zero
        ## TODO: for reliability measure, exclude last sec of blankscreen in each trial:
        LFPSNR = np.abs(lfpmean) / lfpstd # something like S/N ratio
        LFPSNRs[statei].append(LFPSNR)
        #LFPSNR = lfpmean**2 / lfpstd**2 # something like (S/N)**2 ratio
        LFPSNRa.plot(lfpt, LFPSNR, fmt)
        # MUA:
        MUAa.plot(muat, muatrials.T, fmt, alpha=alpha)
        MUAa.plot(muat, muamean, 'w-', alpha=1)
        MUAa.plot(muat, (muamean+muastd), 'k-', alpha=1)
        MUAa.plot(muat, (muamean-muastd), 'k-', alpha=1)
        MUAa.set_xlim(xmax=5.5)
        maxmua = intround(MAXMUA[rec.absname])
        MUAa.set_ylim(0, maxmua) # Hz/unit
        MUAa.set_yticks([0, maxmua])
        MUAa.set_ylabel("MUA (Hz/unit)")
        MUAa.yaxis.set_label_coords(YLABELX, 0.5)
        # Fano-factor and CV don't work very well when mean approaches zero
        ## TODO: for reliability measure, exclude last sec of blankscreen in each trial:
        MUASNR = np.abs(muamean) / muastd # something like S/N ratio
        MUASNRs[statei].append(MUASNR)
        #MUASNR = muamean**2 / muastd**2 # something like (S/N)**2 ratio
        MUASNRa.plot(muat, MUASNR, fmt)
        print('LFP sparseness:', sparseness(np.abs(lfpmean)))
        print('mean tLFP S/N:', LFPSNR.mean())
        print('MUA sparseness:', sparseness(np.abs(muamean)))
        print('mean tMUA S/N:', MUASNR.mean())
    LFPSNRa.set_xlim(xmax=5.5)
    LFPSNRa.set_ylim(0, 3)
    LFPSNRa.set_yticks([0, 1, 2, 3])
    LFPSNRa.set_xlabel("time (s)")
    LFPSNRa.set_ylabel(r"LFP S/N ($|\mu|/\sigma)$")
    LFPSNRa.yaxis.set_label_coords(YLABELX+0.0012, 0.5) # has a descender
    lfpf.canvas.manager.set_window_title("LFP reliability %s" % rec.absname)
    lfpf.tight_layout(pad=0.3) # crop figure to contents
    MUASNRa.set_xlim(xmax=5.5)
    MUASNRa.set_ylim(0, 5)
    MUASNRa.set_yticks([0, 5])
    MUASNRa.set_xlabel("time (s)")
    MUASNRa.set_ylabel(r"MUA S/N ($|\mu|/\sigma)$")
    MUASNRa.yaxis.set_label_coords(YLABELX+0.0012, 0.5) # has a descender
    muaf.canvas.manager.set_window_title("MUA reliability %s" % rec.absname)
    muaf.tight_layout(pad=0.3) # crop figure to contents

    pl.show() # ensures figures pop up in order

# collapse SNRs across recs:
for statei in range(2):
    LFPSNRs[statei] = np.hstack(LFPSNRs[statei])
    MUASNRs[statei] = np.hstack(MUASNRs[statei])

# calculate significance
u, lfpp = mannwhitneyu(LFPSNRs[0], LFPSNRs[1]) # 1-sided
#if lfpp < np.finfo(np.float64).eps: # prevent 0 p-value due to underflow
#    lfpp = np.finfo(np.float64).eps
print('lfpp = %.2g' % lfpp)
u, muap = mannwhitneyu(MUASNRs[0], MUASNRs[1]) # 1-sided
#if muap < np.finfo(np.float64).eps: # prevent 0 p-value due to underflow
#    muap = np.finfo(np.float64).eps
print('muap = %.2g' % muap)

# calculate LFP PDFs:
SNRBINW = 0.1
snrbins = np.arange(0, 3+SNRBINW, SNRBINW) # left edges + rightmost edge
nd = np.histogram(LFPSNRs[0], bins=snrbins, density=True)[0]
ns = np.histogram(LFPSNRs[1], bins=snrbins, density=True)[0]
#nmax = max(np.hstack([nd, ns]))
# normalize density to arbitrary units:
#nd = nd / nmax
#ns = ns / nmax
nd = np.hstack([nd, nd[-1]]) # repeat last point for right edge
ns = np.hstack([ns, ns[-1]])
figure(figsize=(3, 3))
plot(snrbins, nd, 'b-', lw=2)
plot(snrbins, ns, 'r-', lw=2)
# draw arrows at means:
ah = 0.4 # arrow height
dmean, smean = LFPSNRs[0].mean(), LFPSNRs[1].mean()
arrow(dmean, 3.5, 0, -ah, head_width=0.1, head_length=ah/2,
      length_includes_head=True, color='b')
arrow(smean, 3.5, 0, -ah, head_width=0.1, head_length=ah/2,
      length_includes_head=True, color='r')
xlabel(r"LFP S/N ($|\mu|/\sigma)$")
ylabel('probability density')
xticks([0, 1, 2, 3])
yticks([0, 1, 2, 3])
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(lfpp, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('LFP_SNR_hist')
tight_layout(pad=0.3)

# calculate MUA PDFs:
SNRBINW = 0.1
snrbins = np.arange(0, 5+SNRBINW, SNRBINW) # left edges + rightmost edge
nd = np.histogram(MUASNRs[0], bins=snrbins, density=True)[0]
ns = np.histogram(MUASNRs[1], bins=snrbins, density=True)[0]
#nmax = max(np.hstack([nd, ns]))
# normalize density to arbitrary units:
#nd = nd / nmax
#ns = ns / nmax
nd = np.hstack([nd, nd[-1]]) # repeat last point for right edge
ns = np.hstack([ns, ns[-1]])
figure(figsize=(3, 3))
plot(snrbins, nd, 'b-', lw=2)
plot(snrbins, ns, 'r-', lw=2)
# draw arrows at means:
ah = 0.4 # arrow height
dmean, smean = MUASNRs[0].mean(), MUASNRs[1].mean()
arrow(dmean, 1, 0, -ah, head_width=0.3, head_length=ah/2,
      length_includes_head=True, color='b')
arrow(smean, 1, 0, -ah, head_width=0.3, head_length=ah/2,
      length_includes_head=True, color='r')
xlabel(r"MUA S/N ($|\mu|/\sigma)$")
ylabel('probability density')
#xticks([0, 1, 2, 3])
#yticks([0, 1, 2, 3])
# display means and p value:
text(0.98, 0.98, '$\mu$ = %.1f' % smean, # synched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.98, 0.90, '$\mu$ = %.1f' % dmean, # desynched
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.98, 0.82, 'p < %.1g' % ceilsigfig(muap, 1),
                 horizontalalignment='right', verticalalignment='top',
                 transform=gca().transAxes, color='k')
gcfm().window.setWindowTitle('MUA_SNR_hist')
tight_layout(pad=0.3)

pl.show()
