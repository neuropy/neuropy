"""Calculate MUA coupling and its effect on natscene response reliability and precision. Run
from within neuropy using `run -i scripts/MUA_coupling.py`"""

from __future__ import division, print_function

from scipy.stats import mannwhitneyu

import core
from core import sparseness, intround, ceilsigfig, floorsigfig

from psth_funcs import get_nids_psths

BLANK = False # consider blank periods between trials?
WEIGHT = False # weight trials by spike count for reliability measure?
BINW, TRES = 0.02, 0.001 # PSTH and MUA time bins, sec
BINWMS = '%dms' % intround(BINW * 1000)
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
if GAUSS:
    KERNEL = 'gauss'
else:
    KERNEL = 'square'
KIND = 'responsive' # which type of neurons to use? 'responsive' or 'active'
MEDIANX = 2 # PSTH median multiplier, Hz
MINTHRESH = 3 # peak detection thresh, Hz

# plotting params
#LOGNULLREL = -3
#NULLREL = 10**LOGNULLREL
FIGSIZE = 3, 3
COUPMIN, COUPMAX = -0.25, 1

# sort recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted
nrecs = len(urecs)

stateis = [0, 1] # desynched, then synched
slabels = ['desynch', 'synch'] # state labels
colours = ['b', 'r'] # corresponding state colours

rels, spars, coups = [[], []], [[], []], [[], []]
nreplacedbynullrel = 0
for rec in urecs:
    stranges = REC2STATETRANGES[rec.absname]
    for statei, strange in zip(stateis, stranges):
        print(rec.absname, slabels[statei], strange)
        nids, psths = get_nids_psths(rec, strange, kind=KIND, blank=BLANK,
                                     binw=BINW, tres=TRES, gauss=GAUSS,
                                     medianx=MEDIANX, minthresh=MINTHRESH)
        # n2count is needed for calculating reliability:
        n2count = rec.bintraster(nids=nids, blank=BLANK, strange=strange,
                                 binw=BINW, tres=TRES, gauss=GAUSS)[0]
        ## tmua calc should really be inside the
        ## loop below, and exclude the nid it's being correlated against
        ## That will be slow, need to make a new tmua mean f'n that collapses
        ## spikes *before* binning and convolving with gauss
        tmuas = rec.tmuas(width=BINW, tres=TRES, gauss=GAUSS, blank=BLANK,
                          trange=strange, plot=False)[1] # Hz/unit
        tmua = tmuas.mean(axis=0) # average across trials
        for nid, psth in zip(nids, psths):
            assert len(psth) == len(tmua)
            # calculate reliability of this PSTH:
            cs = n2count[nid] # 2D array of spike counts over trial time, one row per trial
            rhos, weights = core.pairwisecorr(cs, weight=WEIGHT, invalid='ignore')
            # set rho to 0 for trial pairs with undefined rho (one or both trials
            # with 0 spikes):
            nanis = np.isnan(rhos)
            rhos[nanis] = 0.0
            # for log plotting convenience, replace any mean rhos < NULLREL with NULLREL
            rel = np.mean(rhos)
            if rel < NULLREL:
                rel = NULLREL
                nreplacedbynullrel += 1
            rels[statei].append(rel)
            # calculate sparseness of this PSTH:
            spars[statei].append(sparseness(psth))
            # calculate coupling of this PSTH with tMUA:
            coup = core.corrcoef(psth, tmua)
            coups[statei].append(coup)
        print('\n') # two newlines

for statei in stateis:
    rels[statei] = np.asarray(rels[statei])
    spars[statei] = np.asarray(spars[statei])
    coups[statei] = np.asarray(coups[statei])

# plot MUA coupling histogram:
dmean = coups[0].mean()
smean = coups[1].mean()
u, p = mannwhitneyu(coups[0], coups[1]) # 1-sided
if p < ALPHA:
    pstring = 'p < %g' % ceilsigfig(p)
else:
    pstring = 'p > %g' % floorsigfig(p)
figure(figsize=FIGSIZE)
coupbins = np.arange(COUPMIN, COUPMAX+0.1, 0.1) # left edges + rightmost edge
nd = hist(coups[0], bins=coupbins, histtype='step', color='b')[0]
ns = hist(coups[1], bins=coupbins, histtype='step', color='r')[0]
nmax = max(np.hstack([nd, ns]))
axvline(x=0, c='e', ls='-', alpha=0.5, zorder=-1) # draw vertical grey line at x=0
# draw arrows at means:
ah = nmax / 8 # arrow height
arrow(dmean, nmax, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='b')
arrow(smean, nmax, 0, -ah, head_width=0.05, head_length=ah/2, length_includes_head=True,
      color='r')
xlim(xmin=COUPMIN, xmax=COUPMAX)
ylim(ymax=nmax*1.01)
# remove unnecessary decimal places:
#coupticks = [-0.25, 0, 0.25, 0.5, 0.75, 1], ['-0.25', '0', '0.25', '0.5', '0.75', '1']
#xticks(*coupticks)
yticks([0, nmax]) # turn off y ticks to save space
xlabel('MUA coupling')
ylabel('unit pair count')
text(0.98, 0.98, r'$\mu$ = %.2g' % smean, color='r',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
text(0.98, 0.90, r'$\mu$ = %.2g' % dmean, color='b',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
text(0.98, 0.82, '%s' % pstring, color='k',
     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
titlestr = 'MUA_coupling_hist'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)
