"""Plot 2D matrices of natural scene movie PSTH correlations of all responsive cells in all 6
natural scene 5s movie clip recordings with cortical state changes in them. Also plot PSTH
correlation distributions, and correlations as a function of cell pair separation. Run from
within neuropy using `run -i scripts/psthcorr.py`"""

#do i just want to plot two distribs, or do i want to try scatter plots? probably both

from __future__ import division, print_function
import pylab as pl
import numpy as np
from scipy.stats import mannwhitneyu, ttest_1samp #, ttest_ind, ks_2samp

import core
from core import ceilsigfig, floorsigfig, scatterbin

from psth_funcs import get_psth_peaks_gac, get_seps

BINW, TRES = 0.02, 0.0001 # PSTH time bins, sec
GAUSS = True # calculate PSTH and single trial rates by convolving with Gaussian kernel?
BLANK = False # consider blank periods between trials?
KIND = 'active' # 'responsive' or 'active'
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz

# plotting params:
FIGSIZE = (3, 3)
PLOTRHOMATRICES = False
SHOWCOLORBAR = False # show colorbar for rho matrices?
SEPBINW = 200 # separation bin width, um
RHOMIN, RHOMAX = -0.4, 1
SEPMAX = 1375 # max pairwise separation, um

ALPHA = 0.05 # for comparing the means of psthcorr distribs to 0
VMIN, VMAX = -1, 1 # rho limits for correlation matrices


def get_nids_psths(rec, strange, kind='responsive'):
    """Return responsive or active nids and corresponding PSTHs for strange in rec"""
    if kind == 'responsive': # start with all nids, then weed out unresponsive ones
        nids = rec.get_nids(tranges=[strange], kind='all')
    else: # kind is 'active'
        nids = rec.get_nids(tranges=[strange], kind='active')
    t, psths, spikets = rec.psth(nids=nids, natexps=False, blank=BLANK, strange=strange,
                                 plot=False, binw=BINW, tres=TRES, gauss=GAUSS, norm='ntrials')
    if kind == 'responsive':
        rnids, rpsths = [], [] # nids and PSTHS to return
        for nid, psth, ts in zip(nids, psths, spikets):
            # run PSTH peak detection on this PSTH:
            baseline = MEDIANX * np.median(psth)
            thresh = baseline + MINTHRESH # peak detection threshold
            print("n%d" % nid, end='')
            peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh)
            npeaks = len(peakis)
            if npeaks == 0:
                continue # this PSTH has no peaks, this nid is unresponsive
            rnids.append(nid)
            rpsths.append(psth)
        print()
        return np.asarray(rnids), np.asarray(rpsths)
    else: # kind is 'active'
        return np.asarray(nids), np.asarray(psths)


if __name__ == "__main__":

    # mapping of recording to list of desynched and synched trange, in that order, copied from
    # psth_precision.py:
    rec2tranges = {ptc17.tr2b.r58: [(0, 700e6), # desynched trange, 66 Hz refresh rate
                                    (800e6, 1117e6)], # synched trange, 66 Hz refresh rate
                   ptc18.tr1.r38:  [(0, 425e6), # desynched trange, ends ~ trial 76
                                    (550e6, 2243e6)], # synched trange, starts ~ trial 98
                   ptc18.tr2c.r58: [(0, 750e6), # desynched trange
                                    (1000e6, 2248e6)], # synched trange
                   ptc22.tr1.r08:  [(0, 1500e6), # desynched trange
                                    (1550e6, 2329e6)], # synched trange
                   ptc22.tr1.r10:  [(1480e6, 2331e6), # desynched trange
                                    (0, 1400e6)], # synched trange
                   ptc22.tr4b.r49: [(0, 1475e6), # desynched trange
                                    (1500e6, 2331e6)], # synched trange
                  }
    # compare and sort recordings by their absname:
    reccmp = lambda reca, recb: cmp(reca.absname, recb.absname)
    urecs = sorted(rec2tranges, cmp=reccmp) # unique recordings, no repetition, sorted
    urecnames = ' '.join([rec.absname for rec in urecs])

    slabels = ['desynch', 'synch'] # state labels
    colours = ['b', 'r'] # corresponding state colours
    rhoslist = {'desynch': [], 'synch': []}
    sepslist = {'desynch': [], 'synch': []}
    for rec in urecs:
        stranges = rec2tranges[rec]
        for slabel, strange in zip(slabels, stranges):
            print()
            print(rec.absname, slabel, strange)
            nids, psths = get_nids_psths(rec, strange, kind=KIND)
            nn = len(nids)
            rho = np.corrcoef(psths) # rho matrix, defaults to bias=1
            rho[np.diag_indices(nn)] = np.nan # nan the diagonal, which imshow plots as white

            # collect rho values:
            lti = np.tril_indices(nn, -1) # lower triangle indices of rho matrix
            rhoslist[slabel].append(rho[lti])

            # collect corresponding pairwise neuron separation distances:
            sepslist[slabel].append(get_seps(nids, rec.alln))

            if PLOTRHOMATRICES:
                # plot rho matrix:
                figure(figsize=FIGSIZE)
                imshow(rho, vmin=VMIN, vmax=VMAX, cmap='jet') # cmap='gray' is too bland
                nidticks = np.arange(0, nn, 10)
                xticks(nidticks)
                yticks(nidticks)
                if SHOWCOLORBAR:
                    colorbar(ticks=[-1, 0, 1])
                titlestr = rec.absname + '_' + ('%s' % slabel) + '_rho_mat'
                gcfm().window.setWindowTitle(titlestr)
                tight_layout(pad=0.3)

    # concatenate rho and sep lists into arrays:
    rhos, seps = {}, {}
    for slabel in slabels:
        rhos[slabel] = np.hstack(rhoslist[slabel])
        seps[slabel] = np.hstack(sepslist[slabel])

    # plot rho histogram:
    dmean = rhos['desynch'].mean()
    smean = rhos['synch'].mean()
    u, p = mannwhitneyu(rhos['desynch'], rhos['synch']) # 1-sided
    if p < ALPHA:
        pstring = '$p<%g$' % ceilsigfig(p)
    else:
        pstring = '$p>%g$' % floorsigfig(p)
    '''
    # T-tests of both distribs relative to 0, not so useful:
    dt, dp = ttest_1samp(rhos['desynch'], 0) # 2-sided ttest relative to 0
    st, sp = ttest_1samp(rhos['synch'], 0) # 2-sided ttest relative to 0
    print('dmean=%g, t=%g, p=%g' % (dmean, dt, dp))
    print('smean=%g, t=%g, p=%g' % (smean, st, sp))
    if dp < ALPHA:
        dpstring = '$p<%g$' % ceilsigfig(dp)
    else:
        dpstring = '$p>%g$' % floorsigfig(dp)
    if sp < ALPHA:
        spstring = '$p<%g$' % ceilsigfig(sp)
    else:
        spstring = '$p>%g$' % floorsigfig(sp)
    '''
    figure(figsize=FIGSIZE)
    rhobins = np.arange(RHOMIN, RHOMAX+0.0333, 0.0333) # left edges + rightmost edge
    nd = hist(rhos['desynch'], bins=rhobins, histtype='step', color='b')[0]
    ns = hist(rhos['synch'], bins=rhobins, histtype='step', color='r')[0]
    nmax = max(nd.max(), ns.max())
    axvline(x=0, c='e', ls='-') # draw vertical grey line at x=0
    # draw arrows at means:
    arrow(dmean, 162, 0, -20, head_width=0.05, head_length=10, length_includes_head=True,
          color='b')
    arrow(smean, 162, 0, -20, head_width=0.05, head_length=10, length_includes_head=True,
          color='r')
    # draw vertical lines at means
    #axvline(x=dmean, c='b', ls='--')
    #axvline(x=smean, c='r', ls='--')
    xlim(xmin=RHOMIN, xmax=RHOMAX)
    ylim(ymax=nmax) # effectively normalizes the histogram
    # remove unnecessary decimal places:
    rhoticks = [-0.25, 0, 0.25, 0.5, 0.75, 1], ['-0.25', '0', '0.25', '0.5', '0.75', '1']
    #rhoticks = np.arange(-0.4, 1+0.2, 0.2)
    xticks(*rhoticks)
    yticks([0, nmax]) # turn off y ticks to save space
    xlabel(r'$\rho$')
    ylabel('cell pair count')
    text(0.98, 0.98, r'$\mu$=%.2g' % dmean, color='b',
         transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
    text(0.98, 0.90, r'$\mu$=%.2g' % smean, color='r',
         transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
    text(0.98, 0.82, '%s' % pstring, color='k',
         transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
    #text(0.98, 0.82, '%s' % dpstring, color='b',
    #     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
    #text(0.98, 0.74, '%s' % spstring, color='r',
    #     transform=gca().transAxes, horizontalalignment='right', verticalalignment='top')
    gcfm().window.setWindowTitle('rho_hist')
    tight_layout(pad=0.3)

    # plot rho vs separation:
    figure(figsize=FIGSIZE)
    #pl.plot(sepmeans, rhomeans, 'r.-', ms=10, lw=2)
    for slabel, c in zip(slabels, colours):
        # scatter plot:
        pl.plot(seps[slabel], rhos[slabel], c+'.', alpha=0.5, ms=2)
        # bin seps and plot mean rho in each bin:
        sepbins = np.arange(0, SEPMAX+SEPBINW, SEPBINW) # left edges
        midseps, rhomeans, rhostds = scatterbin(seps[slabel], rhos[slabel], sepbins,
                                                xaverage=None)
        errorbar(midseps, rhomeans, yerr=rhostds, fmt=c+'.-', ms=10, lw=2, zorder=9999)
    xlim(xmin=0, xmax=SEPMAX)
    ylim(ymin=RHOMIN, ymax=RHOMAX)
    septicks = np.arange(0, SEPMAX+100, 500)
    xticks(septicks)
    yticks(*rhoticks)
    xlabel(r'cell pair separation (${\mu}m$)')
    ylabel(r'$\rho$')
    gcfm().window.setWindowTitle('rho_sep')
    tight_layout(pad=0.3)

    # plot rho histograms for ptc22.tr1.r08 and ptc22.tr1.r10: same track, different movies:
    d3mean = rhoslist['desynch'][3].mean()
    s3mean = rhoslist['synch'][3].mean()
    d4mean = rhoslist['desynch'][4].mean()
    s4mean = rhoslist['synch'][4].mean()
    p3 = mannwhitneyu(rhoslist['desynch'][3], rhoslist['synch'][3])[1] # 1-sided
    p4 = mannwhitneyu(rhoslist['desynch'][4], rhoslist['synch'][4])[1] # 1-sided
    p34d = mannwhitneyu(rhoslist['desynch'][3], rhoslist['desynch'][4])[1] # 1-sided
    p34s = mannwhitneyu(rhoslist['synch'][3], rhoslist['synch'][4])[1] # 1-sided
    print('p3=%.2g, p4=%.2g, p34d=%.2g, p34s=%.2g' % (p3, p4, p34d, p34s))
    #if p < ALPHA:
    #    pstring = '$p<%g$' % ceilsigfig(p)
    #else:
    #    pstring = '$p>%g$' % floorsigfig(p)
    figure(figsize=FIGSIZE)
    rhobins = np.arange(RHOMIN, RHOMAX+0.0667, 0.0667) # left edges + rightmost edge
    # ptc22.tr1.r08 and ptc22.tr1.r10 are recording indices 3 and 4, respectively
    nd3 = hist(rhoslist['desynch'][3], bins=rhobins, histtype='step', color='b')[0]
    ns3 = hist(rhoslist['synch'][3], bins=rhobins, histtype='step', color='r')[0]
    nd4 = hist(rhoslist['desynch'][4], bins=rhobins, histtype='step', color='b', alpha=0.5)[0]
    ns4 = hist(rhoslist['synch'][4], bins=rhobins, histtype='step', color='r', alpha=0.5)[0]
    nmax = max(np.hstack([nd3, ns3, nd4, ns4]))
    axvline(x=0, c='e', ls='-') # draw vertical grey line at x=0
    # draw arrows at means:
    arrow(d3mean, nmax/2, 0, -10, head_width=0.05, head_length=5, length_includes_head=True,
          color='b')
    arrow(s3mean, nmax, 0, -10, head_width=0.05, head_length=5, length_includes_head=True,
          color='r')
    arrow(d4mean, nmax/2, 0, -10, head_width=0.05, head_length=5, length_includes_head=True,
          color='b', alpha=0.5)
    arrow(s4mean, nmax, 0, -10, head_width=0.05, head_length=5, length_includes_head=True,
          color='r', alpha=0.5)
    xlim(xmin=RHOMIN, xmax=RHOMAX)
    ylim(ymax=nmax) # effectively normalizes the histogram
    xticks(*rhoticks)
    yticks([0, nmax]) # turn off y ticks to save space
    xlabel(r'$\rho$')
    ylabel('cell pair count')
    gcfm().window.setWindowTitle('rho_hist_r08_r10')
    tight_layout(pad=0.3)

    show()
