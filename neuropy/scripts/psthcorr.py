"""Plot 2D matrix of natural scene movie PSTH correlations of all active cells. Run from
within neuropy using `run -i scripts/psthcorr.py`"""

from __future__ import division
import pylab as pl
from scipy.stats import ttest_1samp, ttest_ind, ks_2samp

import core
from core import get_ssnids

FIGSIZE = (3, 3)
SHOWCOLORBAR = False # show colorbar
SEPBINW = 200 # separation bin width, um
RHOMIN, RHOMAX = -0.3, 1
RHODIFFMIN, RHODIFFMAX = -0.65, 0.65 # symmetric about 0 for the delta rhos

ptc15tr7crecs = [ptc15.tr7c.r74, ptc15.tr7c.r95b]
nateids = [3, 4, 10, 12] # for both recs in ptc15.tr7c
etrangesr74 = [ ptc15.tr7c.r74.e[nateid].trange for nateid in nateids ] # us
etrangesr95b = [ ptc15.tr7c.r74.e[nateid].trange for nateid in nateids ] # us

# copied to psth_precision.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

ptc22tr1recs = [ptc22.tr1.r05, ptc22.tr1.r08, ptc22.tr1.r10, ptc22.tr1.r19]
ptc22tr2recs = [ptc22.tr2.r28, ptc22.tr2.r33]

celltype2int = {'fast':0, 'slow':1, 'fastasym':2, 'slowasym':3,
                'simple':4, 'complex':5, 'LGN':6, None: 7}
typelabels = ['fast', 'slow', 'fast asym', 'slow asym',
              'simple', 'complex', 'LGN aff', 'unknown']
spiketypelabels = typelabels[:4]
rftypelabels = typelabels[4:]

def psthcorr(rec, nids=None, ssnids=None, ssseps=None, natexps=False, strange=None, plot=True):
    if nids == None:
        nids = sorted(rec.n) # use active neurons
    if ssnids == None:
        ssnids = nids # use nids as the superset
    nn = len(nids)
    nnss = len(ssnids)
    midbins, psths = rec.traster(nids=nids, natexps=natexps, strange=strange, plot=False,
                                 psth=True, binw=0.02, tres=0.005, norm=True)
    rho = np.corrcoef(psths) # defaults to bias=1
    rho[np.diag_indices(nn)] = np.nan # nan the diagonal, which imshow plots as white
    ssrho = np.zeros((nnss, nnss)) # superset rho matrix
    ssrho.fill(np.nan) # init with nans
    # load up values into appropriate spots in superset rho matrix:
    for i in range(nn):
        for j in range(nn):
            ssi, ssj = ssnids.searchsorted([nids[i], nids[j]])
            ssrho[ssi, ssj] = rho[i, j]

    if plot == False:
        return ssrho

    # plot superset rho matrix:
    figure(figsize=FIGSIZE)
    imshow(ssrho, vmin=-1, vmax=1, cmap='jet') # cmap='gray' is too bland
    ssnidticks = np.arange(0, nnss, 10)
    xticks(ssnidticks)
    yticks(ssnidticks)
    if SHOWCOLORBAR:
        colorbar()
    basetitle = rec.absname
    if strange != None:
        strange_sec = tuple(np.array(strange)/1e6) # convert to sec for display
        basetitle += '_strange=(%.f, %.f)' % strange_sec
    gcfm().window.setWindowTitle(basetitle + '_rho_mat')
    tight_layout(pad=0.3)

    # plot rho histogram:
    lti = np.tril_indices(nnss, -1) # lower triangle (below diagonal) indices off ssrho
    ssrhol = ssrho[lti]
    notnanis = np.logical_not(np.isnan(ssrhol)) # indices of non-nan values
    fssrhol = ssrhol[notnanis] # ssrhol filtered out for nans
    figure(figsize=FIGSIZE)
    rhobins = np.arange(RHOMIN, RHOMAX+0.0333, 0.0333) # left edges + rightmost edge
    n = hist(fssrhol, bins=rhobins, color='k')[0]
    axvline(x=fssrhol.mean(), c='r', ls='--') # draw vertical red line at mean fssrhol
    axvline(x=0, c='e', ls='--') # draw vertical grey line at x=0
    xlim(xmin=RHOMIN, xmax=RHOMAX)
    ylim(ymax=n.max()) # effectively normalizes the histogram
    rhoticks = np.arange(-0.2, 1+0.2, 0.2) # excluding the final 1
    xticks(rhoticks)
    yticks([n.max()]) # turn off y ticks to save space
    #yticks([0, n.max()])
    gcfm().window.setWindowTitle(basetitle + '_rho_hist')
    tight_layout(pad=0.3)

    # plot rho vs separation:
    fssseps = ssseps[notnanis] # ssseps filtered out for nans
    figure(figsize=FIGSIZE)
    # scatter plot:
    pl.plot(fssseps, fssrhol, 'k.')
    # bin seps and plot mean rho in each bin:
    sortis = np.argsort(fssseps)
    seps = fssseps[sortis]
    rhos = fssrhol[sortis]
    sepbins = np.arange(0, seps.max()+SEPBINW, SEPBINW) # left edges
    sepis = seps.searchsorted(sepbins)
    sepmeans, rhomeans, rhostds = [], [], []
    for sepi0, sepi1 in zip(sepis[:-1], sepis[1:]): # iterate over sepbins
        sepmeans.append(seps[sepi0:sepi1].mean()) # mean sep of all points in this sepbin
        rhoslice = rhos[sepi0:sepi1] # rhos in this sepbin
        rhomeans.append(rhoslice.mean()) # mean rho of all points in this sepbin
        rhostds.append(rhoslice.std()) # std of rho in this sepbin
    #pl.plot(sepmeans, rhomeans, 'r.-', ms=10, lw=2)
    errorbar(sepmeans, rhomeans, yerr=rhostds, fmt='r.-', ms=10, lw=2, zorder=9999)
    xlim(xmin=0, xmax=SEPXMAX)
    ylim(ymin=RHOMIN, ymax=RHOMAX)
    septicks = np.arange(0, seps.max()+100, 500)
    xticks(septicks)
    yticks(rhoticks)
    gcfm().window.setWindowTitle(basetitle + '_rho_sep')
    tight_layout(pad=0.3)
    return ssrho

def psthcorrdiff(rhos, seps, basetitle):
    """Plot difference of a pair of rho matrices (rhos[0] - rhos[1]). seps is the
    corresponding distance matrix"""
    assert len(rhos) == 2

    # calc rho difference matrix:
    rhod = rhos[0] - rhos[1]
    assert rhod.shape[0] == rhod.shape[1] # square
    nn = rhod.shape[0]
    
    # plot rho diff matrix:
    figure(figsize=FIGSIZE)
    imshow(rhod, vmin=-1, vmax=1, cmap='jet') # cmap='gray' is too bland
    ssnidticks = np.arange(0, nn, 10)
    xticks(ssnidticks)
    yticks(ssnidticks)
    if SHOWCOLORBAR:
        colorbar()
    gcfm().window.setWindowTitle(basetitle + '_rhod_mat')
    tight_layout(pad=0.3)

    # plot rho histogram:
    lti = np.tril_indices(nn, -1) # lower triangle (below diagonal) indices
    rhol = rhod[lti]
    notnanis = np.logical_not(np.isnan(rhol)) # indices of non-nan values
    frhol = rhol[notnanis] # rhol filtered out for nans
    figure(figsize=FIGSIZE)
    rhobins = np.arange(RHODIFFMIN, RHODIFFMAX+0.0333, 0.0333) # left edges + rightmost edge
    n = hist(frhol, bins=rhobins, color='k')[0]
    axvline(x=frhol.mean(), c='r', ls='--') # draw vertical red line at mean frhol
    axvline(x=0, c='e', ls='--') # draw vertical grey line at x=0
    xlim(xmin=RHODIFFMIN, xmax=RHODIFFMAX)
    ylim(ymax=n.max()) # effectively normalizes the histogram
    rhoticks = np.arange(-0.6, 0.6+0.2, 0.2)
    xticks(rhoticks)
    yticks([n.max()]) # turn off y ticks to save space
    #yticks([0, n.max()])
    gcfm().window.setWindowTitle(basetitle + '_rhod_hist')
    tight_layout(pad=0.3)
    
    # plot rho vs separation:
    fseps = seps[notnanis] # seps filtered out for nans
    figure(figsize=FIGSIZE)
    # scatter plot:
    pl.plot(fseps, frhol, 'k.')
    # bin seps and plot mean rho in each bin:
    sortis = np.argsort(fseps)
    seps = fseps[sortis]
    rhos = frhol[sortis]
    sepbins = np.arange(0, seps.max()+SEPBINW, SEPBINW) # left edges
    sepis = seps.searchsorted(sepbins)
    sepmeans, rhomeans, rhostds = [], [], []
    for sepi0, sepi1 in zip(sepis[:-1], sepis[1:]): # iterate over sepbins
        sepmeans.append(seps[sepi0:sepi1].mean()) # mean sep of all points in this sepbin
        rhoslice = rhos[sepi0:sepi1] # rhos in this sepbin
        rhomeans.append(rhoslice.mean()) # mean rho of all points in this sepbin
        rhostds.append(rhoslice.std()) # std of rho in this sepbin
    #pl.plot(sepmeans, rhomeans, 'r.-', ms=10, lw=2)
    errorbar(sepmeans, rhomeans, yerr=rhostds, fmt='r.-', ms=10, lw=2, zorder=9999)
    xlim(xmin=0, xmax=SEPXMAX)
    ylim(ymin=RHODIFFMIN, ymax=RHODIFFMAX)
    septicks = np.arange(0, seps.max()+100, 500)
    xticks(septicks)
    yticks(rhoticks)
    gcfm().window.setWindowTitle(basetitle + '_rhod_sep')
    tight_layout(pad=0.3)

def psthcorrtype(trackrecs, pool=False, alpha=0.0005, vmin=0, vmax=1, separatetypeplots=True):
    """Plot mean PSTH correlation (rho) 2D histograms, indexed by spike and RF type. Plot one
    for each set of recordings in trackrecs (ostensibly, one per track). If pool, plot only
    one rho celltype histogram pooled across all trackrecs. If pool, return rhotype matrix
    filled only with those distributions whose mean is significantly different from 0."""
    ntracks = len(trackrecs)
    tracknames = [ trackrec[0].tr.absname for trackrec in trackrecs ]
    rhotype = listarr(np.empty((8, 8))) # init rho cell type 2D array of lists
    npairs = 0 # init npairs
    for tracki, recs in enumerate(trackrecs):
        track = recs[0].tr
        natexps = False
        trackname = recs[0].tr.absname
        if trackname == 'ptc15.tr7c':
            natexps = True
        ssnids, recsecnids = get_ssnids(recs)
        ssrhos = []
        for rec in recs:
            ssrho = psthcorr(rec, nids=None, ssnids=ssnids, natexps=natexps, plot=False)
            ssrhos.append(ssrho)
        ssrhos = np.asarray(ssrhos) # convert to 3D array
        if pool == False:
            listarr(rhotype) # reset between tracks
            npairs = 0 # reset between tracks
        nn = len(ssnids)
        nanis = np.isnan(ssrhos) # indices of nan values
        ssrhos[nanis] = 0 # replace nans with 0s
        maxabsssrhos = core.maxabs(ssrhos, axis=0) # keep only the max rho of each cell pair
        alln = track.alln
        for i in range(nn):
            ni = alln[ssnids[i]] # neuron i
            si = celltype2int[ni.spiketype]
            ri = celltype2int[ni.rftype]
            for j in range(i+1, nn): # use only upper triangle, don't double count cell pairs
                nj = alln[ssnids[j]] # neuron j
                sj = celltype2int[nj.spiketype]
                rj = celltype2int[nj.rftype]
                rho = maxabsssrhos[i, j]
                if rho == 0:
                    # ignore this cell pair's rho (they were never simultaneously active) so it
                    # doesn't mess up the celltype stats
                    continue
                # fill in the upper triangle of rhotype matrix:
                rhotype[si, sj].append(rho)
                rhotype[ri, rj].append(rho)
                # these cross terms should best be left disabled, because they conflate the
                # correlations between spiketype and rftype:
                #rhotype[ri, sj].append(rho)
                #rhotype[si, rj].append(rho)
                npairs += 1
        rhotypemeans = np.zeros(rhotype.shape); rhotypemeans.fill(nan)
        rhotypestds = np.zeros(rhotype.shape); rhotypestds.fill(nan)
        rhotypeps = np.zeros(rhotype.shape); rhotypeps.fill(np.inf)
        rhotypesigmeans = np.zeros(rhotype.shape); rhotypesigmeans.fill(nan)
        # calculate rho stats for each combination of cell type:
        for i in range(8):
            for j in range(i, 8): # use only upper triangle, don't double count celltype stats
                if len(rhotype[i, j]) > 1:
                    rhotypemeans[i, j] = np.mean(rhotype[i, j])
                    rhotypestds[i, j] = np.std(rhotype[i, j])
                    # 2-sided sample mean ttest relative to 0:
                    t, p = ttest_1samp(rhotype[i, j], 0)
                    rhotypeps[i, j] = p
        sigis = rhotypeps < alpha # indices of significant deviations of mean from 0
        rhotypesigmeans[sigis] = rhotypemeans[sigis]
        #arrs = [rhotypemeans, rhotypestds, rhotypeps, rhotypesigmeans]
        #plottypes = ['mean', 'stdev', 'pval', 'sigmean']
        arrs = [rhotypesigmeans]
        plottypes = ['sigmean']
        if pool:
            if tracki < ntracks-1:
                continue # only plot once all tracks have been iterated over
            trackname = ', '.join(tracknames)
        for arr, plottype in zip(arrs, plottypes):
            # get symmetric arr by copying upper triangle, transposing to get lower triangle,
            # and adding to arr:
            symarr = nansum([arr, np.triu(arr, k=1).T], axis=0)
            thisvmin, thisvmax = nanmin(symarr), nanmax(symarr)
            vmin = min(vmin, thisvmin) # set to manual vmin at most
            vmax = max(vmax, thisvmax) # set to manual vmax at least
            if separatetypeplots:
                figure(figsize=(8, 3))
                # plot spiketypes:
                subplot(121)
                imshow(symarr[:4, :4], vmin=vmin, vmax=vmax, origin='upper', cmap='jet')
                xticks(np.arange(4), spiketypelabels, rotation=90)
                yticks(np.arange(4), spiketypelabels)
                colorbar(ticks=(vmin, vmax), format='%.2f')
                # plot rftypes:
                subplot(122)
                imshow(symarr[4:, 4:], vmin=vmin, vmax=vmax, origin='upper', cmap='jet')
                xticks(np.arange(4), rftypelabels, rotation=90)
                yticks(np.arange(4), rftypelabels)
                colorbar(ticks=(vmin, vmax), format='%.2f')
                plottype += ' separate'
            else: # plot spike and rf types in the same matrix
                figure(figsize=(4, 4))
                imshow(symarr, vmin=vmin, vmax=vmax, origin='upper', cmap='jet')
                xticks(np.arange(8), typelabels, rotation=90)
                yticks(np.arange(8), typelabels)
                colorbar(ticks=(vmin, vmax), format='%.2f')
                plottype += ' combined'
            titlestr = (trackname + ' psthcorrtype ' + plottype +
                        ' alpha=%.4f, npairs=%d' % (alpha, npairs))
            gcfm().window.setWindowTitle(titlestr)
            tight_layout(pad=0.4)
    if pool: # only return rhotype if pooling across all tracks
        insigis = np.logical_not(sigis)
        rhotype[insigis] = listarr(rhotype[insigis]) # set insig entries to empty lists
        return rhotype # only significant entires aren't empty

def psthcorrtypestats(rhotype, sigiss=None, test=ttest_ind, alpha=0.01):
    """Run some statistical tests on rhotype matrix. First, test if suspected significantly
    high entries of spike types and RF types (sigiss) in rhotype array do indeed have
    significantly different means from all the other entries pooled together. Also, test each
    entry against each other entry for significant differences in mean. Test can be ttest_ind
    or ks_2samp (or probably many others). Best is probably ttest_ind, because we're
    interested in testing whether differences in means are significant. KS tests if overall
    differences in distributions are significant. T-test seems to return more stringent p
    values, and therefore seems to be the more conservative choice in this case."""
    rhosptype = rhotype[:4, :4] # spike type, top left part of upper triangle
    rhorftype = rhotype[4:, 4:] # RF type, right bottom part of upper triangle
    rhosubtypes = [rhosptype, rhorftype]
    typelabelss = [spiketypelabels, rftypelabels]
    spsigis, rfsigis = sigiss # unpack indices of suspected significant entries
    spinsigis, rfinsigis = np.logical_not(spsigis), np.logical_not(rfsigis)
    insigiss = [spinsigis, rfinsigis]
    print('alpha=%s' % alpha)
    for typelabels, rhosubtype, sigis, insigis in zip(typelabelss, rhosubtypes,
                                                      sigiss, insigiss):
        # collect distributions of suspected insignificantly low entries in rhotype,
        # to then compare to the suspected significantly high entries.
        insig = []
        for i in range(4):
            for j in range(i, 4): # iterate over upper triangle of insigis
                if insigis[i, j]: # suspected insignificantly low entry
                    # some entires in rhosubtype may be empty due to T-test in psthcorrtype():
                    insig.append(rhosubtype[i, j])
        insig = np.hstack(insig)
        # test each suspected significantly high entry against insig distrib
        for i in range(4):
            for j in range(i, 4): # iterate over upper triangle of sigis
                if sigis[i, j]: # suspected significantly high entry
                    sig = rhosubtype[i, j]
                    t, p = test(sig, insig)
                    ti, tj = typelabels[i], typelabels[j]
                    if p < alpha:
                        result = 'significantly different'
                    else:
                        result = 'NOT significantly different'
                    print('%s-%s %s from pooled others: p=%s' % (ti, tj, result, p))
    # for more statistical power, instead of pooling over all insigis, compare each
    # entry to every other entry:
    for typelabels, rhosubtype in zip(typelabelss, rhosubtypes):
        print('--')
        for ai in range(4):
            for aj in range(ai, 4): # iterate over upper triangle of rhosubtype
                a = rhosubtype[ai, aj]
                if len(a) == 0: continue # skip empty entries
                tai, taj = typelabels[ai], typelabels[aj]
                for bi in range(4):
                    for bj in range(bi, 4): # iterate over upper triangle of rhosubtype
                        b = rhosubtype[bi, bj]
                        if len(b) == 0: continue # skip empty entries
                        tbi, tbj = typelabels[bi], typelabels[bj]
                        t, p = test(a, b)
                        if p < alpha:
                            result = 'significantly different'
                        else:
                            result = 'NOT significantly different'
                        print('%s-%s %s from %s-%s: p=%s, N=%r'
                              % (tai, taj, result, tbi, tbj, p, (len(a), len(b))))

def get_seps(ssnids, nd):
    """Build flattened array of distances between all unique pairs in ssnids, given neuron
    dict nd"""
    nnss = len(ssnids)
    lti = np.tril_indices(nnss, -1) # lower triangle (below diagonal) indices, ie unique pairs
    seps = []
    for nidii0, nidii1 in np.asarray(lti).T:
        sep = dist(nd[ssnids[nidii0]].pos, nd[ssnids[nidii1]].pos)
        seps.append(sep)
    seps = np.hstack(seps)
    return seps

listarr = np.frompyfunc(lambda x: [], 1, 1) # take 1 input array, return 1 list in each entry


# ptc15.tr7c:
SEPXMAX = 1675
# get superset of active nids for all natexps of both recs in ptc15tr7crecs:
stranges = etrangesr74 + etrangesr95b # 8 stranges in total
recs = [ptc15.tr7c.r74]*4 + [ptc15.tr7c.r95b]*4 # 8 recs corresponding to 8 stranges
ssnids, recsecnids = get_ssnids(recs, stranges)
ssseps = get_seps(ssnids, ptc15.tr7c.alln)
# get separate supersets of active nids for all 4 natexps in each recording:
ptc15tr7crecsecnids = [np.unique(np.hstack(recsecnids[:4])),
                       np.unique(np.hstack(recsecnids[4:]))]
# do psthcorr plots and collect ssrho matrices:
ssrhos = []
for rec, nids in zip(ptc15tr7crecs, ptc15tr7crecsecnids):
    ssrho = psthcorr(rec, nids=nids, ssnids=ssnids, ssseps=ssseps, natexps=True) # in sec
    ssrhos.append(ssrho)
# plot differences in superset rho matrices for the two recordings:
psthcorrdiff(ssrhos, ssseps, 'r74-r95b')

## rho for ns1 figure:
#In [124]: np.where(ssnids == 328)
#Out[124]: (array([31]),)
#In [125]: np.where(ssnids == 345)
#Out[125]: (array([34]),)
#In [126]: ssrhos[0][31,34]
#Out[126]: -0.11554163740884685

## rho for ns2 figure:
#In [127]: np.where(ssnids == 87)
#Out[127]: (array([8]),)
#In [128]: np.where(ssnids == 93)
#Out[128]: (array([10]),)
#In [129]: ssrhos[0][8,10]
#Out[129]: 0.82907446734056678

'''
# ptc22.tr1.r08 sections:
SEPXMAX = 1200
ssnids, recsecnids = get_ssnids(ptc22tr1r08s, strangesr08s)
ssseps = get_seps(ssnids, ptc22.tr1.alln)
for rec, nids, strange in zip(ptc22tr1r08s, recsecnids, strangesr08s):
    psthcorr(rec, nids=nids, ssnids=ssnids, ssseps=ssseps, natexps=False, strange=strange)

# ptc22.tr1.r10 sections:
SEPXMAX = 1200
ssnids, recsecnids = get_ssnids(ptc22tr1r10s, strangesr10s)
ssseps = get_seps(ssnids, ptc22.tr1.alln)
for rec, nids, strange in zip(ptc22tr1r10s, recsecnids, strangesr10s):
    psthcorr(rec, nids=nids, ssnids=ssnids, ssseps=ssseps, natexps=False, strange=strange)
'''
# ptc22.tr1.r08 + ptc22.tr1.r10 sections:
plotpsthcorr = True
plotpsthcorrdiff = True
SEPXMAX = 1200
ptc22tr1s = ptc22tr1r08s+ptc22tr1r10s
stranges = strangesr08s+strangesr10s
ssnids, recsecnids = get_ssnids(ptc22tr1s, stranges)
ssseps = get_seps(ssnids, ptc22.tr1.alln)
# do psthcorr plots and collect ssrho matrices:
ssrhos = []
for rec, nids, strange in zip(ptc22tr1s, recsecnids, stranges):
    ssrho = psthcorr(rec, nids=nids, ssnids=ssnids, ssseps=ssseps, natexps=False,
                     strange=strange, plot=plotpsthcorr)
    ssrhos.append(ssrho)
ssrhos = np.asarray(ssrhos) # convert to 3D array
if plotpsthcorrdiff:
    # plot differences in superset rho matrices for various pairs of recording sections:
    psthcorrdiff([ssrhos[1], ssrhos[0]], ssseps, 'B-A')
    psthcorrdiff([ssrhos[2], ssrhos[1]], ssseps, 'C-B')
    psthcorrdiff([ssrhos[3], ssrhos[2]], ssseps, 'D-C')
    #psthcorrdiff([ssrhos[3], ssrhos[0]], ssseps, 'D-A')
    #psthcorrdiff([ssrhos[3], ssrhos[1]], ssseps, 'D-B')
    #psthcorrdiff([ssrhos[2], ssrhos[0]], ssseps, 'C-A')


# run psthcorrtype and psthcorrtypestats on ptc15.tr7c:
trackrecs = [ptc15tr7crecs]
rhotype = psthcorrtype(trackrecs, pool=True, alpha=0.0005, vmin=0, vmax=0.13,
                       separatetypeplots=True)
spsigis, rfsigis = np.zeros((4,4), dtype=bool), np.zeros((4,4), dtype=bool)
print('\nptc15.tr7c')
psthcorrtypestats(rhotype, sigiss=[spsigis, rfsigis], test=ttest_ind, alpha=0.01)

# run psthcorrtype and psthcorrtypestats on ptc22.tr1:
trackrecs = [ptc22tr1recs]
rhotype = psthcorrtype(trackrecs, pool=True, alpha=0.0005, vmin=0, vmax=0.13,
                       separatetypeplots=True)
spsigis, rfsigis = np.zeros((4,4), dtype=bool), np.zeros((4,4), dtype=bool)
spsigis[0, 0] = True
rfsigis[1, 1] = True; #rfsigis[0, 2] = True
print('\nptc22.tr1')
psthcorrtypestats(rhotype, sigiss=[spsigis, rfsigis], test=ttest_ind, alpha=0.01)

# run psthcorrtype and psthcorrtypestats on ptc22.tr2:
trackrecs = [ptc22tr2recs]
rhotype = psthcorrtype(trackrecs, pool=True, alpha=0.0005, vmin=0, vmax=0.13,
                       separatetypeplots=True)
spsigis, rfsigis = np.zeros((4,4), dtype=bool), np.zeros((4,4), dtype=bool)
spsigis[0, 0] = True
rfsigis[1, 1] = True; #rfsigis[0, 2] = True
print('\nptc22.tr2')
psthcorrtypestats(rhotype, sigiss=[spsigis, rfsigis], test=ttest_ind, alpha=0.01)

# run psthcorrtype and psthcorrtypestats on ptc22:
trackrecs = [ptc22tr1recs, ptc22tr2recs]
rhotype = psthcorrtype(trackrecs, pool=True, alpha=0.0005, vmin=0, vmax=0.13,
                       separatetypeplots=True)
spsigis, rfsigis = np.zeros((4,4), dtype=bool), np.zeros((4,4), dtype=bool)
spsigis[0, 0] = True
rfsigis[1, 1] = True; #rfsigis[0, 2] = True
print('\nptc22')
psthcorrtypestats(rhotype, sigiss=[spsigis, rfsigis], test=ttest_ind, alpha=0.01)

# run psthcorrtype and psthcorrtypestats on all tracks:
trackrecs = [ptc15tr7crecs, ptc22tr1recs, ptc22tr2recs]
rhotype = psthcorrtype(trackrecs, pool=True, alpha=0.0005, vmin=0, vmax=0.13,
                       separatetypeplots=True)
spsigis, rfsigis = np.zeros((4,4), dtype=bool), np.zeros((4,4), dtype=bool)
print('\nall tracks pooled')
psthcorrtypestats(rhotype, sigiss=[spsigis, rfsigis], test=ttest_ind, alpha=0.01)


show()

