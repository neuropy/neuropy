"""Plot 2D matrix of natural scene movie PSTH correlations of all active cells. Run from
within neuropy using `run -i scripts/psthcorr.py`"""

from __future__ import division

figsize = (3, 3)
showcolorbar = False # show colorbar
sepbinw = 200 # separation bin width, um
rhomin, rhomax = -0.3, 1

ptc15tr7crecs = [ptc15.tr7c.r74, ptc15.tr7c.r95b]
nateids = [3, 4, 10, 12] # for both recs in ptc15.tr7c
etrangesr74 = [ ptc15.tr7c.r74.e[nateid].trange for nateid in nateids ] # us
etrangesr95b = [ ptc15.tr7c.r74.e[nateid].trange for nateid in nateids ] # us

ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500), # r08 desynched
                (1550, np.inf)] # r08 synched, end is ~ 2300
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400), # r10 synched
                (1480, np.inf)] # r10 desynched, end is ~ 2300

#ptc22tr2recs  = [ptc22.tr2.r33, ptc22.tr2.r28] # 28 is a 5 min movie

def psthcorr(rec, nids=None, ssnids=None, natexps=False, strange=None):
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
    # load up relevant values into superset rho matrix:
    for i in range(nn):
        for j in range(nn):
            ssi, ssj = ssnids.searchsorted([nids[i], nids[j]])
            ssrho[ssi, ssj] = rho[i, j]

    # plot superset rho matrix:
    figure(figsize=figsize)
    imshow(ssrho, vmin=-1, vmax=1, cmap='jet') # cmap='gray' is too bland
    ssnidticks = np.arange(0, nnss, 10)
    xticks(ssnidticks)
    yticks(ssnidticks)
    if showcolorbar:
        colorbar()
    basetitle = rec.absname
    if strange != None:
        basetitle += '_strange=%s' % (strange,)
    gcfm().window.setWindowTitle(basetitle + '_rho_mat')
    tight_layout(pad=0.3)

    # plot rho histogram:
    lti = np.tril_indices(nn, -1) # lower triangle (below diagonal) indices
    rhol = rho[lti]
    figure(figsize=figsize)
    rhobins = np.arange(rhomin, rhomax+0.0333, 0.0333) # left edges + rightmost edge
    n = hist(rhol, bins=rhobins, color='k')[0]
    axvline(x=rhol.mean(), c='r', ls='--') # draw vertical red line at mean rho
    axvline(x=0, c='e', ls='--') # draw vertical grey line at x=0
    xlim(xmin=rhomin, xmax=rhomax)
    ylim(ymax=n.max()) # effectively normalizes the histogram
    rhoticks = np.arange(-0.2, 1+0.2, 0.2) # excluding the final 1
    xticks(rhoticks)
    yticks([n.max()]) # turn off y ticks to save space
    #yticks([0, n.max()])
    gcfm().window.setWindowTitle(basetitle + '_rho_hist')
    tight_layout(pad=0.3)

    # plot rho vs separation:
    seps = []
    for nidii0, nidii1 in np.asarray(lti).T:
        sep = dist(rec.alln[nids[nidii0]].pos, rec.alln[nids[nidii1]].pos)
        seps.append(sep)
    seps = np.hstack(seps)
    figure(figsize=figsize)
    # scatter plot:
    plot(seps, rhol, 'k.')
    # bin seps and plot mean rho in each bin:
    sortis = np.argsort(seps)
    seps = seps[sortis]
    rhos = rhol[sortis]
    sepbins = np.arange(0, seps.max()+sepbinw, sepbinw) # left edges
    sepis = seps.searchsorted(sepbins)
    sepmeans, rhomeans, rhostds = [], [], []
    for sepi0, sepi1 in zip(sepis[:-1], sepis[1:]): # iterate over sepbins
        sepmeans.append(seps[sepi0:sepi1].mean()) # mean sep of all points in this sepbin
        rhoslice = rhos[sepi0:sepi1] # rhos in this sepbin
        rhomeans.append(rhoslice.mean()) # mean rho of all points in this sepbin
        rhostds.append(rhoslice.std()) # std of rho in this sepbin
    #plot(sepmeans, rhomeans, 'r.-', ms=10, lw=2)
    errorbar(sepmeans, rhomeans, yerr=rhostds, fmt='r.-', ms=10, lw=2, zorder=9999)
    xlim(xmin=0, xmax=sepxmax)
    ylim(ymin=rhomin, ymax=rhomax)
    septicks = np.arange(0, seps.max()+100, 500)
    xticks(septicks)
    yticks(rhoticks)
    gcfm().window.setWindowTitle(basetitle + '_rho_sep')
    tight_layout(pad=0.3)

# ptc15.tr7c:
sepxmax = 1675
recsecnids = []
# get superset of active nids for all natexps of both recs in ptc15tr7crecs:
strangesus = etrangesr74 + etrangesr95b # 8 stranges in total
#stranges = np.asarray(strangesus) / 1e6 # in sec
recs = [ptc15.tr7c.r74]*4 + [ptc15.tr7c.r95b]*4 # 8 recs corresponding to 8 stranges
for rec, strangeus in zip(recs, strangesus):
    recsecnids.append(rec.get_nids(tranges=[strangeus])) # 8 recsecnids in total, in us
ssnids = np.unique(np.hstack(recsecnids)) # superset of active nids from all rec sections
# separate supersets of active nids for all 4 natexps in each recording:
ptc15tr7crecsecnids = [np.unique(np.hstack(recsecnids[:4])),
                       np.unique(np.hstack(recsecnids[4:]))]
for rec, nids in zip(ptc15tr7crecs, ptc15tr7crecsecnids):
    psthcorr(rec, nids=nids, ssnids=ssnids, natexps=True) # in sec

'''
# ptc22.tr1.r08 sections:
sepxmax = 1200
recsecnids = [] # holds arrays of active nids of each recording section
for rec, strange in zip(ptc22tr1r08s, strangesr08s):
    recsecnids.append(rec.get_nids(tranges=[np.asarray(strange) * 1000000])) # convert to us
ssnids = np.unique(np.hstack(recsecnids)) # superset of active nids from rec sections
for rec, nids, strange in zip(ptc22tr1r08s, recsecnids, strangesr08s):
    psthcorr(rec, nids=nids, ssnids=ssnids, natexps=False, strange=strange)

# ptc22.tr1.r10 sections:
sepxmax = 1200
recsecnids = [] # holds arrays of active nids of each recording section
for rec, strange in zip(ptc22tr1r10s, strangesr10s):
    recsecnids.append(rec.get_nids(tranges=[np.asarray(strange) * 1000000])) # convert to us
ssnids = np.unique(np.hstack(recsecnids)) # superset of active nids from rec sections
for rec, nids, strange in zip(ptc22tr1r10s, recsecnids, strangesr10s):
    psthcorr(rec, nids=nids, ssnids=ssnids, natexps=False, strange=strange)
'''
'''
# ptc22.tr1.r08 + ptc22.tr1.r10 sections:
sepxmax = 1200
recsecnids = [] # holds arrays of active nids of each recording section
ptc22tr1s = ptc22tr1r08s+ptc22tr1r10s
stranges = strangesr08s+strangesr10s
for rec, strange in zip(ptc22tr1s, stranges):
    recsecnids.append(rec.get_nids(tranges=[np.asarray(strange) * 1000000])) # convert to us
ssnids = np.unique(np.hstack(recsecnids)) # superset of active nids from rec sections
for rec, nids, strange in zip(ptc22tr1s, recsecnids, stranges):
    psthcorr(rec, nids=nids, ssnids=ssnids, natexps=False, strange=strange)
'''
"""
## TODO: take difference between pairs of PSTH corr matrices. To include only those cell pairs
that don't have nan as a value, for each difference matrix, take a[npnot(isnan(a))] to get
just such pairs. Then take the mean of this nan filtered difference matrix.
"""

show()
