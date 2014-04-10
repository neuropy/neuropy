"""Plot 2D matrix of natural scene movie PSTH correlations of all active cells. Run from
within neuropy using `run -i scripts/psthcorr.py`"""

from __future__ import division

figsize = (3, 3)
common = False # use common set of active nids within a recording list?
showcolorbar = False # show colorbar
sepbinw = 200 # separation bin width, um

#ptc15tr7crecs = [ptc15.tr7c.r74, ptc15.tr7c.r95b]

ptc22tr1recs = [ptc22.tr1.r08, ptc22.tr1.r08, ptc22.tr1.r10, ptc22.tr1.r10]
stranges = [(0, 1500), # r08 desynched
            (1550, np.inf), # r08 synched, end is ~ 2300
            (0, 1400), # r10 synched
            (1480, np.inf)] # r10 desynched, end is ~ 2300

#ptc22tr2recs  = [ptc22.tr2.r33, ptc22.tr2.r28] # 28 is a 5 min movie
"""
## TODO: make rec method that takes strange and returns list of nids that are active within
that trange. The, do this for every rec section to be plotted here. For the rho matrix, use
the superset of all nids to be compared, and leave set to -1 for those rec sections for which
a nid is inactive. This will result in entire rows and columns coloured dark blue, but will
allow direct comparison of as many cell pairs as possible across rec sections. Would be nice
to color them black or white instead of dark blue. Hard to do with the cmap? Maybe used masked
arrays. What would a cmap do with an entry that's masked? Leave it white like the background?
"""

def psthcorr(rec, nids=None, natexps=False, strange=None):
    if nids == None:
        nids = sorted(rec.n)
    midbins, psths = rec.traster(nids=nids, natexps=natexps, strange=strange, plot=False,
                                 psth=True, binw=0.02, tres=0.005, norm=True)
    nn = len(nids)
    nidticks = np.arange(0, nn, 10)
    rho = np.corrcoef(psths) # defaults to bias=1
    rho[np.diag_indices(nn)] = 0.0 # null the diagonal

    # plot rho matrix:
    figure(figsize=figsize)
    imshow(rho, vmin=-1, vmax=1, cmap='jet') # cmap='gray' is too bland
    xticks(nidticks)
    yticks(nidticks)
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
    rhobins = np.arange(-0.3, 1+0.0333, 0.0333) # left edges + rightmost edge
    n = hist(rhol, bins=rhobins, color='k')[0]
    axvline(x=rhol.mean(), c='r', ls='--') # draw vertical line at mean rho
    xlim(xmin=-0.3, xmax=1)
    ylim(ymax=n.max()) # effectively normalizes the histogram
    rhoticks = np.arange(-0.2, 1, 0.2) # excluding the final 1
    xticks(rhoticks)
    yticks([]) # turn off y ticks for space
    #yticks([0, n.max()])
    gcfm().window.setWindowTitle(basetitle + '_rho_hist')
    tight_layout(pad=0.3)

    # plot rho vs separation:
    seps = []
    for nidii0, nidii1 in np.asarray(lti).T:
        sep = dist(rec.n[nids[nidii0]].pos, rec.n[nids[nidii1]].pos)
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
    sepmeans, rhomeans, rhosems = [], [], []
    for sepi0, sepi1 in zip(sepis[:-1], sepis[1:]): # iterate over sepbins
        sepmeans.append(seps[sepi0:sepi1].mean()) # mean sep of all points in this sepbin
        rhoslice = rhos[sepi0:sepi1] # rhos in this sepbin
        rhomeans.append(rhoslice.mean()) # mean rho of all points in this sepbin
        rhosems.append(rhoslice.std() / np.sqrt(len(rhoslice))) # SEM of rho in this sepbin
    #plot(sepmeans, rhomeans, 'r.-', ms=10, lw=2)
    errorbar(sepmeans, rhomeans, yerr=rhosems, fmt='r.-', ms=10, lw=2)
    ylim(ymin=-0.3, ymax=1)
    septicks = np.arange(0, seps.max()+100, 500)
    xticks(septicks)
    yticks(rhoticks)
    gcfm().window.setWindowTitle(basetitle + '_rho_sep')
    tight_layout(pad=0.3)

'''
nids = None
if common:
    nids = ptc15.tr7c.get_nids(['74', '95b']) # sorted
for rec in ptc15tr7crecs:
    psthcorr(rec, nids=nids, natexps=True)
'''
nids = None
if common:
    nids = ptc22.tr1.get_nids(['08', '10']) # sorted
for rec, strange in zip(ptc22tr1recs, stranges):
    psthcorr(rec, nids=nids, natexps=False, strange=strange)


show()
