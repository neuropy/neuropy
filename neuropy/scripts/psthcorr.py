"""Plot 2D matrix of natural scene movie PSTH correlations of all active cells. Run from
within neuropy using `run -i scripts/psthcorr.py`"""

from __future__ import division

figsize = (3, 3)
common = False # use common set of active nids within a recording list?
showcolorbar = False # show colorbar
sepbinw = 200 # separation bin width, um

ptc15tr7crecs = [ptc15.tr7c.r74, ptc15.tr7c.r95b]
ptc22tr2recs  = [ptc22.tr2.r33, ptc22.tr2.r28] # 28 is a 5 min movie
#synched = [ptc22.tr1.r08s, ptc22.tr1.r10s]
#desynched = [ptc22.tr1.r08d, ptc22.tr1.r10d]

if common:
    nids = ptc15.tr7c.get_nids(['74', '95b']) # sorted
    nn = len(nids)
    npairs = nCr(nn, 2)

for rec in ptc15tr7crecs:
    if not common:
        nids = sorted(rec.n)
        nn = len(nids)
        npairs = nCr(nn, 2)
    nidticks = np.arange(0, nn, 10)
    midbins, psths = rec.traster(nids=nids, natexps=True, plot=False, psth=True, binw=0.02,
                                 tres=0.005, norm=True)
    rho = np.corrcoef(psths) # defaults to bias=1

    # plot rho matrix:
    #rho[np.diag_indices(nn)] = 0.0 # null the diagonal
    figure(figsize=figsize)
    imshow(rho, vmin=-1, vmax=1, cmap='jet') # cmap='gray' is too bland
    xticks(nidticks)
    yticks(nidticks)
    if showcolorbar:
        colorbar()
    gcfm().window.setWindowTitle(rec.absname + '_rho_mat')
    tight_layout(pad=0.3)

    # plot rho histogram:
    lti = np.tril_indices(nn, -1) # lower triangle (below diagonal) indices
    rhol = rho[lti]
    figure(figsize=figsize)
    hist(rhol, bins=30, color='k')
    gcfm().window.setWindowTitle(rec.absname + '_rho_hist')
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
    sepmeans = []
    rhomeans = []
    for sepi0, sepi1 in zip(sepis[:-1], sepis[1:]): # iterate over sepbins
        sepmeans.append(seps[sepi0:sepi1].mean()) # average sep of all points in this bin
        rhomeans.append(rhos[sepi0:sepi1].mean()) # average rho of all points in this bin
    plot(sepmeans, rhomeans, 'r.-')
    septicks = np.arange(0, seps.max()+500, 500)
    xticks(septicks)
    gcfm().window.setWindowTitle(rec.absname + '_rho_sep')
    tight_layout(pad=0.3)

show()
