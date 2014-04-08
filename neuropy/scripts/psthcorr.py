"""Plot 2D matrix of natural scene movie PSTH correlations of all active cells. Run from
within neuropy using `run -i scripts/psthcorr.py`"""

from __future__ import division

figsize = (3, 3)
common = False # use common set of active nids within a recording list?
showcolorbar = False # show colorbar

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
    # pull out lower triangle (below diagonal), upper would work just as well:
    rhol = rho[np.tril_indices(nn, -1)]
    figure(figsize=figsize)
    hist(rhol.flatten(), bins=30, color='k')
    gcfm().window.setWindowTitle(rec.absname + '_rho_hist')
    tight_layout(pad=0.3)

show()
