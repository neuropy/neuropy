"""Examine trial-averaged LFP during synched and desynched state of the 2 natural scene movies
in ptc22.tr1. Run from within neuropy using `run -i scripts/lfp_precision.py`"""

from __future__ import division, print_function

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

figsize=(10, 3)
for recs, stranges, fmts in zip((ptc22tr1r08s, ptc22tr1r10s),
                                (strangesr08s, strangesr10s),
                                (('b-', 'r-'), ('r-', 'b-'))):
    f = pl.figure(figsize=figsize)
    a = f.add_subplot(111)
    for rec, strange, fmt in zip(recs, stranges, fmts):
        t, lfpmean, lfpstd = rec.tlfp(trange=strange, plot=False)
        #plot_tlfp(t, lfpmean, f=f)
        #plot_tlfp(t, lfpmean+lfpstd, fmt='r-', f=f)
        #plot_tlfp(t, lfpmean-lfpstd, fmt='r-', f=f)
        # Fano-factor and CV don't work very well when mean approaches zero
        SN = lfpmean**2 / lfpstd**2 # something like S/N ratio
        a.plot(t, SN, fmt)
        a.set_xlim(xmax=5.5)
        a.set_xlabel("time (sec)")
        a.set_ylabel("trial-averaged LFP $\mu^{2}/\sigma^{2}$")
        gcfm().window.setWindowTitle("tLFP SN %s" % rec.absname)
        print('mean tLFP S/N:', SN.mean())
        f.tight_layout(pad=0.3) # crop figure to contents

pl.show()
