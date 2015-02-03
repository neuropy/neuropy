"""Examine trial-averaged LFP during synched and desynched state of the 2 natural scene movies
in ptc22.tr1. Run from within neuropy using `run -i scripts/lfp_precision.py`"""

from __future__ import division, print_function

from core import sparseness

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
    LFPf = pl.figure(figsize=figsize)
    LFPa = LFPf.add_subplot(111)
    SNf = pl.figure(figsize=figsize)
    SNa = SNf.add_subplot(111)
    for rec, strange, fmt, delta in zip(recs, stranges, fmts, (0.5, -0.5)):
        t, lfptrials = rec.tlfp(trange=strange, plot=False)
        lfpmean, lfpstd = lfptrials.mean(axis=0), lfptrials.std(axis=0)
        ntrials = len(lfptrials)
        # to make saturation represent LFP reliability, scale transparency inversely
        # with ntrials:
        alpha = 10 / ntrials
        # plot all trials for this rec, in mV
        LFPa.plot(t, lfptrials.T/1e3+delta, fmt, alpha=alpha)
        LFPa.plot(t, lfpmean/1e3+delta, 'w-', alpha=1)
        LFPa.plot(t, (lfpmean+lfpstd)/1e3+delta, 'k-', alpha=1)
        LFPa.plot(t, (lfpmean-lfpstd)/1e3+delta, 'k-', alpha=1)
        LFPa.set_xlim(xmax=5.5)
        LFPa.set_ylim(-1, 1) # mV
        LFPa.set_xlabel("time (sec)")
        LFPa.set_ylabel("LFP (mV)")
        LFPf.canvas.manager.set_window_title("tLFP %s" % rec.absname)
        LFPf.tight_layout(pad=0.3) # crop figure to contents
        # Fano-factor and CV don't work very well when mean approaches zero
        #SN = np.abs(lfpmean) / lfpstd # something like S/N ratio
        SN = lfpmean**2 / lfpstd**2 # something like S/N ratio
        SNa.plot(t, SN, fmt)
        SNa.set_xlim(xmax=5.5)
        SNa.set_xlabel("time (sec)")
        #SNa.set_ylabel("trial-averaged LFP $\mu^{2}/\sigma^{2}$")
        SNa.set_ylabel("trial-averaged LFP S/N")
        SNf.canvas.manager.set_window_title("tLFP SN %s" % rec.absname)
        print('sparseness:', sparseness(np.abs(lfpmean)))
        print('mean tLFP S/N:', SN.mean())
        SNf.tight_layout(pad=0.3) # crop figure to contents

pl.show()
