"""Examine trial-averaged LFP during synched and desynched state of the specified natural
scene movies. Run from within neuropy using `run -i scripts/lfp_precision.py`"""

from __future__ import division, print_function

from core import sparseness

# copied from psth_precision.py:
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

FIGSIZE = (10, 3)

reccmp = lambda reca, recb: cmp(reca.absname, recb.absname)
urecs = sorted(rec2tranges, cmp=reccmp) # unique recordings, no repetition, sorted
#urecnames = ' '.join([rec.absname for rec in urecs])
fmts = ('b-', 'r-') # desynched and synched

for rec in urecs:
    print(rec.absname)
    LFPf = pl.figure(figsize=FIGSIZE)
    LFPa = LFPf.add_subplot(111)
    SNRf = pl.figure(figsize=FIGSIZE)
    SNRa = SNRf.add_subplot(111)
    stranges = rec2tranges[rec]
    for strange, fmt, delta in zip(stranges, fmts, (0.5, -0.5)): # desynched, then synched
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
        SNR = np.abs(lfpmean) / lfpstd # something like S/N ratio
        #SNR = lfpmean**2 / lfpstd**2 # something like S/N ratio
        SNRa.plot(t, SNR, fmt)
        SNRa.set_xlim(xmax=5.5)
        SNRa.set_xlabel("time (sec)")
        SNRa.set_ylabel(r"LFP S/N ($|\mu|/\sigma)$")
        SNRf.canvas.manager.set_window_title("tLFP SNR %s" % rec.absname)
        print('sparseness:', sparseness(np.abs(lfpmean)))
        print('mean tLFP S/N:', SNR.mean())
        SNRf.tight_layout(pad=0.3) # crop figure to contents

pl.show()
