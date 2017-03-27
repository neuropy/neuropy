"""Plot PSTHs of specific neurons in specific natscene recordings, separately for synched and
desynched states. Run from within neuropy using `run -i scripts/psth_natstate_plots.py `"""

from psth_funcs import get_psth_peaks_gac

BINW, TRES = 0.020, 0.0001 # s
MINTHRESH = 3 # peak detection thresh, Hz
MEDIANX = 2 # PSTH median multiplier, Hz
figsize = 3, 1.5
ms = 8
sigma = 0.020

ns = [ ptc22.tr1.r08.alln[17], ptc22.tr1.r08.alln[20], ptc22.tr1.r08.alln[94], # natstate1
       ptc22.tr1.r10.alln[17], ptc22.tr1.r10.alln[20], ptc22.tr1.r10.alln[94],
       ptc22.tr1.r08.alln[3], ptc22.tr1.r08.alln[128], ptc22.tr1.r08.alln[175], # natstate2
       ptc22.tr1.r10.alln[3], ptc22.tr1.r10.alln[128], ptc22.tr1.r10.alln[175],
       ptc18.tr1.r38.alln[103], ptc18.tr1.r38.alln[178], ptc18.tr1.r38.alln[203] ] # natstate3

for n in ns:
    nid = n.id
    rec = n.sort.r
    stranges = REC2STATETRANGES[rec.absname] # [desynched, synched]
    figure(figsize=figsize)
    for strange, c in zip(stranges[::-1], ['r', 'b']): # [synched, desynched]
        # calculate PSTH:
        t, psths, spikets = rec.psth(nids=[nid], strange=strange, binw=BINW, tres=TRES,
                                     gauss=True, norm='ntrials', plot=False)
        psth, ts = psths[0], spikets[0] # lists of len 1
        # run PSTH peak detection:
        baseline = MEDIANX * np.median(psth)
        thresh = baseline + MINTHRESH # peak detection threshold
        print("n%d" % nid, end='')
        peakis, lis, ris = get_psth_peaks_gac(ts, t, psth, thresh, sigma=sigma)
        plot(t, psth, c=c, marker=None, ls='-')
        if len(peakis) > 0:
            plot(t[peakis], psth[peakis], c=c, marker='.', linestyle='none', ms=ms, mec='none')
        '''
        if len(lis) > 0:
            plot(t[lis], psth[lis], c='e', marker='.', linestyle='none', ms=ms, mec='none')
        if len(ris) > 0:
            plot(t[ris], psth[ris], c='k', marker='.', linestyle='none', ms=ms, mec='none')
        '''
    xlim(xmin=0, xmax=t[-1])
    #ylim(0, 1) # arbitrary units
    #yticks([0, 1])
    yticks([])
    gcfm().window.setWindowTitle(rec.absname+'_n'+str(nid)+'_PSTHs')
    tight_layout(pad=0.3) # crop figure to contents

    show()
