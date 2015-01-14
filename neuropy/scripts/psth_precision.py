"""Detect response events within PSTHs, count them, measure the FWHM of each one, and plot
their distributions as a function of cortical state within each of the 2 natural scene movies
in ptc22.tr1"""

from __future__ import division, print_function
from scipy.signal import argrelextrema
from core import argfwhm, get_ssnids

# copied from psthcorr.py:
ptc22tr1r08s = [ptc22.tr1.r08, ptc22.tr1.r08]
strangesr08s = [(0, 1500e6), # r08 desynched, us
                (1550e6, np.inf)] # r08 synched, us, end is ~ 2300s
ptc22tr1r10s = [ptc22.tr1.r10, ptc22.tr1.r10]
strangesr10s = [(0, 1400e6), # r10 synched, us
                (1480e6, np.inf)] # r10 desynched, us, end is ~ 2300s

recs = ptc22tr1r08s + ptc22tr1r10s
stranges = strangesr08s + strangesr10s

NIDSKIND = 'all'

BINW, TRES = 0.02, 0.0002 # PSTH time bins, sec
MINTHRESH = 3 # peak detection thresh, Hz
#BASELINEX = 5 # PSTH baseline multiplier, Hz
FWFRACTION = 0.5 # full width fraction of max

# plotting params:
FWHMMIN, FWHMMAX, FWHMSTEP, XTICKSTEP = 0, 100, 2.5, 25
PLOTPSTH = False
SPARSTEP = 0.1
figsize = (3.5, 3.5) # inches


def sparseness(psth):
    """Sparseness measure, from Vinje and Gallant, 2000. This is basically 1 minus the ratio
    of the square of the sums over the sum of the squares of the values in the PSTH"""
    n = len(psth) # number of bins
    return (1 - (psth.sum()/n)**2 / np.sum((psth**2)/n)) / (1 - 1/n)

def plot_psth(psthparams, nid, fmt='k-', ms=10):
    psth, thresh, baseline, peakis, lis, ris, fwhms = psthparams[nid]
    figure(figsize=(23, 7))
    pl.plot(np.arange(0, len(psth))*TRES, psth, fmt)
    # plot thresh and baseline levels:
    axhline(y=thresh, c='r', ls='--')
    axhline(y=baseline, c='e', ls='--')
    # mark peaks and their edges:
    pl.plot(lis*TRES, psth[lis], 'co', ms=ms, mec='none') # left edges
    pl.plot(ris*TRES, psth[ris], 'bo', ms=ms, mec='none') # rigth edges
    pl.plot(peakis*TRES, psth[peakis], 'ro', ms=ms, mec='none') # peaks
    gcfm().window.setWindowTitle('n%d, thresh=%g, baseline=%g' % (nid, thresh, baseline))
    gcf().tight_layout(pad=0.3) # crop figure to contents

def get_psth_peaks(psth, nid):
    """Extract peaks from PSTH"""
    baseline = np.median(psth)
    thresh = baseline + MINTHRESH # peak detection threshold
    #thresh = max(MINTHRESH, BASELINEX*baseline) # peak detection threshold

    # find all local peaks above thresh:
    allpeakis, = argrelextrema(psth, np.greater_equal) # indices of all local maxima in psth
    threshis, = np.where(psth >= thresh) # indices of all thresh exceeding points in psth
    peakis = np.intersect1d(threshis, allpeakis) # indices of thresh exceeding maxima
    '''
    # alternatively, find only those local peaks above thresh that are separated from each
    # other by at least one point below threshold:
    splitis = np.where(np.diff(threshis) > 1)[0] + 1 # indices of thresh exceeding ranges
    splitthreshis = np.array_split(threshis, splitis) # list of arrays of contiguous indices
    peakis = []
    for threshis in splitthreshis: # for each contiguous thresh-exceeding range of points
        localpeakis = np.intersect1d(threshis, allpeakis) # indices of thresh exceeding maxima
        if len(localpeakis) == 0:
            continue # skip this range of points
        # keep only biggest maximum in this contiguous range:
        peakii = np.argmax(psth[localpeakis])
        peaki = localpeakis[peakii]
        peakis.append(peaki)
    peakis = np.asarray(peakis)
    '''
    ## TODO: or as alternative to above, divide peaks up according to where psth falls to
    ## baseline instead of to threshold. This may remove need for fwhm overlap test below,
    ## which is rather complicated to explain. Within each range above threshold, if there's
    ## a peak there, use the earliest value for li and the latest value for ri?

    ## TODO: maybe exclude peaks wider than some threshold, like 150 ms, as invalid

    if len(peakis) == 0:
        print('x%d' % nid, end='')
        return
    else:
        print("n%d" % nid, end='')


    # sort peakis in decreasing order of peak amplitude:
    peakis = peakis[psth[peakis].argsort()[::-1]]

    # collect left and right edges of all peaks:
    lis, ris, rmpeakiis = [], [], []
    for peakii, peaki in enumerate(peakis):
        ## TODO: maybe FWHM should be taken relative to baseline instead of 0.
        ## Or, is it more fair to measure FWHM relative to 0, not to some baseline that's
        ## different for each cell?
        try:
            li, ri = argfwhm(psth, peaki, fraction=FWFRACTION) # left and right indices
        except ValueError: # peaki has no FWHM
            rmpeakiis.append(peakii) # mark peaki for removal
            continue # skip to next peaki
        lis.append(li)
        ris.append(ri)
    lis = np.asarray(lis)
    ris = np.asarray(ris)
    peakis = np.delete(peakis, rmpeakiis) # delete peaks marked for removal

    # For each peak, check if any others overlap with it in width. If so, discard them.
    # This rejects smaller nearby spurious peaks:
    skippeakis, keeppeakis = [], []
    keeplis, keepris = [], []
    fwhms = []
    for peaki, li, ri in zip(peakis, lis, ris):
        if peaki in skippeakis:
            continue # was trumped by a larger peak with overlapping width, skip it
        # skip peaks which fall within width of this peak:
        #skippeakiis = (peakis >= li) & (peakis <= ri) & (peakis != peaki) # bool array
        # skip peaks whose right edge overlaps with this peak's left edge, and whose
        # left edge overlaps with this peak's right edge
        skippeakiis = ((li <= ris) & (lis <= ri)) & (peakis != peaki)
        skippeakis.extend(peakis[skippeakiis])
        fwhm = (ri - li)
        fwhms.append(fwhm)
        keeppeakis.append(peaki)
        keeplis.append(li)
        keepris.append(ri)
        print('.', end='')

    #skippeakis = np.asarray(skippeakis)
    fwhms = np.asarray(fwhms) * TRES * 1000 # sec
    peakis = np.asarray(keeppeakis)
    lis = np.asarray(keeplis)
    ris = np.asarray(keepris)

    return psth, thresh, baseline, peakis, lis, ris, fwhms

## TODO: use only neurons qualitatively deemed responsive?
# get active or all neuron ids for each section of r08:
#ssnids, recsecnids = get_ssnids(ptc22tr1r08s, strangesr08s, kind=NIDSKIND)
# get active or all neuron ids for each section of both r08 and r10:
ssnids, recsecnids = get_ssnids(recs, stranges, kind=NIDSKIND)

# calculate PSTHs for both sections of r08:
#midbinss = []
psthss = []
for rec, nids, strange in zip(recs, recsecnids, stranges):
    midbins, psths = rec.traster(nids=nids, natexps=False, strange=strange, plot=False,
                                 psth=True, binw=BINW, tres=TRES, norm='ntrials')
    psthss.append(psths)
    #midbinss.append(midbins)
    #figure()
    #pl.plot(midbins, psths.T, '-')

psthparamsrecsec = [] # params returned for each PSTH, for each recording section
fwhmsrecsec = [] # fwhm values, for each recording section
sparsrecsec = [] # sparseness values of cells with at least 1 peak, for each recording section
nidsrecsec = [] # nids of cells with at least 1 peak, for each recording section
psthsrecsec = [] # psths with at least 1 peak, for each recording section
for nids, psths, fmt in zip(recsecnids, psthss, ['b-', 'r-', 'r-', 'b-']):
    psthparams = {} # params returned for each PSTH by get_psth_peaks
    psthsfwhms = [] # fwhm values from all PSTHs from this recording section
    sparsenesses = [] # sparseness values for each PSTH in this recording section
    peaknids = [] # nids with peaks in this recording section
    peakpsths = [] # psths with peaks in this recording section
    for nid, psth in zip(nids, psths):
        #sparsenesses.append(sparseness(psth)) # for all PSTHS, even those without peaks
        psthparams[nid] = get_psth_peaks(psth, nid)
        if PLOTPSTH:
            plot_psth(psthparams, fmt)
        if psthparams[nid] == None:
            continue # this psth has no peaks, skip it
        psthsfwhms.append(psthparams[nid][-1])
        sparsenesses.append(sparseness(psth)) # for only those PSTHS with peaks
        peaknids.append(nid) # this nid had at least one peak
        peakpsths.append(psth)
    psthparamsrecsec.append(psthparams)
    fwhmsrecsec.append(np.hstack(psthsfwhms))
    sparsrecsec.append(np.hstack(sparsenesses))
    nidsrecsec.append(np.hstack(peaknids))
    psthsrecsec.append(np.hstack(peakpsths))
    print('\n') # two newlines

fwhms = [np.hstack([fwhmsrecsec[0], fwhmsrecsec[3]]), # desynched
         np.hstack([fwhmsrecsec[1], fwhmsrecsec[2]])] # synched

spars = [np.hstack([sparsrecsec[0], sparsrecsec[3]]), # desynched
         np.hstack([sparsrecsec[1], sparsrecsec[2]])] # synched

 
# plot FWHM distributions:
ticks = np.arange(FWHMMIN, FWHMMAX, XTICKSTEP)
bins = np.arange(FWHMMIN, FWHMMAX+FWHMSTEP, FWHMSTEP)
figure(figsize=figsize)
## TODO: try plotting distributions on log scale
n1 = hist(fwhms[1], bins=bins, color='r')[0] # second section in r08, synched
n0 = hist(fwhms[0], bins=bins, color='b')[0] # first section in r08, desynched
n = np.hstack([n0, n1])
xlim(xmin=FWHMMIN, xmax=FWHMMAX)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('FWHM (ms)')
ylabel('PSTH peak count')
gcfm().window.setWindowTitle('PSTH FWHM ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

# plot sparseness distributions:
ticks = np.arange(0, 1, 0.25)
bins = np.arange(0, 1+SPARSTEP, SPARSTEP)
figure(figsize=figsize)
n1 = hist(spars[1], bins=bins, color='r')[0] # second section in r08, synched
n0 = hist(spars[0], bins=bins, color='b')[0] # first section in r08, desynched
n = np.hstack([n0, n1])
xlim(xmin=0, xmax=1)
ylim(ymax=n.max()) # effectively normalizes the histogram
xticks(ticks)
xlabel('sparseness')
ylabel('PSTH count')
gcfm().window.setWindowTitle('PSTH sparseness ptc22.tr1.r08 ptc22.tr1.r10')
tight_layout(pad=0.3)

show()
