"""Calculate and plot power spectral density from the LFP of desired tracks/recordings"""

from __future__ import division, print_function

import numpy as np
import matplotlib as mpl
import pylab as pl
from pylab import get_current_fig_manager as gcfm

from core import intround
import filter

#tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
tracks = [ptc15.tr7c, ptc17.tr1, ptc17.tr2b, ptc18.tr1, ptc18.tr2c, ptc20.tr1, ptc20.tr2,
          ptc20.tr3, ptc21.tr2, ptc21.tr5c, ptc21.tr6b, ptc22.tr1, ptc22.tr2, ptc22.tr4b,
          ptc22.tr5b]

F0, F1 = 0.2, 110 # Hz
P0, P1 = None, None
chanis = -1
width, tres = 10, 5 # sec
figsize = (5, 5)
XSCALE = 'log'

uns = get_ipython().user_ns
if width == None: # window width
    width = uns['LFPWIDTH'] # sec
if tres == None: # window tres
    tres = uns['LFPTRES'] # sec
assert tres <= width

SAMPFREQ = 1000 # Hz, should be the same for all LFPs
NFFT = intround(width * SAMPFREQ)
NOVERLAP = intround(NFFT - tres * SAMPFREQ)

def plot_psd(data, titlestr):
    data = filter.notch(data)[0] # remove 60 Hz mains noise, as for SI calc
    # convert data from uV to mV. I think P is in mV^2?:
    P, freqs = mpl.mlab.psd(data/1e3, NFFT=NFFT, Fs=SAMPFREQ, noverlap=NOVERLAP)
    # keep only freqs between F0 and F1:
    f0, f1 = F0, F1 # need to set different local names, since they're not read-only
    if f0 == None:
        f0 = freqs[0]
    if f1 == None:
        f1 = freqs[-1]
    lo, hi = freqs.searchsorted([f0, f1])
    P, freqs = P[lo:hi], freqs[lo:hi]
    # check for and replace zero power values (ostensibly due to gaps in recording)
    # before attempting to convert to dB:
    zis = np.where(P == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        P[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float
        minnzval = P.min() # get minimum nonzero value
        P[zis] = minnzval # replace with min nonzero values
    P = 10. * np.log10(P) # convert power to dB wrt 1 mV^2?
    # for better visualization, clip power values to within (P0, P1) dB
    if P0 != None:
        P[P < P0] = P0
    if P1 != None:
        P[P > P1] = P1
    f = pl.figure(figsize=figsize)
    a = f.add_subplot(111)
    a.plot(freqs, P, 'k-')
    # add SI frequency band limits:
    a.axvline(x=LFPSILOWBAND[0], c='r', ls='--')
    a.axvline(x=LFPSILOWBAND[1], c='r', ls='--')
    a.axvline(x=LFPSIHIGHBAND[0], c='b', ls='--')
    a.axvline(x=LFPSIHIGHBAND[1], c='b', ls='--')
    a.axis('tight')
    a.set_xscale(XSCALE)
    a.set_ylim(ymin=P[-1]) # use last power value to set ymin
    a.set_xlabel("frequency (Hz)")
    a.set_ylabel("power (dB)")
    gcfm().window.setWindowTitle(titlestr+' '+XSCALE)
    f.tight_layout(pad=0.3) # crop figure to contents


# collect LFP data from all recordings in all tracks:
alldata = []
for track in tracks:
    print(track.absname)
    trackdata = []
    #for rec in recs:
    #for rid in BSRIDS[track.absname]:
    for rec in track.r.values():
        #rec = track.r[rid]
        lfp = rec.lfp
        data = lfp.get_data()
        assert lfp.sampfreq == SAMPFREQ
        if iterable(chanis):
            data = data[chanis].mean(axis=0) # take mean of data on chanis
        else:
            data = data[chanis] # get single row of data at chanis
        trackdata.append(data)
        # keeping all 10 chans for each LFP for all tracks gives MemoryError:
        del data
        del lfp.data
        print('.', end='')
    trackdata = np.hstack(trackdata)
    alldata.append(trackdata)
    del trackdata
    nMiB = sum([data.nbytes for data in alldata])/(1024**2)
    print("%g MiB" % nMiB)
    #plot_psd(trackdata, track.absname+' PSD')
alldata = np.hstack(alldata)

plot_psd(alldata, 'all tracks PSD')

pl.show()
