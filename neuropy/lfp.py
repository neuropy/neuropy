"""Defines the LFP class"""

from __future__ import division
from __future__ import print_function

import numpy as np

import matplotlib as mpl
import pylab as pl
from pylab import get_current_fig_manager as gcfm
from matplotlib.collections import LineCollection

from core import intround, issorted, iterable, lastcmd, split_tranges, tolist
import filter


class LFP(object):
    """Holds LFP data loaded from a numpy .npz-compatible .lfp.zip file"""
    def __init__(self, recording, fname):
        """
        self.chanpos: array of (x, y) LFP channel positions on probe, in order of
                      increasing zero-based channel IDs
        self.chans: channel IDs of rows in self.chanpos and self.data, in vertical
                    spatial order
        self.data: LFP voltage values, channels in rows (in vertical spatial order),
                   timepoints in columns
        self.t0: time in us of first LFP timepoint, from start of recording acquisition
        self.t1: time in us of last LFP timepoint
        self.tres: temporal resolution in us of each LFP timepoint
        self.uVperAd: number of uV per AD voltage value in LFP data
        """
        self.r = recording
        self.fname = fname # with full path

    def load(self):
        with open(self.fname, 'rb') as f:
            d = np.load(f)
            assert sorted(d.keys()) == ['chanpos', 'chans', 'data', 't0', 't1', 'tres',
                                        'uVperAD']
            # bind arrays in .lfp.zip file to self:
            for key, val in d.iteritems():
                # pull some singleton vals out of their arrays:
                if key in ['t0', 't1', 'tres']: # should all be us
                    val = int(val)
                elif key == 'uVperAD':
                    val = float(val)
                self.__setattr__(key, val)
        # make sure chans are in vertical spatial order:
        assert issorted(self.chanpos[self.chans][1])
        self.sampfreq = intround(1e6 / self.tres) # in Hz
        assert self.sampfreq == 1000 # should be 1000 Hz
        self.data = self.data * self.uVperAD # convert to float uV
        self.UV2UM = 0.05 # transforms LFP voltage in uV to position in um

    def save(self):
        ## TODO: option to overwrite original .lfp.zip file from spyke with filtered data,
        ## add filteredfreqs and filteredbws keywords when resaving to indicate what exactly
        ## was filtered out. Also, convert data back to int16?
        raise NotImplementedError

    def get_data(self):
        """Return data, testing first to see if it's been loaded"""
        try:
            self.data
        except AttributeError:
            self.load()
        return self.data

    def get_tssec(self):
        """Return full set of timestamps, in sec"""
        return np.arange(self.t0/1e6, self.t1/1e6, self.tres/1e6)

    def plot(self, t0=None, t1=None, chanis=None, gain=1, c='k', alpha=1.0, yunits='um',
             title=True, xlabel=True, figsize=(20, 6.5)):
        """Plot chanis of LFP data between t0 and t1 in sec. Unfortunatley, setting an alpha <
        1 doesn't seem to reveal detail when a line obscures itself, such as when plotting a
        very long time series"""
        self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        if t0 == None:
            t0, t1 = ts[0], ts[-1]
        if t1 == None:
            t1 = t0 + 10 # 10 sec window
        if chanis == None:
            chanis = range(len(self.chans)) # all chans
        t0i, t1i = ts.searchsorted((t0, t1))
        ts = ts[t0i:t1i] # constrained set of timestamps, in sec
        chanis = tolist(chanis)
        nchans = len(chanis)
        # grab desired channels and time range, convert uV to um:
        totalgain = self.UV2UM * gain
        data = self.data[chanis][:, t0i:t1i] * totalgain
        nt = len(ts)
        assert nt == data.shape[1]
        x = np.tile(ts, nchans)
        x.shape = nchans, nt
        segments = np.zeros((nchans, nt, 2)) # x vals in col 0, yvals in col 1
        segments[:, :, 0] = x
        segments[:, :, 1] = -data # set to -ve here because of invert_yaxis() below
        # add y offsets:
        for chanii, chani in enumerate(chanis):
            chan = self.chans[chani]
            ypos = self.chanpos[chan][1] # in um
            segments[chanii, :, 1] += ypos # vertical distance below top of probe
        if yunits == 'mm':
            segments[:, :, 1] /= 1000
        lc = LineCollection(segments, linewidth=1, linestyle='-', colors=c, alpha=alpha,
                            antialiased=True, visible=True)
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.add_collection(lc) # add to axes' pool of LCs
        a.autoscale(enable=True, tight=True)
        a.invert_yaxis()
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        if xlabel:
            a.set_xlabel("time (sec)")
        if yunits == 'um':
            a.set_ylabel("depth ($\mu$m)")
        elif yunits == 'mm':
            a.set_ylabel("depth (mm)")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
            a.text(0.998, 0.99, '%s' % self.r.name, transform=a.transAxes,
                   horizontalalignment='right', verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def std(self, t0=None, t1=None, chani=-1, width=None, tres=None, fmt='k-', title=True,
            figsize=(20, 3.5)):
        """Plot standard deviation of LFP signal from t0 to t1 on chani, using bins of width
        and tres"""
        uns = get_ipython().user_ns
        self.get_data()
        data = self.data[chani]
        ts = self.get_tssec()
        if t0 == None:
            t0 = ts[0]
        if t1 == None:
            t1 = ts[-1]
        if width == None:
            width = uns['LFPSIWIDTH'] # sec
        if tres == None:
            tres = uns['LFPSITRES'] # sec
        tranges = split_tranges([(t0, t1)], width, tres)
        stds = []
        for trange in tranges:
            ti0, ti1 = ts.searchsorted(trange)
            stds.append(data[ti0:ti1].std())
        stds = np.hstack(stds)
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.plot(tranges[:, 0], stds, fmt)
        a.autoscale(enable=True, tight=True)
        a.set_xlim(xmin=0) # acquisition starts at t=0
        a.set_xlabel('time (sec)')
        a.set_ylabel('LFP $\sigma$ ($\mu$V)')
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
        f.tight_layout(pad=0.3)
        self.f = f
        return stds

    def psd(self, t0=None, t1=None, f0=0.2, f1=110, p0=None, p1=None, chanis=-1,
            width=None, tres=None, xscale='log', figsize=(5, 5)):
        """Plot power spectral density from t0 to t1 in sec, from f0 to f1 in Hz, and clip
        power values from p0 to p1 in dB, based on channel index chani of LFP data. chanis=0
        uses most superficial channel, chanis=-1 uses deepest channel. If len(chanis) > 1,
        take mean of specified chanis. width and tres are in sec."""
        uns = get_ipython().user_ns
        self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        if t0 == None:
            t0, t1 = ts[0], ts[-1] # full duration
        if t1 == None:
            t1 = t0 + 10 # 10 sec window
        if width == None:
            width = uns['LFPWIDTH'] # sec
        if tres == None:
            tres = uns['LFPTRES'] # sec
        assert tres <= width
        NFFT = intround(width * self.sampfreq)
        noverlap = intround(NFFT - tres * self.sampfreq)
        t0i, t1i = ts.searchsorted((t0, t1))
        #ts = ts[t0i:t1i] # constrained set of timestamps, in sec
        data = self.data[:, t0i:t1i] # slice data
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        if iterable(chanis):
            data = data[chanis].mean(axis=0) # take mean of data on chanis
        else:
            data = data[chanis] # get single row of data at chanis
        #data = filter.notch(data)[0] # remove 60 Hz mains noise
        # convert data from uV to mV. I think P is in mV^2?:
        P, freqs = mpl.mlab.psd(data/1e3, NFFT=NFFT, Fs=self.sampfreq, noverlap=noverlap)
        # keep only freqs between f0 and f1:
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
        # for better visualization, clip power values to within (p0, p1) dB
        if p0 != None:
            P[P < p0] = p0
        if p1 != None:
            P[P > p1] = p1
        #self.P = P
        a.plot(freqs, P, 'k-')
        # add SI frequency band limits:
        LFPSILOWBAND, LFPSIHIGHBAND = uns['LFPSILOWBAND'], uns['LFPSIHIGHBAND']
        a.axvline(x=LFPSILOWBAND[0], c='r', ls='--')
        a.axvline(x=LFPSILOWBAND[1], c='r', ls='--')
        a.axvline(x=LFPSIHIGHBAND[0], c='b', ls='--')
        a.axvline(x=LFPSIHIGHBAND[1], c='b', ls='--')
        a.axis('tight')
        a.set_xscale(xscale)
        a.set_xlabel("frequency (Hz)")
        a.set_ylabel("power (dB)")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        a.text(0.998, 0.99, '%s' % self.r.name, color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return P, freqs
        
    def specgram(self, t0=None, t1=None, f0=0.1, f1=100, p0=-60, p1=None, chanis=-1,
                 width=None, tres=None, cm=None, colorbar=False, title=True,
                 figsize=(20, 6.5)):
        """Plot a spectrogram from t0 to t1 in sec, from f0 to f1 in Hz, and clip power values
        from p0 to p1 in dB, based on channel index chani of LFP data. chanis=0 uses most
        superficial channel, chanis=-1 uses deepest channel. If len(chanis) > 1, take mean of
        specified chanis. width and tres are in sec. As an alternative to cm.jet (the
        default), cm.gray, cm.hsv cm.terrain, and cm.cubehelix_r colormaps seem to bring out
        the most structure in the spectrogram"""
        uns = get_ipython().user_ns
        self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        if t0 == None:
            t0, t1 = ts[0], ts[-1] # full duration
        if t1 == None:
            t1 = t0 + 10 # 10 sec window
        if width == None:
            width = uns['LFPWIDTH'] # sec
        if tres == None:
            tres = uns['LFPTRES'] # sec
        assert tres <= width
        NFFT = intround(width * self.sampfreq)
        noverlap = intround(NFFT - tres * self.sampfreq)
        t0i, t1i = ts.searchsorted((t0, t1))
        #ts = ts[t0i:t1i] # constrained set of timestamps, in sec
        data = self.data[:, t0i:t1i] # slice data
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        if iterable(chanis):
            data = data[chanis].mean(axis=0) # take mean of data on chanis
        else:
            data = data[chanis] # get single row of data at chanis
        #data = filter.notch(data)[0] # remove 60 Hz mains noise
        # convert data from uV to mV, returned t is midpoints of time bins in sec from
        # start of data. I think P is in mV^2?:
        P, freqs, t = mpl.mlab.specgram(data/1e3, NFFT=NFFT, Fs=self.sampfreq,
                                        noverlap=noverlap)
        # convert t to time from start of acquisition:
        t += t0
        # keep only freqs between f0 and f1:
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
        # for better visualization, clip power values to within (p0, p1) dB
        if p0 != None:
            P[P < p0] = p0
        if p1 != None:
            P[P > p1] = p1
        #self.P = P
        # Label far left, right, top and bottom edges of imshow image. imshow interpolates
        # between these to place the axes ticks. Time limits are
        # set from midpoints of specgram time bins
        extent = t[0], t[-1], freqs[0], freqs[-1]
        #print('specgram extent: %r' % (extent,))
        # flip P vertically for compatibility with imshow:
        im = a.imshow(P[::-1], extent=extent, cmap=cm)
        a.autoscale(enable=True, tight=True)
        a.axis('tight')
        a.set_xlim(xmin=0) # acquisition starts at t=0
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        a.set_xlabel("time (sec)")
        a.set_ylabel("frequency (Hz)")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
            a.text(0.998, 0.99, '%s' % self.r.name, color='w', transform=a.transAxes,
                   horizontalalignment='right', verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        if colorbar:
            f.colorbar(im, pad=0) # creates big whitespace to the right for some reason
        self.f = f
        return P, freqs, t

    def notch(self, chanis=None, freq=60, bw=0.25, gpass=0.01, gstop=30, ftype='ellip'):
        """Filter out frequencies centered on freq (Hz), of bandwidth +/- bw (Hz) on
        data row indices chanis, in-place.

        ftype: 'ellip', 'butter', 'cheby1', 'cheby2', 'bessel'
        """
        data = self.get_data()
        if chanis == None:
            chanis = np.arange(len(data))
        data = data[chanis]
        data, b, a = filter.notch(data, self.sampfreq, freq, bw, gpass, gstop, ftype)
        self.data[chanis] = data
        return b, a

    def naivenotch(self, freqs=60, bws=1):
        """Filter out frequencies in data centered on freqs (Hz), of bandwidths bws (Hz),
        in-place. Filtering out by setting components to 0 is probably naive"""
        data = self.get_data()
        self.data = filter.naivenotch(data, self.sampfreq, freqs, bws)

    def filter(self, chanis=None, f0=0, f1=7, fr=0.5, gpass=0.01, gstop=30, ftype='ellip'):
        """Bandpass filter data on row indices chanis, between f0 and f1 (Hz), with filter
        rolloff (?) fr (Hz). Done in-place.

        ftype: 'ellip', 'butter', 'cheby1', 'cheby2', 'bessel'
        """
        data = self.get_data()
        if chanis == None:
            chanis = np.arange(len(data))
        data = data[chanis]
        data, b, a = filter.filter(data, self.sampfreq, f0, f1, fr, gpass, gstop, ftype)
        self.data[chanis] = data
        return b, a

    def filterord(self, chanis=None, f0=300, f1=None, order=4, rp=None, rs=None,
                  btype='highpass', ftype='butter'):
        """Bandpass filter data in-place by specifying filter order and btype, instead of
        gpass and gstop"""
        data = self.get_data()
        if chanis == None:
            chanis = np.arange(len(data))
        data = data[chanis]
        data, b, a = filter.filterord(data, self.sampfreq, f0, f1, order, rp, rs, btype, ftype)
        self.data[chanis] = data
        return b, a

    def si(self, kind=None, chani=-1, width=None, tres=None,
           lfpwidth=None, lfptres=None, lowband=None, highband=None, plot=True,
           showxlabel=True, showylabel=True, showtitle=True, showtext=True,
           figsize=(20, 3.5), swapaxes=False):
        """Calculate an LFP synchrony index, using potentially overlapping windows of width
        and tres, in sec, from the LFP spectrogram, itself composed of bins of lfpwidth and
        lfptres. Note that for power ratio methods (kind: L/(L+H) or L/H), width and tres are
        not used, only lfpwidth and lfptres. Options for kind are:

        'L/(L+H)': fraction of power in low band vs total power (Saleem2010)

        'L/H': low to highband power ratio (Li, Poo, Dan 2009)

        'cv': coefficient of variation (std / mean) of all power

        'ncv': normalized CV: (std - mean) / (std + mean)

        'nstdmed': normalized stdmed: (std - med) / (std + med)

        'n2stdmean': normalized 2stdmean: (2*std - mean) / (2*std + mean)

        'n3stdmean': normalized 3stdmean: (3*std - mean) / (3*std + mean)

        """
        uns = get_ipython().user_ns
        if kind == None:
            kind = uns['LFPSIKIND']
        if kind.startswith('L/'):
            pratio = True
        else:
            pratio = False

        data = self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        x = data[chani] / 1e3 # convert from uV to mV
        x = filter.notch(x)[0] # remove 60 Hz mains noise
        try:
            rr = self.r.e0.I['REFRESHRATE']
        except AttributeError: # probably a recording with no experiment
            rr = 200 # assume 200 Hz refresh rate
        if rr <= 100: # CRT was at low vertical refresh rate
            print('filtering out %d Hz from LFP in %s' % (intround(rr), self.r.name))
            x = filter.notch(x, freq=rr)[0] # remove CRT interference

        if width == None:
            width = uns['LFPSIWIDTH'] # sec
        if tres == None:
            tres = uns['LFPSITRES'] # sec
        if lfpwidth == None:
            lfpwidth = uns['LFPWIDTH'] # sec
        if lfptres == None:
            lfptres = uns['LFPTRES'] # sec
        if lowband == None:
            lowband = uns['LFPSILOWBAND']
        f0, f1 = lowband
        if highband == None:
            highband = uns['LFPSIHIGHBAND']
        f2, f3 = highband

        assert lfptres <= lfpwidth
        NFFT = intround(lfpwidth * self.sampfreq)
        noverlap = intround(NFFT - lfptres * self.sampfreq)
        #print('len(x), NFFT, noverlap: %d, %d, %d' % (len(x), NFFT, noverlap))
        # t is midpoints of timebins in sec from start of data. P is in mV^2?:
        P, freqs, Pt = mpl.mlab.specgram(x, NFFT=NFFT, Fs=self.sampfreq, noverlap=noverlap)
        # don't convert power to dB, just washes out the signal in the ratio:
        #P = 10. * np.log10(P)
        # convert t to time from start of acquisition:
        Pt += ts[0]
        nfreqs = len(freqs)

        # keep only freqs between f0 and f1, and f2 and f3:
        f0i, f1i, f2i, f3i = freqs.searchsorted([f0, f1, f2, f3])
        lP = P[f0i:f1i] # nsubfreqs x nt
        hP = P[f2i:f3i] # nsubfreqs x nt
        lP = lP.sum(axis=0) # nt
        hP = hP.sum(axis=0) # nt

        if pratio:
            t = Pt
            ylim = 0, 1
            ylabel = 'SI (%s)' % kind
        else:
            # potentially overlapping bin time ranges:
            trange = Pt[0], Pt[-1]
            tranges = split_tranges([trange], width, tres) # in sec
            ntranges = len(tranges)
            tis = Pt.searchsorted(tranges) # ntranges x 2 array
            # number of timepoints to use for each trange, almost all will be the same width:
            binnt = intround((tis[:, 1] - tis[:, 0]).mean())
            binhP = np.zeros((ntranges, binnt)) # init appropriate array
            for trangei, t0i in enumerate(tis[:, 0]):
                binhP[trangei] = hP[t0i:t0i+binnt]
            # get midpoint of each trange:
            t = tranges.mean(axis=1)

        #old_settings = np.seterr(all='ignore') # suppress div by 0 errors
        # plot power signal to be analyzed
        #self.si_plot(Pt, hP, t0=0, t1=t[-1], ylim=None, ylabel='highband power',
        #             title=lastcmd()+' highband power', text=self.r.name)
        hlines = []
        if kind[0] == 'n':
            ylim = -1, 1
            hlines = [0]
        # calculate some metric of each column, ie each width:
        if kind == 'L/(L+H)':
            si = lP/(lP + hP)
        elif kind == 'L/H':
            si = lP/hP
        elif kind == 'nLH':
            t = Pt
            si = (lP - hP) / (lP + hP)
            ylabel = 'LFP (L - H) / (L + H)'
        elif kind == 'cv':
            si = binhP.std(axis=1) / binhP.mean(axis=1)
            ylim = 0, 2
            ylabel = 'LFP power CV'
        elif kind == 'ncv':
            s = binhP.std(axis=1)
            mean = binhP.mean(axis=1)
            si = (s - mean) / (s + mean)
            ylabel = 'LFP power (std - mean) / (std + mean)'
            #pl.plot(t, s)
            #pl.plot(t, mean)
        elif kind == 'n2stdmean':
            s2 = 2 * binhP.std(axis=1)
            mean = binhP.mean(axis=1)
            si = (s2 - mean) / (s2 + mean)
            ylabel = 'LFP power (2*std - mean) / (2*std + mean)'
            hlines = [-0.1, 0, 0.1] # demarcate desynched and synched thresholds
            #pl.plot(t, s2)
            #pl.plot(t, mean)
        elif kind == 'n3stdmean':
            s3 = 3 * binhP.std(axis=1)
            mean = binhP.mean(axis=1)
            si = (s3 - mean) / (s3 + mean)
            ylabel = 'LFP power (3*std - mean) / (3*std + mean)'
            hlines = [-0.1, 0, 0.1] # demarcate desynched and synched thresholds
            #pl.plot(t, s3)
            #pl.plot(t, mean)
        elif kind == 'n4stdmean':
            s4 = 4 * binhP.std(axis=1)
            mean = binhP.mean(axis=1)
            si = (s4 - mean) / (s4 + mean)
            ylabel = 'LFP power (4*std - mean) / (4*std + mean)'
            #pl.plot(t, s4)
            #pl.plot(t, mean)
        elif kind == 'nstdmed':
            s = binhP.std(axis=1)
            med = np.median(binhP, axis=1)
            si = (s - med) / (s + med)
            ylabel = 'LFP power (std - med) / (std + med)'
            hlines = [-0.1, 0, 0.1] # demarcate desynched and synched thresholds
            #pl.plot(t, s)
            #pl.plot(t, med)
        elif kind == 'n2stdmed':
            s2 = 2 * binhP.std(axis=1)
            med = np.median(binhP, axis=1)
            si = (s2 - med) / (s2 + med)
            ylabel = 'LFP power (2*std - med) / (2*std + med)'
            hlines = [-0.1, 0, 0.1] # demarcate desynched and synched thresholds
            #pl.plot(t, s2)
            #pl.plot(t, med)
        elif kind == 'n3stdmed':
            s3 = 3 * binhP.std(axis=1)
            med = np.median(binhP, axis=1)
            si = (s3 - med) / (s3 + med)
            ylabel = 'LFP power (3*std - med) / (3*std + med)'
            hlines = [-0.1, 0, 0.1] # demarcate desynched and synched thresholds
            #pl.plot(t, s3)
            #pl.plot(t, med)
        elif kind == 'nstdmin':
            s = binhP.std(axis=1)
            min = binhP.min(axis=1)
            si = (s - min) / (s + min)
            ylabel = 'LFP power (std - min) / (std + min)'
            #pl.plot(t, s)
            #pl.plot(t, min)
        elif kind == 'nmadmean':
            mean = binhP.mean(axis=1)
            mad = (np.abs(binhP - mean[:, None])).mean(axis=1)
            si = (mad - mean) / (mad + mean)
            ylabel = 'MUA (MAD - mean) / (MAD + mean)'
            #pl.plot(t, mad)
            #pl.plot(t, mean)
        elif kind == 'nmadmed':
            med = np.median(binhP, axis=1)
            mad = (np.abs(binhP - med[:, None])).mean(axis=1)
            si = (mad - med) / (mad + med)
            ylabel = 'MUA (MAD - median) / (MAD + median)'
            #pl.plot(t, mad)
            #pl.plot(t, med)
        elif kind == 'nvarmin':
            v = binhP.var(axis=1)
            min = binhP.min(axis=1)
            si = (v - min) / (v + min)
            ylabel = 'LFP power (std - min) / (std + min)'
            #pl.plot(t, v)
            #pl.plot(t, min)
        elif kind == 'nptpmean':
            ptp = binhP.ptp(axis=1)
            mean = binhP.mean(axis=1)
            si = (ptp - mean) / (ptp + mean)
            ylabel = 'MUA (ptp - mean) / (ptp + mean)'
            #pl.plot(t, ptp)
            #pl.plot(t, mean)
        elif kind == 'nptpmed':
            ptp = binhP.ptp(axis=1)
            med = np.median(binhP, axis=1)
            si = (ptp - med) / (ptp + med)
            ylabel = 'MUA (ptp - med) / (ptp + med)'
            #pl.plot(t, ptp)
            #pl.plot(t, med)
        elif kind == 'nptpmin':
            ptp = binhP.ptp(axis=1)
            min = binhP.min(axis=1)
            si = (ptp - min) / (ptp + min)
            ylabel = 'MUA (ptp - min) / (ptp + min)'
            #pl.plot(t, ptp)
            #pl.plot(t, min)
        elif kind == 'nmaxmin':
            max = binhP.max(axis=1)
            min = binhP.min(axis=1)
            si = (max - min) / (max + min)
            ylabel = 'MUA (max - min) / (max + min)'
            #pl.plot(t, max)
            #pl.plot(t, min)
        else:
            raise ValueError('unknown kind %r' % kind)
        if plot:
            self.si_plot(t, si, t0=0, t1=t[-1], ylim=ylim, showxlabel=showxlabel,
                         showylabel=showylabel, ylabel=ylabel, showtitle=showtitle,
                         title=lastcmd(), showtext=showtext, text=self.r.name, hlines=hlines,
                         figsize=figsize, swapaxes=swapaxes)
        #np.seterr(**old_settings) # restore old settings
        return si, t # t are midpoints of bins, from start of acquisition
    '''
    def si_hilbert(self, chani=-1, lowband=None, highband=None, ratio='L/(L+H)',
                   plot=True):
        """Return synchrony index, i.e. power ratio of low vs high bands, as measured by
        Hilbert transform (Saleem2010). Use either L/(L+H) ratio (Saleem2010) or L/H ratio
        (Li, Poo, Dan 2009)"""
        if lowband == None:
            lowband = uns['LFPSILOWBAND']
        f0, f1 = lowband
        if highband == None:
            highband = uns['LFPSIHIGHBAND']
        f2, f3 = highband
        data = self.get_data()
        t = self.get_tssec() # full set of timestamps, in sec
        t0, t1 = t[0], t[-1] # full duration
        x = data[chani] / 1e3 # convert from uV to mV
        x = filter.notch(x)[0] # remove 60 Hz mains noise
        rr = self.r.e0.I['REFRESHRATE']
        if rr <= 100: # CRT was at low vertical refresh rate
            print('filtering out %d Hz from LFP in %s' % (intround(rr), self.r.name))
            x = filter.notch(x, freq=rr)[0] # remove CRT interference
        # remove everything below f0:
        #x = filter.filterord(data=x, f0=f0, order=4, btype='highpass')[0]
        #x = filter.filter(data=x, f0=0.5, f1=0, fr=0.1, ftype='ellip', gstop=20)
        # remove everything above f3:
        x = filter.filterord(data=x, f0=f3, order=10, btype='lowpass')[0]
        l = filter.filterord(data=x, f0=f1, order=4, btype='lowpass')[0]
        h = filter.filterord(data=x, f0=f2, order=11, btype='highpass')[0]
        lP, lPh, lE, lA = filter.hilbert(l)
        hP, hPh, hE, hA = filter.hilbert(h)

        if ratio == 'L/(L+H)':
            r = lP/(hP + lP)
        elif ratio == 'L/H':
            r = lP/hP
        else:
            raise ValueError
        if plot:
            ylabel = 'LFP synchrony index (%s)' % ratio
            self.si_plot(t, r, t0, t1, ylabel, title=lastcmd(), text=self.r.name)
        return r, t
    '''
    def si_plot(self, t, si, t0=None, t1=None, ylim=None, ylabel=None,
                showxlabel=True, showylabel=True, showtitle=True,
                title=None, showtext=True, text=None, hlines=[0], figsize=(20, 6.5),
                swapaxes=False):
        """Plot synchrony index as a function of time, with hopefully the same
        temporal scale as some of the other plots in self"""
        if figsize == None:
            f = pl.gcf()
            a = pl.gca()
        else:
            f = pl.figure(figsize=figsize)
            a = f.add_subplot(111)

        xlim = t0, t1
        xlabel = "time (sec)"
        if ylabel == None:
            ylabel = "synchrony index (AU?)"

        if swapaxes:
            t, si = si, t # swap t and si
            xlim, ylim = ylim, xlim
            ylim = ylim[1], ylim[0] # swap new ylimits so t=0 is at top
            xlabel, ylabel = ylabel, xlabel # swap labels
            showxlabel, showylabel = showylabel, showxlabel # swap flags
            # underplot vertical lines:
            for hline in hlines:
                a.axvline(x=hline, c='e', ls='--', marker=None)
        else:
            # underplot horizontal lines:
            for hline in hlines:
                a.axhline(y=hline, c='e', ls='--', marker=None)

        a.plot(t, si, 'k-')
        a.set_xlim(xlim) # low/high limits are unchanged if None
        a.set_ylim(ylim)
        if showxlabel:
            a.set_xlabel(xlabel)
        if showylabel:
            a.set_ylabel(ylabel)
        #a.autoscale(axis='x', enable=True, tight=True)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        if title:
            gcfm().window.setWindowTitle(title)
            if showtitle:
                a.set_title(title)
        if showtext:
            a.text(0.998, 0.01, '%s' % text, color='k', transform=a.transAxes,
                   horizontalalignment='right', verticalalignment='bottom')
        f.tight_layout(pad=0.3) # crop figure to contents

    def filterwavelet(self, chanis=None, wname="db4", maxlevel=6):
        """Filter data using wavelet multi-level decomposition and reconstruction (WMLDR).
        See Wiltschko2008"""
        data = self.get_data()
        if chanis == None:
            chanis = np.arange(len(data))
        data = data[chanis]
        data = filter.wavelet(data, wname, maxlevel)
        self.data[chanis] = data
