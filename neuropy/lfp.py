"""Defines the LFP class"""

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
        self.t0: time in us of first LFP timepoint, from start of ADC clock
        self.t1: time in us of last LFP timepoint
        self.tres: temporal resolution in us of each LFP timepoint
        self.uVperAd: number of uV per AD voltage value in LFP data
        """
        self.r = recording
        self.fname = fname # with full path

    def load(self):
        with open(self.fname, 'rb') as f:
            d = np.load(f)
            stdnames = ['chanpos', 'chans', 'data', 't0', 't1', 'tres', 'uVperAD']
            optnames = ['chan0', 'probename']
            # bind standard array names in .lfp.zip file to self:
            for key in stdnames:
                assert key in d
                val = d[key]
                # pull some singleton vals out of their arrays:
                if key in ['t0', 't1', 'tres']: # should all be us
                    val = int(val)
                elif key == 'uVperAD':
                    val = float(val)
                self.__setattr__(key, val)
            # bind optional array names in .lfp.zip file to self:
            for key in optnames:
                if key in d:
                    val = d[key]
                    # pull some singleton vals out of their arrays:
                    if key == 'chan0':
                        val = int(val)
                    elif key == 'probename':
                        val = val.tostring().decode() # convert from bytes to py3 unicode str
                    self.__setattr__(key, val)
        try:
            self.chan0
        except AttributeError: # try and figure out base of channel numbering
            nchans, nprobechans = len(self.chans), len(self.chanpos)
            if nchans < nprobechans:
                # it's probably from a .srf recording with only a subset of chans selected for
                # separate analog LFP recording
                self.chan0 = 0 # base of channel numbering, always 0-based for .srf recordings
                print("Found %d LFP channels, assuming 0-based channel numbering from .srf "
                      "recording" % nchans)
            elif nchans == nprobechans: # all probe channels have LFP
                self.chan0 = min(self.chans) # base of channel numbering
            else: # nchans > nprobechans
                raise ValueError("don't know how to handle nchans=%d > nprobechans=%d" %
                                 (nchans, nprobechans))
        assert self.chan0 in [0, 1] # either 0- or 1-based
        # make sure chans are in vertical spatial order:
        ypos = self.chanpos[self.chans - self.chan0][:, 1]
        if not issorted(ypos):
            print("LFP chans in %s aren't sorted by depth, sorting them now" % self.fname)
            sortis = ypos.argsort()
            self.chans = self.chans[sortis]
            self.data = self.data[sortis]
            newypos = self.chanpos[self.chans - self.chan0][:, 1]
            assert issorted(newypos)
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

    def get_ts(self):
        """Return full set of timestamps, in us"""
        return np.arange(self.t0, self.t1, self.tres)

    def get_tssec(self):
        """Return full set of timestamps, in s"""
        return self.get_ts() / 1e6

    def apply_lim2stim(self, t0, t1):
        """Limit t0 and t1 (in sec) to start and end of first and last trial"""
        assert len(self.r.e) == 1
        print("original trange: %r" % ((t0, t1),))
        ttranges = self.r.e0.ttranges
        t0 = max(t0, ttranges[0, 0]/1e6)
        t1 = min(t1, ttranges[-1, 1]/1e6)
        print("lim2stim trange: %r" % ((t0, t1),))
        return t0, t1

    def plot(self, t0=None, t1=None, chanis=None, gain=1, c='k', alpha=1.0, yunits='um',
             yticks=None, title=True, xlabel=True, relative2t0=False, lim2stim=False,
             scalebar=True, lw=4, figsize=(20, 6.5)):
        """Plot chanis of LFP data between t0 and t1 in sec. Unfortunatley, setting an alpha <
        1 doesn't seem to reveal detail when a line obscures itself, such as when plotting a
        very long time series. relative2t0 controls whether to plot relative to t0, or
        relative to start of ADC clock. lim2stim limits the time range only to when a stimulus
        was on screen, i.e. to the outermost times of non-NULL din. If only one chan is
        requested, it's plotted on a mV scale instead of a spatial scale."""
        self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        if t0 == None:
            t0, t1 = ts[0], ts[-1]
        if t1 == None:
            t1 = t0 + 10 # 10 sec window
        if chanis == None:
            chanis = range(len(self.chans)) # all chans
        if lim2stim:
            t0, t1 = self.apply_lim2stim(t0, t1)
        t0i, t1i = ts.searchsorted((t0, t1))
        ts = ts[t0i:t1i] # constrained set of timestamps, in sec
        chanis = tolist(chanis)
        nchans = len(chanis)
        # grab desired channels and time range:
        data = self.data[chanis][:, t0i:t1i]
        if nchans > 1: # convert uV to um:
            totalgain = self.UV2UM * gain
            data = data * totalgain
        else: # convert uV to mV:
            data = data / 1000
            yunits = 'mV'
        nt = len(ts)
        assert nt == data.shape[1]
        if relative2t0:
            # convert ts to time from t0, otherwise plot time from start of ADC clock:
            ts -= t0
        x = np.tile(ts, nchans)
        x.shape = nchans, nt
        segments = np.zeros((nchans, nt, 2)) # x vals in col 0, yvals in col 1
        segments[:, :, 0] = x
        if nchans > 1:
            segments[:, :, 1] = -data # set to -ve here because of invert_yaxis() below
        else:
            segments[:, :, 1] = data
        if nchans > 1: # add y offsets:
            maxypos = 0
            for chanii, chani in enumerate(chanis):
                chan = self.chans[chani]
                ypos = self.chanpos[chan][1] # in um
                segments[chanii, :, 1] += ypos # vertical distance below top of probe
                maxypos = max(maxypos, ypos)
            if yunits == 'mm': # convert from um to mm
                segments[:, :, 1] /= 1000
                maxypos = maxypos / 1000 # convert from int to float
                totalgain = totalgain / 1000
        lc = LineCollection(segments, linewidth=1, linestyle='-', colors=c, alpha=alpha,
                            antialiased=True, visible=True)
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.add_collection(lc) # add to axes' pool of LCs
        if scalebar: # add vertical scale bar at end of last channel to represent 1 mV:
            if nchans > 1:
                ymin, ymax = maxypos-500*totalgain, maxypos+500*totalgain # +/- 0.5 mV
            else:
                ymin, ymax = -0.5, 0.5 # mV
            a.vlines(ts.max()*0.99, ymin, ymax, lw=lw, colors='e')
        a.autoscale(enable=True, tight=True)
        # depending on relative2t0 above, x=0 represents either t0 or time ADC clock started:
        a.set_xlim(xmin=0)
        if nchans > 1:
            a.invert_yaxis() # for spatial scale
        if yticks != None:
            a.set_yticks(yticks)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        if xlabel:
            a.set_xlabel("time (s)")
        if yunits == 'um':
            a.set_ylabel("depth ($\mu$m)")
        elif yunits == 'mm':
            a.set_ylabel("depth (mm)")
        elif yunits == 'mV':
            a.set_ylabel("LFP (mV)")
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
        a.set_xlim(xmin=0) # ADC clock starts at t=0
        a.set_xlabel('time (s)')
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
            width = uns['LFPSPECGRAMWIDTH'] # sec
        if tres == None:
            tres = uns['LFPSPECGRAMTRES'] # sec
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
        LFPPRLOBAND, LFPPRHIBAND = uns['LFPPRLOBAND'], uns['LFPPRHIBAND']
        a.axvline(x=LFPPRLOBAND[0], c='r', ls='--')
        a.axvline(x=LFPPRLOBAND[1], c='r', ls='--')
        a.axvline(x=LFPPRHIBAND[0], c='b', ls='--')
        a.axvline(x=LFPPRHIBAND[1], c='b', ls='--')
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
                 width=None, tres=None, cm='jet', colorbar=False,
                 showstates=False, lw=4, alpha=1, relative2t0=False, lim2stim=False,
                 title=True, reclabel=True, swapaxes=False, figsize=None):
        """Plot a spectrogram from t0 to t1 in sec, from f0 to f1 in Hz, and clip power values
        from p0 to p1 in dB, based on channel index chani of LFP data. chanis=0 uses most
        superficial channel, chanis=-1 uses deepest channel. If len(chanis) > 1, take mean of
        specified chanis. width and tres are in sec. As an alternative to cm.jet (the
        default), cm.gray, cm.hsv cm.terrain, and cm.cubehelix_r colormaps seem to bring out
        the most structure in the spectrogram. showstates controls whether to plot lines
        demarcating desynchronized and synchronized periods. relative2t0 controls whether to
        plot relative to t0, or relative to start of ADC clock. lim2stim limits the time range
        only to when a stimulus was on screen, i.e. to the outermost times of non-NULL din"""
        uns = get_ipython().user_ns
        self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        if t0 == None:
            t0, t1 = ts[0], ts[-1] # full duration
        if t1 == None:
            t1 = t0 + 10 # 10 sec window
        if lim2stim:
            t0, t1 = self.apply_lim2stim(t0, t1)
        dt = t1 - t0
        if width == None:
            width = uns['LFPSPECGRAMWIDTH'] # sec
        if tres == None:
            tres = uns['LFPSPECGRAMTRES'] # sec
        assert tres <= width
        NFFT = intround(width * self.sampfreq)
        noverlap = intround(NFFT - tres * self.sampfreq)
        t0i, t1i = ts.searchsorted((t0, t1))
        #ts = ts[t0i:t1i] # constrained set of timestamps, in sec
        data = self.data[:, t0i:t1i] # slice data
        if figsize == None:
            # convert from recording duration time to width in inches, 0.87 accommodates
            # padding around the specgram:
            figwidth = (dt / 1000) * 5 + 0.87
            figheight = 2.5 # inches
            figsize = figwidth, figheight
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
        if not relative2t0:
            t += t0 # convert t to time from start of ADC clock:
        # keep only freqs between f0 and f1:
        if f0 == None:
            f0 = freqs[0]
        if f1 == None:
            f1 = freqs[-1]
        df = f1 - f0
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

        # plot horizontal bars over time demarcating different ranges of SI values,
        # or manually defined desynched and synched periods:
        statelinepos = f0 - df*0.015 # plot horizontal bars just below x axis
        if showstates:
            if showstates in [True, 'auto']:
                print("TODO: there's an offset plotting bug for 'auto', compare with 'manual'")
                si, t = self.si(plot=False)
                stranges, states = self.si_split(si, t) # sec
                STATECOLOURS = uns['LFPPRBINCOLOURS']
            elif showstates == 'manual':
                stranges, states = [], []
                for state in uns['MANUALSTATES']:
                    for strange in uns['REC2STATE2TRANGES'][self.r.absname][state]:
                        stranges.append(strange)
                        states.append(state)
                stranges = np.vstack(stranges) # 2D array
                STATECOLOURS = uns['MANUALSTATECOLOURS']
            else:
                raise ValueError('invalid value showstates=%r' % showstates)
            # clip stranges to t0, t1:
            stranges[0, 0] = max(stranges[0, 0], t0)
            stranges[-1, 1] = min(stranges[-1, 1], t1)
            if swapaxes:
                lines = a.vlines
            else:
                lines = a.hlines
            for strange, state in zip(stranges, states):
                clr = STATECOLOURS[state]
                lines(statelinepos, strange[0], strange[1], colors=clr, lw=lw, alpha=alpha,
                      clip_on=False)

        # Label far left, right, top and bottom edges of imshow image. imshow interpolates
        # between these to place the axes ticks. Time limits are
        # set from midpoints of specgram time bins
        extent = t[0], t[-1], freqs[0], freqs[-1]
        #print('specgram extent: %r' % (extent,))
        # flip P vertically for compatibility with imshow:
        im = a.imshow(P[::-1], extent=extent, cmap=cm)
        a.autoscale(enable=True, tight=True)
        a.axis('tight')
        # depending on relative2t0 above, x=0 represents either t0 or time ADC clock started:
        a.set_xlim(xmin=0, xmax=t[-1])
        a.set_ylim(ymin=freqs[0], ymax=freqs[-1])
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        a.set_xlabel("time (s)")
        a.set_ylabel("frequency (Hz)")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
        if reclabel:
            a.text(0.994, 0.95, '%s' % self.r.absname, color='w', transform=a.transAxes,
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
           lfpwidth=None, lfptres=None, loband=None, hiband=None, plot=True,
           showstates='auto', statelinepos=[0.2], lw=4, alpha=1, relative2t0=False,
           lim2stim=False, showxlabel=True, showylabel=True, showtitle=True, title=None,
           reclabel=True, swapaxes=False, figsize=None):
        """Calculate an LFP synchrony index, using potentially overlapping windows of width
        and tres, in sec, from the LFP spectrogram, itself composed of bins of lfpwidth and
        lfptres. relative2t0 controls whether to plot relative to t0, or relative to start of
        ADC clock. lim2stim limits the time range only to when a stimulus was presented, i.e.
        to the outermost times of non-NULL din.

        Note that for power ratio methods (kind: L/(L+H) or L/H),
        width and tres are not used, only lfpwidth and lfptres. Options for kind are:

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
        if kind in ['L/(L+H)', 'L/H', 'nLH']: # it's a power ratio measure
            pr = True
        else:
            pr = False

        data = self.get_data()
        ts = self.get_tssec() # full set of timestamps, in sec
        t0, t1 = ts[0], ts[-1]
        if lim2stim:
            t0, t1 = self.apply_lim2stim(t0, t1)
        dt = t1 - t0
        if figsize == None:
            # convert from recording duration time to width in inches, 0.87 accommodates
            # padding around the SI plot:
            figwidth = (dt / 1000) * 5 + 0.87
            figheight = 2.5 # inches
            figsize = figwidth, figheight

        t0i, t1i = ts.searchsorted((t0, t1))
        x = data[chani, t0i:t1i] / 1e3 # slice data, convert from uV to mV
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
            lfpwidth = uns['LFPPRWIDTH'] if pr else uns['LFPSPECGRAMWIDTH'] # sec
        if lfptres == None:
            lfptres = uns['LFPPRTRES'] if pr else uns['LFPSPECGRAMTRES'] # sec
        if loband == None:
            loband = uns['LFPPRLOBAND']
        f0, f1 = loband
        if hiband == None:
            hiband = uns['LFPPRHIBAND']
        f2, f3 = hiband

        assert lfptres <= lfpwidth
        NFFT = intround(lfpwidth * self.sampfreq)
        noverlap = intround(NFFT - lfptres * self.sampfreq)
        #print('len(x), NFFT, noverlap: %d, %d, %d' % (len(x), NFFT, noverlap))
        # t is midpoints of timebins in sec from start of data. P is in mV^2?:
        P, freqs, Pt = mpl.mlab.specgram(x, NFFT=NFFT, Fs=self.sampfreq, noverlap=noverlap)
        # don't convert power to dB, just washes out the signal in the ratio:
        #P = 10. * np.log10(P)
        if not relative2t0:
            Pt += t0 # convert t to time from start of ADC clock:
        nfreqs = len(freqs)

        # keep only freqs between f0 and f1, and f2 and f3:
        f0i, f1i, f2i, f3i = freqs.searchsorted([f0, f1, f2, f3])
        lP = P[f0i:f1i] # nsubfreqs x nt
        hP = P[f2i:f3i] # nsubfreqs x nt
        lP = lP.sum(axis=0) # nt
        hP = hP.sum(axis=0) # nt

        if pr:
            t = Pt
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
        #self.si_plot(hP, Pt, t0=0, t1=t[-1], ylim=None, ylabel='highband power',
        #             title=lastcmd()+' highband power', text=self.r.name)

        # set some plotting defaults:
        hlines = []
        if pr:
            ylim = 0, 1
            yticks = 0, 0.2, 0.4, 0.6, 0.8, 1
        else:
            ylim = -1, 1
            yticks = -1, 0, 1
            hlines = [0]

        # calculate some metric of each column, i.e. each bin:
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
            ytiks = 0, 1, 2
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
            # calculate xlim, always start from 0, add half a bin width to xmax:
            if pr:
                xlim = (0, t[-1]+lfpwidth/2)
            else:
                xlim = (0, t[-1]+width/2)
            self.si_plot(si, t, t0=t0, t1=t1, xlim=xlim, ylim=ylim, yticks=yticks,
                         ylabel=ylabel, showxlabel=showxlabel, showylabel=showylabel,
                         showtitle=showtitle, title=title,
                         reclabel=reclabel, hlines=hlines,
                         showstates=showstates, statelinepos=statelinepos, lw=lw,
                         alpha=alpha, swapaxes=swapaxes, figsize=figsize)
        #np.seterr(**old_settings) # restore old settings
        return si, t # t are midpoints of bins, offset depends on relative2t0
    '''
    def si_hilbert(self, chani=-1, loband=None, hiband=None, ratio='L/(L+H)',
                   plot=True):
        """Return synchrony index, i.e. power ratio of low vs high bands, as measured by
        Hilbert transform (Saleem2010). Use either L/(L+H) ratio (Saleem2010) or L/H ratio
        (Li, Poo, Dan 2009)"""
        if loband == None:
            loband = uns['LFPPRLOBAND']
        f0, f1 = loband
        if hiband == None:
            hiband = uns['LFPPRHIBAND']
        f2, f3 = hiband
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
            self.si_plot(r, t, t0, ylabel, title=lastcmd(), text=self.r.name)
        return r, t
    '''
    def si_plot(self, si, t, t0=None, t1=None, xlim=None, ylim=None, yticks=None,
                ylabel=None, showxlabel=True, showylabel=True, showtitle=True,
                title=None, reclabel=True, hlines=[0], showstates=False,
                statelinepos=None, lw=4, alpha=1,
                swapaxes=False, figsize=None):
        """Plot synchrony index as a function of time, with hopefully the same
        temporal scale as some of the other plots in self"""
        uns = get_ipython().user_ns
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)

        xlabel = "time (s)"
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

        # plot horizontal bars over time demarcating different ranges of SI values,
        # or manually defined desynched and synched periods:
        if showstates in [True, 'auto']:
            stranges, states = self.si_split(si, t)
            if swapaxes:
                lines = a.vlines
            else:
                lines = a.hlines
            slposs = statelinepos
            if len(slposs) == 1: # use same statelinepos for all states
                nstranges = len(stranges)
                slposs = slposs * nstranges
            for strange, state, slpos in zip(stranges, states, slposs):
                clr = uns['LFPPRBINCOLOURS'][state]
                lines(slpos, strange[0], strange[1], colors=clr, lw=lw, alpha=alpha)
        elif showstates == 'manual':
            REC2STATETRANGES = uns['REC2STATETRANGES']
            dtrange, strange = np.asarray(REC2STATETRANGES[self.r.absname]) / 1e6
            dtrange = max(dtrange[0], t0), min(dtrange[1], t1) # clip desynch trange to t0, t1
            strange = max(strange[0], t0), min(strange[1], t1) # clip synch trange to t0, t1
            if swapaxes:
                lines = a.vlines
            else:
                lines = a.hlines
            slposs = statelinepos
            if len(slposs) == 1: # use same statelinepos for both states
                slposs = slposs * 2
            lines(slposs[0], dtrange[0], dtrange[1], colors='b', lw=lw, alpha=alpha)
            lines(slposs[1], strange[0], strange[1], colors='r', lw=lw, alpha=alpha)

        a.plot(t, si, 'k-')
        # depending on relative2t0 in si(), x=0 represents either t0 or time ADC clock started:
        a.set_xlim(xlim) # low/high limits are unchanged if None
        a.set_ylim(ylim)
        if yticks != None:
            a.set_yticks(yticks)
        if showxlabel:
            a.set_xlabel(xlabel)
        if showylabel:
            a.set_ylabel(ylabel)
        #a.autoscale(axis='x', enable=True, tight=True)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        if title == None:
            title = lastcmd()
        gcfm().window.setWindowTitle(title)
        if showtitle:
            a.set_title(title)
        if reclabel:
            a.text(0.994, 0.01, '%s' % self.r.absname, color='k', transform=a.transAxes,
                   horizontalalignment='right', verticalalignment='bottom')
        f.tight_layout(pad=0.3) # crop figure to contents

    def si_split(self, si, t):
        """Return contiguous time ranges of SI values that fall between bin
        edges in LFPPRBINLEDGES. Return an array of tranges and the associated states"""
        nt = len(t)
        assert len(si) == nt
        bw = np.diff(t).mean() # time bin width
        # t from self.si() represents mid time points for each value of SI,
        # calculate left and right time points around each SI value instead:
        tedges = np.asarray(list(t-bw/2) + [t[-1]+bw/2])
        stateis = []
        stateis = np.zeros(nt, dtype=int)
        uns = get_ipython().user_ns
        ledges = uns['LFPPRBINLEDGES']
        redges = ledges[1:] + [np.inf]
        for statei, (ledge, redge), in enumerate(zip(ledges, redges)):
            tis = (ledge <= si) & (si < redge)
            stateis[tis] = statei
        tranis = list(np.diff(stateis).nonzero()[0] + 1) # state transition indices
        ltis = [0] + tranis # left tedges indices
        rtis = tranis + [-1] # right tedges indices
        tranges, states = [], []
        for lti, rti in zip(ltis, rtis):
            trange = tedges[lti], tedges[rti]
            state = stateis[lti]
            tranges.append(trange)
            states.append(state)
        tranges = np.asarray(tranges)
        states = np.asarray(states)
        return tranges, states

    def filterwavelet(self, chanis=None, wname="db4", maxlevel=6):
        """Filter data using wavelet multi-level decomposition and reconstruction (WMLDR).
        See Wiltschko2008"""
        data = self.get_data()
        if chanis == None:
            chanis = np.arange(len(data))
        data = data[chanis]
        data = filter.wavelet(data, wname, maxlevel)
        self.data[chanis] = data
