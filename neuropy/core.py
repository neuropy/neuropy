"""Miscellaneous functions and classes"""

from __future__ import division
from __future__ import print_function

import os
import sys
import time
import types
import __main__
import struct
import re
import random
import math
import datetime

from copy import copy
from pprint import pprint
printraw = sys.stdout.write # useful for raw printing

from PyQt4 import QtGui
from PyQt4.QtGui import QPixmap, QImage, QPalette, QColor
from PyQt4.QtCore import Qt, QSize

import numpy as np
# make overflow, underflow, div by zero, and invalid all raise errors
# this really should be the default in numpy...
np.seterr(all='raise', under='ignore') # raise all except float underflow
import scipy.signal
from scipy.special import cbrt # real cube root
import scipy.stats
from scipy.spatial.distance import pdist#, squareform
from scipy.stats import linregress

import matplotlib as mpl
import matplotlib.cm
import matplotlib.pyplot as plt
import pylab as pl
from pylab import get_current_fig_manager as gcfm
from matplotlib.collections import LineCollection

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

import filter
from colour import CCWHITERGBDICT1

TAB = '    ' # 4 spaces
EPOCH = datetime.datetime(1899, 12, 30, 0, 0, 0) # epoch for datetime stamps in .ptcs


class dictattr(dict):
    """Dictionary with attribute access. Copied from dimstim.Core"""
    def __init__(self, *args, **kwargs):
        super(dictattr, self).__init__(*args, **kwargs)
        for k, v in kwargs.iteritems():
            # call our own __setitem__ so we get keys as attribs even on kwarg init:
            self.__setitem__(k, v)
    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError, '%r object has no attribute %r' % ('dictattr', key)
    
    def __setattr__(self, key, val):
        self[key] = val

    def __getitem__(self, key):
        """On KeyError, see if converting the key from an int to a 1 or 2 digit str
        works instead"""        
        try:
            return super(dictattr, self).__getitem__(key)
        except KeyError, e: # try converting key to str of up to 2 digits in length
            for ndigits in [1, 2]:
                try:
                    #print('key: %r' % key)
                    #print('ndigits: %d' % ndigits)
                    key = pad0s(key, ndigits=ndigits)
                    #print('padded key: %r' % key)
                except ValueError:
                    raise e
                try:
                    return super(dictattr, self).__getitem__(key)
                except KeyError:
                    pass
            raise e
                
    def __setitem__(self, key, val):
        super(dictattr, self).__setitem__(key, val)
        # key isn't a number or a string starting with a number:
        if key.__class__ == str and not key[0].isdigit():
            key = key.replace(' ', '_') # get rid of any spaces
            self.__dict__[key] = val # make the key show up as an attrib upon dir()


class PTCSHeader(object):
    """Polytrode clustered spikes file header"""
    def __init__(self):
        # call the appropriate method:
        self.VER2FUNC = {1: self.read_ver_1, 2: self.read_ver_2, 3: self.read_ver_3}

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def read(self, f):
        """Read in format version, followed by rest according to verison

        formatversion: int64 (currently version 1)
        """
        self.FORMATVERSION = int(np.fromfile(f, dtype=np.int64, count=1)) # formatversion
        self.VER2FUNC[self.FORMATVERSION](f) # call the appropriate method
        
    def read_ver_1(self, f):
        """Read in header of .ptcs file version 1. For text fields, rstrip both null
        and space bytes, since NVS generated .ptcs files mix the two for padding.

        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        descr: ndescrbytes of ASCII text
            (padded with null bytes if needed for 8 byte alignment)

        nneurons: uint64 (number of neurons)
        nspikes: uint64 (total number of spikes)
        nsamplebytes: uint64 (number of bytes per template waveform sample)
        samplerate: uint64 (Hz)

        npttypebytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        pttype: npttypebytes of ASCII text
            (padded with null bytes if needed for 8 byte alignment)
        nptchans: uint64 (total num chans in polytrode)
        chanpos: nptchans * 2 * float64
            (array of (x, y) positions, in um, relative to top of polytrode,
             indexed by 0-based channel IDs)
        nsrcfnamebytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        srcfname: nsrcfnamebytes of ASCII text
            (source file name, probably .srf, padded with null bytes if needed for
             8 byte alignment)
        datetime: float64
            (absolute datetime corresponding to t=0 us timestamp, stored as days since
             epoch: December 30, 1899 at 00:00)
        ndatetimestrbytes: uint64 
        datetimestr: ndatetimestrbytes of ASCII text
            (human readable string representation of datetime, preferrably ISO 8601,
             padded with null bytes if needed for 8 byte alignment)
        """
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr
        try:
            self.descr = eval(self.descr) # should come out as a dict
        except: pass
        
        self.nneurons = int(np.fromfile(f, dtype=np.uint64, count=1)) # nneurons
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        self.nsamplebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsamplebytes
        self.samplerate = int(np.fromfile(f, dtype=np.uint64, count=1)) # samplerate

        self.npttypebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # npttypebytes
        self.pttype = f.read(self.npttypebytes).rstrip('\0 ') # pttype
        self.nptchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nptchans
        self.chanpos = np.fromfile(f, dtype=np.float64, count=self.nptchans*2) # chanpos
        self.chanpos.shape = self.nptchans, 2 # reshape into rows of (x, y) coords
        self.nsrcfnamebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsrcfnamebytes
        self.srcfname = f.read(self.nsrcfnamebytes).rstrip('\0 ') # srcfname
        # maybe convert this to a proper Python datetime object in the Neuron:
        self.datetime = float(np.fromfile(f, dtype=np.float64, count=1)) # datetime (days)
        self.ndatetimestrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndatetimestrbytes
        self.datetimestr = f.read(self.ndatetimestrbytes).rstrip('\0 ') # datetimestr

    def read_ver_2(self, f):
        """Same header as version 1. NVS created some version 1 files incorrectly, and
        incremented to version 2 for the correctly exported ones"""
        return self.read_ver_1(f)

    def read_ver_3(self, f):
        """Same header as version 2, only difference is in PTCSNeuronRecord"""
        return self.read_ver_1(f)


class PTCSNeuronRecord(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, header):
        # call the appropriate method:
        self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2, 3:self.read_ver_3}
        self.header = header
        nsamplebytes = self.header.nsamplebytes
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[nsamplebytes]

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def read(self, f):
        self.VER2FUNC[self.header.FORMATVERSION](f) # call the appropriate method
        
    def read_ver_1(self, f):
        """Read in neuron record of .ptcs file version 1

        nid: int64 (signed neuron id, could be -ve, could be non-contiguous with previous)
        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment, defaults to 0)
        descr: ndescrbytes of ASCII text
            (padded with null bytes if needed for 8 byte alignment)
        clusterscore: float64
        xpos: float64 (um)
        ypos: float64 (um)
        zpos: float64 (um) (defaults to NaN)
        nchans: uint64 (num chans in template waveforms)
        chanids: nchans * uint64 (0 based IDs of channels in template waveforms)
        maxchanid: uint64 (0 based ID of max channel in template waveforms)
        nt: uint64 (num timepoints per template waveform channel)
        nwavedatabytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavedata: nwavedatabytes of nsamplebytes sized floats
            (template waveform data, laid out as nchans * nt, in uV,
             padded with null bytes if needed for 8 byte alignment)
        nwavestdbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavestd: nwavestdbytes of nsamplebytes sized floats
            (template waveform standard deviation, laid out as nchans * nt, in uV,
             padded with null bytes if needed for 8 byte alignment)
        nspikes: uint64 (number of spikes in this neuron)
        spike timestamps: nspikes * uint64 (us, should be sorted)
        """
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr
        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass
        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xpos = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.ypos = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zpos = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um)
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) # chanids
        self.maxchan = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid
        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt
        self.nwavedatabytes, self.wavedata = self.read_wave(f)
        self.nwavestdbytes, self.wavestd = self.read_wave(f)
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        # spike timestamps (us):
        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes)
        # convert from unsigned to signed int for calculating intervals:
        self.spikes = np.asarray(self.spikes, dtype=np.int64)

    def read_wave(self, f):
        """Read wavedata/wavestd bytes"""
        # nwavedata/nwavestd bytes, padded:
        nbytes = int(np.fromfile(f, dtype=np.uint64, count=1))
        fp = f.tell()
        count = nbytes // self.header.nsamplebytes # trunc to ignore any pad bytes
        X = np.fromfile(f, dtype=self.wavedtype, count=count) # wavedata/wavestd (uV)
        if nbytes != 0:
            X.shape = self.nchans, self.nt # reshape
        f.seek(fp + nbytes) # skip any pad bytes
        return nbytes, X

    def read_ver_2(self, f):
        """Same as version 1. NVS created some version 1 files incorrectly, and
        incremented to version 2 for the correctly exported ones"""
        return self.read_ver_1(f)

    def read_ver_3(self, f):
        """Read in neuron record of .ptcs file version 3. 'zpos' field was replaced
        by 'sigma' field.

        nid: int64 (signed neuron id, could be -ve, could be non-contiguous with previous)
        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment, defaults to 0)
        descr: ndescrbytes of ASCII text
            (padded with null bytes if needed for 8 byte alignment)
        clusterscore: float64
        xpos: float64 (um)
        ypos: float64 (um)
        sigma: float64 (um) (Gaussian spatial sigma)
        nchans: uint64 (num chans in template waveforms)
        chanids: nchans * uint64 (0 based IDs of channels in template waveforms)
        maxchanid: uint64 (0 based ID of max channel in template waveforms)
        nt: uint64 (num timepoints per template waveform channel)
        nwavedatabytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavedata: nwavedatabytes of nsamplebytes sized floats
            (template waveform data, laid out as nchans * nt, in uV,
             padded with null bytes if needed for 8 byte alignment)
        nwavestdbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavestd: nwavestdbytes of nsamplebytes sized floats
            (template waveform standard deviation, laid out as nchans * nt, in uV,
             padded with null bytes if needed for 8 byte alignment)
        nspikes: uint64 (number of spikes in this neuron)
        spike timestamps: nspikes * uint64 (us, should be sorted)
        """
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr
        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass
        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xpos = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.ypos = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.sigma = float(np.fromfile(f, dtype=np.float64, count=1)) # sigma (um)
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) # chanids
        self.maxchan = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid
        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt
        self.nwavedatabytes, self.wavedata = self.read_wave(f)
        self.nwavestdbytes, self.wavestd = self.read_wave(f)
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        # spike timestamps (us):
        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes)
        # convert from unsigned to signed int for calculating intervals:
        self.spikes = np.asarray(self.spikes, dtype=np.int64)


class SPKHeader(object):
    """Represents a folder containing neurons in .spk files. Similar to a
    PTCSHeader, but much more impoverished"""
    def __init__(self, path):
        self.path = path
        fnames = [ fname for fname in os.listdir(self.path)
                   if os.path.isfile(os.path.join(self.path, fname))
                   and fname.endswith('.spk') ] # spike filenames
        self.spkfnames = sorted(fnames)
        self.nspikes = 0
        # for compatibility with PTCSHeader:
        self.datetime = 0
        self.pttype = None
        self.chanpos = None
        self.samplerate = 25000 # assumed

    def read(self, neuron):
        neuron.loadspk() # load the neuron
        self.nspikes += neuron.nspikes
        # look for neuron2pos.py file, which contains a dict mapping neuron id to (x, y)
        # position
        if 'neuron2pos.py' in os.listdir(self.path):
            oldpath = os.getcwd()
            os.chdir(self.path)
            from neuron2pos import neuron2pos
            neuron.record.xpos, neuron.record.ypos = neuron2pos[neuron.id]
            os.chdir(oldpath)
            

class SPKNeuronRecord(object):
    """Represents the spike times in a simple .spk file as a record. Similar to a
    PTCSNeuronRecord, but much more impoverished"""
    def __init__(self, fname):
        self.fname = fname
        # for compatibility with PTCSHeader:
        self.descr = None
        self.xpos, self.ypos = None, None
        self.sigma = None
        self.nchans = None
        self.chans = None
        self.maxchan = None
        self.nt = None
        self.wavedata = None
        self.wavestd = None
        
    def parse_id(self):
        """Return everything from just after the last '_t' to the end of the
        fname, should be all numeric"""
        name = os.path.split(self.fname)[-1] # pathless
        name = os.path.splitext(name)[0] # extensionless
        return int(name.rsplit('_t', 1)[-1]) # id

    def read(self):
        self.nid = self.parse_id()
        with open(self.fname, 'rb') as f:
            self.spikes = np.fromfile(f, dtype=np.int64) # spike timestamps (us)
        self.nspikes = len(self.spikes)
    

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
        self.PLOTGAIN = 2

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
        # grab desired channels and time range, and AD values to uV:
        data = self.data[chanis][:, t0i:t1i] * self.uVperAD * self.PLOTGAIN * gain
        nt = len(ts)
        assert nt == data.shape[1]
        x = np.tile(ts, nchans)
        x.shape = nchans, nt
        segments = np.zeros((nchans, nt, 2)) # x vals in col 0, yvals in col 1
        segments[:, :, 0] = x
        segments[:, :, 1] = -data # set to -ve here because of invert_yaxis() below
        # add offsets:
        for chanii, chani in enumerate(chanis):
            chan = self.chans[chani]
            ypos = self.chanpos[chan][1]
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
           figsize=(20, 6.5), swapaxes=False):
        """Calculate an LFP synchrony index, using potentially overlapping windows of
        width and tres, in sec, from the LFP spectrogram, itself composed of bins of
        lfpwidth and lfptres. Options for kind are:

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
            tranges = split_tranges([trange], width, tres) # in us
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


class DensePopulationRaster(object):
    """Population spike raster plot, with dense vertical spacing according to neuron depth
    rank, and colour proportional to neuron depth"""
    def __init__(self, trange=None, neurons=None, norder=None, units='sec', r=None,
                 marker='|', size=None, color=None, alpha=1.0, title=True, figsize=(20, None)):
        """neurons is a dict, trange is time range in us to raster plot over. Raster plot
        is displayed in time units of units"""
        assert len(trange) == 2
        trange = np.asarray(trange)
        if norder != None:
            nids = norder
        else: # sort neurons by their depth rank:
            nids = np.sort(neurons.keys())
            # depth from top of electrode:
            unsorted_ypos = np.array([ neurons[nid].pos[1] for nid in nids ])
            nids = nids[unsorted_ypos.argsort()]
        self.nids = nids
        print(nids)
        # depth of nids from top of electrode
        ypos = np.array([ neurons[nid].pos[1] for nid in nids ])
        supis, midis, deepis = laminarity(ypos, r.tr.absname)
        nn = len(nids)
        t, y, c = [], [], []
        for nidi, nid in enumerate(nids):
            n = neurons[nid]
            lo, hi = n.spikes.searchsorted(trange)
            spikes = n.spikes[lo:hi]
            nspikes = len(spikes)
            if nspikes > 0:
                t.append(spikes)
                y.append(np.tile(nidi, nspikes)) # depth rank below top of electrode
                if supis[nidi]: nc = 'r'
                elif midis[nidi]: nc = 'g'
                elif deepis[nidi]: nc = 'b'
                else: nc = 'y'
                c.append(np.tile(nc, nspikes))

        t = np.hstack(t)
        # spike time multiplier to use for raster labels:
        tx = {'us': 1, 'ms': 1000, 'sec': 1000000}[units]
        if tx != 1:
            t = t / tx # don't do in-place, allow conversion to float
        y = np.hstack(y)
        c = np.hstack(c)

        s = 50
        # manual size or color settings override the automatic values
        if size != None:
            s = size
        if color != None:
            c = color

        if figsize[1] == None:
            figsize = figsize[0], 1 + nn / 7 # ~1/7th vertical inch per neuron
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.scatter(t, y, marker=marker, c=c, alpha=alpha, s=s)
        a.set_xlim(trange/tx)
        a.set_ylim(nn, -1) # invert the y axis
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        a.set_xlabel("time (%s)" % units)
        if norder == None:
            a.set_ylabel("neuron depth rank")
        else:
            a.set_ylabel("neuron order")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            titlestr += ' (%s)' % r.name
            a.set_title(titlestr)
            # add pseudo legend of coloured text:
            tmax = a.get_xlim()[1]
            #rainbow_text(a, 0.905*tmax, -1.5, ['superficial', 'middle', 'deep'],
            #             ['r', 'g', 'b'])
            a.text(0.908*tmax, -1.5, 'superficial', color='r')
            a.text(0.952*tmax, -1.5, 'middle', color='g')
            a.text(0.980*tmax, -1.5, 'deep', color='b')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f


class SpatialPopulationRaster(object):
    """Population spike raster plot, with vertical spacing proportional to neuron depth,
    colour representing neuron id"""
    def __init__(self, trange=None, neurons=None, norder=None, units='sec', r=None,
                 marker='|', size=None, color=None, alpha=1.0, title=True, figsize=(20, 6.5)):
        """neurons is a dict, trange is time range in us to raster plot over. Raster plot
        is displayed in time units of units"""
        assert len(trange) == 2
        trange = np.asarray(trange)
        if norder == None:
            nids = sorted(neurons.keys())
        else:
            nids = norder
            self.norder = norder
            print(norder)
        t, y, c, s = [], [], [], []
        for nidi, nid in enumerate(nids):
            n = neurons[nid]
            lo, hi = n.spikes.searchsorted(trange)
            spikes = n.spikes[lo:hi]
            nspikes = len(spikes)
            if nspikes > 0:
                t.append(spikes)
                if norder == None:
                    ypos = n.pos[1]
                else:
                    ypos = nidi
                y.append(np.tile(ypos, nspikes)) # -ve, distance below top of electrode
                nc = CCWHITERGBDICT1[nid] # for max colour alternation, use DICT0 and nidi
                c.append(np.tile(nc, nspikes))
                # use big points for low rate cells, small points for high rate cells:
                #ms = max(min(10000/nspikes, 50), 5)
                #s.append(np.tile(ms, nspikes))
        t = np.hstack(t)
        # spike time multiplier to use for raster labels:
        tx = {'us': 1, 'ms': 1000, 'sec': 1000000}[units]
        if tx != 1:
            t = t / tx # don't do in-place, allow conversion to float
        y = np.hstack(y)
        c = np.hstack(c)
        c.shape = -1, 3
        #s = np.hstack(s)
        s = 50

        # manual size or color settings override the automatic values
        if size != None:
            s = size
        if color != None:
            c = color

        if figsize[1] == None:
            figsize = figsize[0], 6.5
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        ec = 'none'
        if marker == '|':
            ec = c
        a.scatter(t, y, marker=marker, c=c, edgecolor=ec, alpha=alpha, s=s)
        a.invert_yaxis() # increasingly +ve values down y axis
        a.set_xlim(trange/tx)
        if norder == None: # set y axis limits according to spatial extent of probe
            # grab first neuron's sort.chanpos, should be the same for all:
            chanpos = neurons[nids[0]].sort.chanpos
            ymax = chanpos[:, 1].max() # max chan distance below top of probe
            ymax = np.ceil(ymax / 50) * 50 # round up to nearest multiple of 100 um
            a.set_ylim(ymax, 0)
        else: # autoscale 'nidi' integer y axis
            a.autoscale(enable=True, axis='y', tight=True)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        a.set_xlabel("time (%s)" % units)
        if norder == None:
            a.set_ylabel("depth ($\mu$m)")
        else:
            a.set_ylabel("neuron order")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            titlestr += ' (%s)' % r.name
            a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
    '''
    def _onmotion(self, event):
        """Called during mouse motion over figure. Pops up neuron and
        experiment info in a tooltip when hovering over a neuron row."""
        self.f.canvas.mpl_connect('motion_notify_event', self._onmotion)
        self.f.canvas.mpl_connect('key_press_event', self._onkeypress)

        if event.inaxes: # if mouse is inside the axes
            # use ydata to get index into sorted list of neurons:
            nii = int(math.floor(event.ydata))
            ni = self.nids[nii]
            neuron = self.neurons[ni]
            currentexp = None
            for e in self.e.values(): # for all experiments
                estart = (e.trange[0]-self.t0)/self.tconv
                eend = (e.trange[1]-self.t0)/self.tconv
                if estart < event.xdata  < eend:
                    currentexp = e
                    break # don't need to check any of the other experiments
            # print timepoint down to nearest us, in units of ms:
            tip = 't: %.3f ms\n' % event.xdata
            tip += 'n%d: %d spikes' % (neuron.id, neuron.nspikes)
            if currentexp == None:
                tip += '\nno experiment'
            else:
                tip += '\nexperiment %s: %r' % (currentexp.id, currentexp.name)
            self.tooltip.SetTip(tip) # update the tooltip
            self.tooltip.Enable(True) # make sure it's enabled
        else: # mouse is outside the axes
            self.tooltip.Enable(False) # disable the tooltip

    def _onkeypress(self, event):
        """Called during a figure keypress"""
        key = event.guiEvent.GetKeyCode() # wx dependent
        #print(key)
        # could also just use the backend-neutral event.key, but that doesn't recognize
        # as many keypresses, like pgup, pgdn, etc.
        if not event.guiEvent.ControlDown(): # Ctrl key isn't down, wx dependent
            if key == wx.WXK_RIGHT: # pan right
                self._panx(+0.1)
            elif key == wx.WXK_LEFT: # pan left
                self._panx(-0.1)
            elif key == wx.WXK_UP: # zoom in
                self._zoomx(1.2)
            elif key == wx.WXK_DOWN: # zoom out
                self._zoomx(1/1.2)
            elif key == wx.WXK_NEXT: # PGDN (page right)
                self._panx(+1)
            elif key == wx.WXK_PRIOR: # PGUP (page left)
                self._panx(-1)
            elif key == wx.WXK_HOME: # go to start of first Experiment
                self._panx(left=self.experimentmarkers[0])
            elif key == wx.WXK_END: # go to end of last Experiment
                self._panx(left=self.experimentmarkers[-1]-self.width)
            elif key == ord('['): # skip backwards to previous jump point
                # current position of left edge of the window in jumpts list:
                i = self.jumpts.searchsorted(self.left, side='left')
                i = max(0, i-1) # decrement by 1, do bounds checking
                self._panx(left=self.jumpts[i])
            elif key == ord(']'): # skip forwards to next jump point
                # current position of left edge of the window in jumpts list:
                i = self.jumpts.searchsorted(self.left, side='right')
                i = min(i, len(self.jumpts)-1) # bounds checking
                self._panx(left=self.jumpts[i])
            elif key == wx.WXK_RETURN: # go to position
                self._goto()
            elif key == ord(','): # cycle tick formatter through thousands separators
                self._cyclethousandssep()
            elif key == ord('B'): # toggle plotting of bin edges
                self._togglebinedges()
        else: # Ctrl key is down
            if key == wx.WXK_LEFT: # skip backwards to previous experiment marker
                # current position of left edge of the window in experimentmarkers list:
                i = self.experimentmarkers.searchsorted(self.left, side='left')
                i = max(0, i-1) # decrement by 1, do bounds checking
                self._panx(left=self.experimentmarkers[i])
            elif key == wx.WXK_RIGHT: # skip forwards to next experiment marker
                # current position of left edge of the window in experimentmarkers list:
                i = self.experimentmarkers.searchsorted(self.left, side='right')
                i = min(i, len(self.experimentmarkers)-1) # bounds checking
                self._panx(left=self.experimentmarkers[i])
            elif key == wx.WXK_UP: # zoom in faster
                self._zoomx(3.0)
            elif key == wx.WXK_DOWN: # zoom out faster
                self._zoomx(1/3.0)
    '''

class Codes(object):
    """A 2D array where each row is a neuron code, and each column
    is a binary population word for that time bin, sorted LSB to MSB from top to bottom.
    neurons is a list of Neurons, also from LSB to MSB. Order in neurons is preserved."""
    def __init__(self, neurons=None, tranges=None, shufflecodes=False):
        self.neurons = neurons
        self.tranges = tolist(tranges)
        self.shufflecodes = shufflecodes
        self.nids = [ neuron.id for neuron in self.neurons ]
        self.nneurons = len(self.neurons)
        # make a dict from keys:self.nids, vals:range(self.nneurons). This converts from nids
        # to niis (from neuron indices to indices into the binary code array self.c)
        self.nids2niisdict = dict(zip(self.nids, range(self.nneurons)))

    def nids2niis(self, nids=None):
        """Converts from nids to niis (from neuron indices to indices into the binary code
        array self.c). nids can be a sequence"""
        try:
            return [ self.nids2niisdict[ni] for ni in nids ]
        except TypeError: # iteration over non-sequence, nids is a scalar
            return self.nids2niisdict[nids]

    def calc(self):
        self.c = [] # stores the 2D code array
        # append neurons in their order in self.neurons, store them LSB to MSB from top to
        # bottom
        for neuron in self.neurons:
            codeo = neuron.code(tranges=self.tranges)
            # build up nested list (ie, 2D) of spike times, each row will have different
            # length:
            if self.shufflecodes:
                c = codeo.c.copy() # make a copy (leave the codeo's codetrain untouched)
                np.random.shuffle(c) # shuffle each neuron's codetrain separately, in-place
            else:
                c = codeo.c # just a pointer
            self.c.append(c) # flat list
        # store the bin edges, for reference. All bin times should be the same for all
        # neurons, because they're all given the same trange. use the bin times of the last
        # neuron
        self.t = codeo.t
        nneurons = len(self.neurons)
        nbins = len(self.c[0]) # all entries in the list should be the same length
        self.c = np.vstack(self.c)

    def syncis(self):
        """Returns synch indices, ie the indices of the bins for which all the
        neurons in this Codes object have a 1 in them"""
        # take product down all rows, only synchronous events across all cells will survive:
        return self.c.prod(axis=0).nonzero()[0]

    def syncts(self):
        """Returns synch times, ie times of the left bin edges for which
        all the neurons in this Codes object have a 1 in them"""
        return self.t[self.syncis()]

    def synctsms(self):
        """Returns synch times in ms, to the nearest ms"""
        return np.int32(np.round(self.syncts() / 1e3))

    def copy(self):
        """Returns a copy of the Codes object"""
        return copy(self)
    '''
    # needs some testing:
    def append(self, others):
        """Adds other Codes objects appended in time (horizontally) to this Codes object.
        Useful for appending Codes objects across Recordings ? (don't really need it
        for appending across Experiments)"""
        others = tolist(others)
        for other in others:
            assert other.neurons == self.neurons
        codesos = [self] # list of codes objects
        codesos.extend(others)
        # this tranges potentially holds multiple tranges from each codes objects,
        # times the number of codes objects
        self.tranges = [ trange for codeso in codesos for trange in codeso.tranges ]
        self.calc() # recalculate this code with its new set of tranges
    '''
    
class SpikeCorr(object):
    """Calculate and plot spike correlations of all cell pairs from nids (or of all
    cell pairs within some torus of radii R=(R0, R1) in um) in source, during tranges
    or experiments. If width is not None, calculate self as a function of time, with bin
    widths width sec and time resolution tres sec. For each pair, shift the
    second spike train by shift ms, or shift it by shiftcorrect ms and subtract the
    correlation from the unshifted value."""
    def __init__(self, source, tranges=None, width=None, tres=None,
                 shift=0, shiftcorrect=0, nidskind=None, R=None):
        uns = get_ipython().user_ns
        recs, tracks = parse_source(source)
        nidss = get_nids(recs, tracks, kind=nidskind)
        self.recs, self.tracks, self.nidss = recs, tracks, nidss

        # set some kind of representative name:
        if len(recs) == 1:
            self.name = recs[0].name
        else:
            self.name = ', '.join([ track.absname for track in tracks ])
        
        if tranges != None and len(recs) > 1: # more than one recording in source:
            raise ValueError("can't analyze across recordings when tranges selected")
        
        # calculate Ising code matrix (or matrices) once:
        if len(tracks) > 1: # recordings from multiple tracks, collect lists of code arrs
            codes = []
            for nids, track in zip(nidss, tracks):
                trcodes = [ rec.codes(nids=nids, tranges=tranges).c for rec in recs
                            if rec.tr == track ]
                trcodes = np.hstack(trcodes)
                codes.append(trcodes)
        else: # all recordings from same track
            reccodes = [ rec.codes(nids=nidss[0], tranges=tranges) for rec in recs ]
            codes = [ np.hstack(reccode.c for reccode in reccodes) ]
        self.codes = codes # list of Ising code matrices

        if len(recs) == 1: # grab actual time values for single recording
            ts = [ reccodes[0].t ]
            tranges = [ reccodes[0].tranges ]
        else: # generate fake time values, starting from 0, one set per code array
            ts, tranges = [], []
            bw = uns['CODETRES'] # us
            for code in codes: # iterate over list of code arrays
                nbins = code.shape[1]
                t = np.arange(0, nbins*bw, bw)
                trange = t[0], t[-1]
                ts.append(t)
                tranges.append(trange)
        self.ts = ts # list of timestamp arrs corresponding to bin edges in codes
        self.tranges = tranges # list of tranges

        # width and tres are for calculating a metric of all corrs as f'n of time
        if width != None:
            if tres == None:
                tres = width
            assert tres <= width
            width = intround(width * 1000000) # convert from sec to us
            tres = intround(tres * 1000000) # convert from sec to us
        self.width = width
        self.tres = tres

        if shift or shiftcorrect:
            raise NotImplementedError("shift and shiftcorrect are currently disabled to "
                                      "simplify the logic")
            self.shift = shift # shift spike train of the second of each neuron pair, in ms
            # shift correct spike train of the second of each neuron pair by this much, in ms
            self.shiftcorrect = shiftcorrect
        if R != None:
            raise NotImplementedError("torus code hasn't been tested in a long time")
            assert len(R) == 2 and R[0] < R[1]  # should be R = (R0, R1) torus
            self.R = R

    def calc(self):
        uns = get_ipython().user_ns
        corrs, counts, pairs = [], [], []
        if self.width != None:
            tranges = []
            highval = uns['CODEVALS'][1]
            for code, t, trange in zip(self.codes, self.ts, self.tranges):
                # compute correlation coefficients as a function of time, one value per trange:
                trange = split_tranges(trange, self.width, self.tres)
                corr, count = util.sct(code, t, trange, highval)
                nneurons = len(code)
                corrs.append(corr)
                counts.append(counts)
                pairs.append(np.asarray(np.triu_indices(nneurons, k=1)).T)
                tranges.append(trange)
            self.tranges = tranges # overwrite
        else:
            corrs, counts, pairs = [], [], []
            for code, nids in zip(self.codes, self.nidss):
                # compute correlation coefficients once across entire set of tranges:
                corr, count, pair = self.calc_single(code, nids)
                corrs.append(corr)
                counts.append(count)
                pairs.append(pair)
        self.corrs = corrs
        self.counts = counts
        self.pairs = pairs
        self.npairs = [ len(pair) for pair in pairs ]

    def calc_single(self, code, nids):
        """Calculate one spike correlation value for each cell pair, given one code array
        spanning some subset of self.tranges"""
        nneurons, nbins = code.shape
        print('nneurons, nbins = %d, %d' % (nneurons, nbins))
        code = np.float64(code) # prevent int8 overflow somewhere

        # precalculate mean and std of each cell's codetrain, rows correspond to nids:
        means = code.mean(axis=1)
        stds = code.std(axis=1)

        # precalculate number of high states in each neuron's code:
        nhigh = np.zeros(nneurons, dtype=np.int64)
        uns = get_ipython().user_ns
        if uns['CODEVALS'] != [0, 1]:
            raise RuntimeError("counting of high states assumes CODEVALS = [0, 1]")
        for nii0 in range(nneurons):
            nhigh[nii0] = code[nii0].sum()
        
        alpha = uns['ALPHA']
        nrejected = 0

        #shift, shiftcorrect = self.shift, self.shiftcorrect
        #if shift and shiftcorrect:
        #    raise ValueError("only one of shift or shiftcorrect can be nonzero")

        # iterate over all pairs:
        #n = self.r.n
        #R = self.R
        corrs = []
        counts = []
        pairs = []
        for nii0 in range(nneurons):
            ni0 = nids[nii0]
            for nii1 in range(nii0+1, nneurons):
                ni1 = nids[nii1]
                # skip the pair if a torus is specified and if
                # the pair's separation falls outside bounds of specified torus:
                #if R != None and not (R[0] < dist(n[ni0].pos, n[ni1].pos) < R[1]):
                #    continue # to next pair
                # potentially shift only the second code train of each pair:
                #c0 = self.r.n[ni0].code(tranges=tranges).c
                #c1 = self.r.n[ni1].code(tranges=tranges, shift=shift).c
                c0 = code[nii0]
                c1 = code[nii1]
                # (mean of product - product of means) / product of stds:
                numer = np.dot(c0, c1) / nbins - means[nii0] * means[nii1]
                denom = stds[nii0] * stds[nii1]
                if numer == 0.0:
                    sc = 0.0 # even if denom is also 0
                elif denom == 0.0: # numer is not 0, but denom is 0, prevent div by 0
                    print('skipped pair (%d, %d) in r%s' % (ni0, ni1, self.r.id))
                    continue # skip to next pair
                else:
                    sc = numer / denom
                # potentially shift correct using only the second spike train of each pair:
                #if shiftcorrect:
                #    c1sc = self.r.n[ni1].code(tranges=tranges, shift=shiftcorrect).c
                #    scsc = ((c0 * c1sc).mean() - means[ni0] * means[ni1]) / denom
                #    ## TODO: might also want to try subtracting abs(scsc)?
                #    sc -= scsc
                # calculate t value for pearson correlation, see
                # http://www.vassarstats.net/textbook/ch4apx.html:
                t = sc / np.sqrt((1 - sc**2)/(nbins - 2))
                # calculate corresponding two-sided pval = Prob(abs(t)>tt), see
                # http://docs.scipy.org/doc/scipy/reference/tutorial/stats.html:
                p = 2*scipy.stats.t.sf(np.abs(t), nbins-1)
                # tested and compared to scipy.pearsonr, identical results, but above should
                # be faster because it avoids unnecessarily recomputing means and stdevs:
                #sc2, p2 = scipy.stats.pearsonr(c0, c1)
                #print(sc, sc2, p, p2)
                if p >= alpha:
                    nrejected += 1
                    continue
                corrs.append(sc)
                pairs.append([nii0, nii1])
                # take sum of high code counts of pair. Note that taking the mean wouldn't
                # change results in self.sct(), because it would end up simply normalizing
                # by half the value
                counts.append(nhigh[nii0] + nhigh[nii1])
        print('%d of %d pairs rejected' % (nrejected, nCr(nneurons, 2)))
        corrs = np.asarray(corrs)
        counts = np.asarray(counts)
        pairs = np.asarray(pairs)
        return corrs, counts, pairs

    def clear_codes(self):
        """Delete cached codes from all recordings"""
        for rec in self.recs:
            for n in rec.alln.values():
                try:
                    del n._codes
                except AttributeError:
                    pass

    def norder(self, metric=False, n_init=10, max_iter=1000, verbose=0, eps=-np.inf,
               n_jobs=1, init=None):
        """Return nids of self's recording sorted according to their pairwise correlations.
        This uses multidimensional scaling (MDS) to take N*(N-1)/2 pairwise values and
        return a 1D projection in which cells with the greatest similarity are kept as close
        to each other as possible, with as little "stress" as possible. This might then be
        useful for sorting raster plots to better reveal ensemble activity.

        Note that so far, this generates orderings that are poorly reproducible across runs,
        and seem fairly random. There is some confusion in the docs for sklearn.manifold.MDS
        as to whether to pass it a similarity matrix (high values mean that pair should be
        kept close together) or a dissimilarity matrix (high values mean that pair should be
        kep far apart). Also, the eps and max_iter kwargs are a bit deceiving. Usually, the
        algorithm reaches a minimum stress after only the 2 or 3 iterations, and after that
        starts increasing again, ie the error becomes -ve and the algorithm exits. If
        max_iter is anything greater than about 3, it'll never be reached. If you provide a
        -ve eps however, this forces the algorithm to exit only once it reaches max_iter.
        The final stress value will be higher than the minimum, but at least it actually
        becomes stable. You can see this by setting verbose=2. Stable non-minimum stress
        values seem to provide more consistent neuron sorting across runs ("inits" in
        sklearn's language) than when simply exiting at minimum stress, but I'm not all that
        certain.

        ## TODO: try Isomap instead
        """
        from sklearn.manifold import MDS
        if len(self.nidss) > 1:
            raise ValueError('too many sources defined, which nids should I work on?')
        nids = self.nidss[0]
        self.calc()
        sim = copy(self.corrs[0]) # similarities, 1D array of length npairs
        npairs = len(sim)
        N = len(nids)
        assert npairs == N * (N-1) / 2 # sanity check
        # scale, might not be necessary, although -ve values might be bad:
        sim -= sim.min()
        sim /= sim.max()
        dissim = 1 - sim # dissimilarities
        # maybe keep all values well above 0, MDS ignores 0 values?
        dissim += 1e-5
        ui = np.triu_indices(N, 1) # upper triangle indices
        li = np.tril_indices(N, -1) # lower triangle indices
        dissimm = np.zeros((N, N)) # dissimilarity matrix, 0s on diagonal
        dissimm[ui] = dissim # fill upper triangle, which is filled row major, same as corrs
        dissimm[li] = dissimm.T[li] # make symmetric by filling lower triangle
        mds = MDS(n_components=1, metric=metric, n_init=n_init, max_iter=max_iter,
                  verbose=verbose, eps=eps, n_jobs=n_jobs, dissimilarity='precomputed')
        pos = mds.fit_transform(dissimm, init=init)
        print('lowest stress: %s' % mds.stress_)
        #print('pos:')
        #print(pos.ravel())
        sortis = pos.ravel().argsort()
        #print('sortis:')
        #print(sortis)
        norder = np.asarray(nids)[sortis]
        #print('sorted nids:')
        #print(norder)
        return norder

    def shifts(self, start=-5000, stop=5000, step=50, shiftcorrect=True, figsize=(7.5, 6.5)):
        """Plot shift-corrected, or just shifted, median pairwise spike correlations of all
        cell pairs as a function of shifts, from start to stop in steps of step ms"""
        ## TODO: update for multiple tracks
        assert step > 0
        if stop % step == 0:
            stop += step # make stop end inclusive
        assert start < stop
        shifts = np.arange(start, stop, step) # shift values, in ms
        uns = get_ipython().user_ns
        self.calc() # run it once here to init self.nids and self.pairs
        c, supis, midis, deepis, otheris = self.pair_laminarity(self.nids, self.pairs)
        nsup, nmid, ndeep, nother = len(supis), len(midis), len(deepis), len(otheris)
        allmeds = np.zeros(len(shifts)) # medians of all pairs
        supmeds = np.zeros(len(shifts)) # medians of superficial pairs
        midmeds = np.zeros(len(shifts)) # medians of middle pairs
        deepmeds = np.zeros(len(shifts)) # medians of deep pairs
        othermeds = np.zeros(len(shifts)) # medians of other pairs
        for shifti, shift in enumerate(shifts):
            # calculate corrs for each shift
            if shiftcorrect:
                self.shiftcorrect = shift
            else:
                self.shift = shift
            self.calc()
            allmeds[shifti] = np.median(self.corrs)
            # check for empty *is, which raise FloatingPointErrors in np.median:
            if nsup: supmeds[shifti] = np.median(self.corrs[supis])
            if nmid: midmeds[shifti] = np.median(self.corrs[midis])
            if ndeep: deepmeds[shifti] = np.median(self.corrs[deepis])
            if nother: othermeds[shifti] = np.median(self.corrs[otheris])
            print('%d,' % shift, end='')
        print()
        self.clear_codes() # free memory
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.plot(shifts, allmeds, 'e-o', mec='e', ms=3, label='all')
        if nsup: a.plot(shifts, supmeds, 'r-o', mec='r', ms=3, label='superficial')
        if nmid: a.plot(shifts, midmeds, 'g-o', mec='g', ms=3, label='middle')
        if ndeep: a.plot(shifts, deepmeds, 'b-o', mec='b', ms=3, label='deep')
        if nother: a.plot(shifts, othermeds, 'y-o', mec='y', ms=3, label='other')
        # underplot horizontal line at y=0:
        a.axhline(y=0, c='e', ls='--', marker=None)
        a.set_xlim(shifts[0], shifts[-1]) # override any MPL smarts
        if shiftcorrect:
            a.set_xlabel("shift correction (ms)")
            a.set_ylabel("median shift-corrected correlation coefficient")
            pos = 0.99, 0.01 # put info text in bottom right
            verticalalignment = 'bottom'
            legendloc = 'lower left'
        else:
            a.set_xlabel("shift (ms)")
            a.set_ylabel("median shifted correlation coefficient")
            verticalalignment = 'top'
            pos = 0.99, 0.99 # put info text in top right
            legendloc = 'upper left'
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        # add info text to top/bottom right of plot:
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][self.r.tr.absname]
        a.text(pos[0], pos[1], '%s\n'
                               'tres = %d ms\n'
                               'phase = %d deg\n'
                               'R = %r um\n'
                               'minrate = %.2f Hz\n'
                               'nneurons = %d\n'
                               'npairs = %d\n'
                               'sup = %r um\n'
                               'mid = %r um\n'
                               'deep = %r um\n'
                               'dt = %d min'
                               % (self.r.name, uns['CODETRES']//1000, uns['CODEPHASE'],
                                  self.R, uns['MINRATE'], len(self.nids), self.npairs,
                                  sup, mid, deep, intround(self.r.dtmin)),
                               transform=a.transAxes,
                               horizontalalignment='right',
                               verticalalignment=verticalalignment)
        # add legend:
        a.legend(loc=legendloc, markerscale=2.0, handletextpad=0.5)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def pdf(self, crange=[-0.05, 0.15], figsize=(7.5, 6.5), limitstats=False,
            nbins=30, density=True, title=True, text=True, labels=True, yticks=True):
        """Plot PDF of pairwise spike correlations. If limitstats, the stats displayed
        exclude any corr values that fall outside of crange"""
        self.calc()
        corrs = np.hstack(self.corrs) # all pairwise corrs from all tracks
        npairs = len(corrs)
        nneurons = sum([ len(nids) for nids in self.nidss ])
        dtmin = sum([ rec.dtmin for rec in self.recs ])
        nbins = max(nbins, 2*intround(np.sqrt(npairs)))

        # figure out the bin edges:
        if crange != None:
            bins = np.linspace(start=crange[0], stop=crange[1], num=nbins,
                               endpoint=True)
        else: # let np.histogram() figure out the bin edges
            bins = nbins
        n, c = np.histogram(corrs, bins=bins, density=density)
        binwidth = c[1] - c[0] # take width of first bin in c

        # potentially constrain n and c values for reporting and plotting:
        if limitstats:
            corrs = corrs[(corrs >= crange[0]) * (corrs <= crange[1])]
            n, c = np.histogram(corrs, bins=bins, density=density)
        mean = np.mean(corrs)
        median = np.median(corrs)
        argmode = n.argmax()
        mode = c[argmode] + binwidth / 2 # middle of tallest bin
        stdev = np.std(corrs)

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        # omit last right edge in c:
        a.bar(left=c[:-1], height=n, width=binwidth, bottom=0, color='k', ec='k')
        a.axvline(x=0, c='e', ls='--') # draw vertical grey line at x=0
        if crange != None:
            a.set_xlim(crange)
            a.set_xticks(np.arange(crange[0], crange[1], 0.1))
        if yticks in (False, None): # disable yticks
            a.set_yticks([])
        a.set_ylim(ymax=n.max()) # scale to height of tallest bin
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if title:
            a.set_title(titlestr)
        
        if labels:
            if density:
                a.set_ylabel('probability density')
            else:
                a.set_ylabel('count')
            a.set_xlabel('correlation coefficient')
            
        # add stuff to top right of plot:
        uns = get_ipython().user_ns
        if text:
            a.text(0.99, 0.99, '%s\n'
                               'mean = %.3f\n'
                               'median = %.3f\n'
                               'mode = %.3f\n'
                               'stdev = %.3f\n'
                               'tres = %d ms\n'
                               'phase = %d deg\n'
                               'minrate = %.2f Hz\n'
                               'nneurons = %d\n'
                               'npairs = %d\n'
                               'dt = %d min\n'
                               'alpha = %.2f'
                               % (self.name, mean, median, mode, stdev,
                                  uns['CODETRES']//1000, uns['CODEPHASE'], uns['MINRATE'],
                                  nneurons, npairs, dtmin, uns['ALPHA']),
                               transform = a.transAxes,
                               horizontalalignment='right',
                               verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def sort(self, figsize=(7.5, 6.5)):
        """Plot pairwise spike correlations in decreasing order"""
        ## TODO: update for multiple tracks
        self.calc()
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        corrs = self.corrs
        sortis = corrs.argsort()[::-1] # indices to get corrs in decreasing order
        corrs = corrs[sortis] # corrs in decreasing order
        pairs = self.pairs[sortis] # pairs in decreasing corrs order
        npairs = len(pairs)

        # identify pairs as superficial, middle, deep, or other:
        c, supis, midis, deepis, otheris = self.pair_laminarity(self.nids, pairs)
        # get percentages of each:
        psup = intround(len(supis) / npairs * 100)
        pmid = intround(len(midis) / npairs * 100)
        pdeep = intround(len(deepis) / npairs * 100)
        pother = intround(len(otheris) / npairs * 100)
        
        a.scatter(np.arange(self.npairs), corrs, marker='o', c=c, edgecolor='none',
                  s=10, zorder=100)
        a.set_xlim(left=-10)
        a.set_ylim(bottom=-0.05)
        # underplot horizontal line at y=0:
        a.axhline(y=0, c='e', ls='--', marker=None)
        a.set_xlabel("pair index")
        a.set_ylabel("correlation coefficient")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        self.mean = np.mean(corrs)
        self.median = np.median(corrs)
        self.stdev = np.std(corrs)
        # add stuff to top right of plot:
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][self.r.tr.absname]
        a.text(0.99, 0.99, '%s\n'
                           'mean = %.3f\n'
                           'median = %.3f\n'
                           'stdev = %.3f\n'
                           'tres = %d ms\n'
                           'phase = %d deg\n'
                           'R = %r um\n'
                           'minrate = %.2f Hz\n'
                           'nneurons = %d\n'
                           'npairs = %d\n'
                           'sup = %r um\n'
                           'mid = %r um\n'
                           'deep = %r um\n'
                           'dt = %d min'
                           % (self.r.name, self.mean, self.median, self.stdev,
                              uns['CODETRES']//1000, uns['CODEPHASE'], self.R,
                              uns['MINRATE'], len(self.nids), self.npairs,
                              sup, mid, deep, intround(self.r.dtmin)),
                           transform = a.transAxes,
                           horizontalalignment='right',
                           verticalalignment='top')
        # make proxy line artists for legend:
        sl = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='r', mec='r')
        ml = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='g', mec='g')
        dl = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='b', mec='b')
        ol = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='y', mec='y')
        # add legend:
        a.legend([sl, ml, dl, ol],
                 ['superficial: %d%%' % psup, 'middle: %d%%' % pmid, 'deep: %d%%' % pdeep,
                  'other: %d%%' % pother],
                 numpoints=1, loc='upper center',
                 handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def scat(self, otherrid, nids=None, crange=[-0.05, 0.25], figsize=(7.5, 6.5)):
        """Scatter plot pairwise spike correlations in this recording vs that of
        another recording. If the two recordings are the same, split it in half and scatter
        plot first half against the second half."""
        ## TODO: for case when splitting a code array in half and comparing the two halves,
        ## extend this to work even when len(self.recs) > 1, as long as len(self.codes) == 1

        ## TODO: add interleave flag which generates a sufficiently interleaved, equally sized,
        ## non-overlapping set of tranges to scatter plot against each other, to eliminate
        ## temporal bias inherent in a simple split in time
        if len(self.recs) > 1:
            raise ValueError('for now, this method only works on single recordings')
        r0 = self.recs[0]
        tr = r0.tr
        otherr = tr.r[otherrid] # assume otherrid is from same track
        r1 = otherr
        # make sure they're from the same track, though the above should guarantee it
        assert r0.tr == r1.tr
        if r0 != r1:
            tranges0 = [r0.trange] # use the usual trange for both
            tranges1 = [r1.trange]
            xlabel = 'correlation coefficient (%s)' % r0.name
            ylabel = 'correlation coefficient (%s)' % r1.name
        else: # same recording, split its trange in half
            start, end = r0.trange
            mid = intround(start + r0.dt / 2)
            tranges0 = [(start, mid)]
            tranges1 = [(mid, end)]
            xlabel = 'correlation coefficient (%s, 1st half)' % r0.name
            ylabel = 'correlation coefficient (%s, 2nd half)' % r0.name
        if nids == None:
            if r0 != r1: # find nids active in both recordings
                nids = tr.get_nids([r0.id, r1.id])
            else: # same recording, find nids active during both tranges
                nids0 = r0.get_nids(tranges0)
                nids1 = r1.get_nids(tranges1)
                nids = np.intersect1d(nids0, nids1)

        # given the same nids, calculate corrs for both, constrained to tranges0
        # and tranges1 respectively, and to the torus described by R:
        sc0 = SpikeCorr([r0], tranges=tranges0, nids=nids)
        sc1 = SpikeCorr([r1], tranges=tranges1, nids=nids)
        sc0.calc()
        sc1.calc()
        # just to be sure:
        if sc0.npairs[0] != sc1.npairs[0] or (sc0.pairs[0] != sc1.pairs[0]).any():
            import pdb; pdb.set_trace()
            raise RuntimeError("sc0 and sc1 pairs don't match")
        pairs = sc0.pairs[0]
        npairs = len(pairs)
        corrs0, corrs1 = sc0.corrs[0], sc1.corrs[0]
        
        # identify pairs as superficial, middle, deep, or other:
        ypos = np.array([ r0.n[nid].pos[1] for nid in nids ])
        c, supis, midis, deepis, otheris = pair_laminarity(nids, ypos, tr.absname, pairs)
        # get percentages of each:
        psup = intround(len(supis) / npairs * 100)
        pmid = intround(len(midis) / npairs * 100)
        pdeep = intround(len(deepis) / npairs * 100)
        pother = intround(len(otheris) / npairs * 100)
        supcorrs0, supcorrs1 = corrs0[supis], corrs1[supis]
        midcorrs0, midcorrs1 = corrs0[midis], corrs1[midis]
        deepcorrs0, deepcorrs1 = corrs0[deepis], corrs1[deepis]
        othercorrs0, othercorrs1 = corrs0[otheris], corrs1[otheris]
        
        # create the scatter plot:
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        lim = crange
        if crange == None:
            # fit to nearest 0.05 encompassing all corr values on both axes:
            xlim = min(corrs0), max(corrs0)
            ylim = min(corrs1), max(corrs1)
            xlim = [ roundto(val, 0.05) for val in xlim ]
            ylim = [ roundto(val, 0.05) for val in ylim ]
            minlim = min(xlim[0], ylim[0])
            maxlim = max(xlim[1], ylim[1])
            lim = minlim, maxlim

        a.plot(lim, lim, c='e', ls='--', marker=None) # y=x line
        if psup > 0:
            a.errorbar(supcorrs0.mean(), supcorrs1.mean(),
                       xerr=supcorrs0.std(), yerr=supcorrs1.std(), color='r')
        if pmid > 0:
            a.errorbar(midcorrs0.mean(), midcorrs1.mean(),
                       xerr=midcorrs0.std(), yerr=midcorrs1.std(), color='g')
        if pdeep > 0:
            a.errorbar(deepcorrs0.mean(), deepcorrs1.mean(),
                       xerr=deepcorrs0.std(), yerr=deepcorrs1.std(), color='b')
        if pother > 0:
            a.errorbar(othercorrs0.mean(), othercorrs1.mean(),
                       xerr=othercorrs0.std(), yerr=othercorrs1.std(), color='y')
        a.scatter(corrs0, corrs1, marker='o', c=c, edgecolor='none', s=10, zorder=100)
        a.set_xlim(lim)
        a.set_ylim(lim)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        # add stuff to top left of plot:
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][tr.absname]
        a.text(0.01, 0.99, 'tres = %d ms\n'
                           'phase = %d deg\n'
                           'minrate = %.2f Hz\n'
                           'nneurons = %d\n'
                           'npairs = %d\n'
                           'sup = %r um\n'
                           'mid = %r um\n'
                           'deep = %r um\n'
                           'r%s.dt = %d min\n'
                           'r%s.dt = %d min'
                           % (uns['CODETRES']//1000, uns['CODEPHASE'], uns['MINRATE'],
                              len(nids), sc0.npairs[0],
                              sup, mid, deep,
                              r0.id, intround(r0.dtmin), r1.id, intround(r1.dtmin)),
                           transform = a.transAxes,
                           horizontalalignment='left',
                           verticalalignment='top')
        # make proxy artists for legend:
        sl = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='r', mec='r')
        ml = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='g', mec='g')
        dl = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='b', mec='b')
        ol = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='y', mec='y')
        # add legend:
        a.legend([sl, ml, dl, ol],
                 ['superficial: %d%%' % psup, 'middle: %d%%' % pmid, 'deep: %d%%' % pdeep,
                  'other: %d%%' % pother],
                 numpoints=1, loc='lower right',
                 handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def sep(self, figsize=(7.5, 6.5)):
        """Scatter plot pairwise spike correlations as a f'n of pair separation"""
        ## TODO: update for multiple tracks
        self.calc()
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        nids = self.nids
        corrs = self.corrs
        pairs = self.pairs
        npairs = len(pairs)
        n = self.r.n

        # identify pairs as superficial, middle, deep, or other:
        c, supis, midis, deepis, otheris = self.pair_laminarity(self.nids, pairs)
        # get percentages of each:
        psup = intround(len(supis) / npairs * 100)
        pmid = intround(len(midis) / npairs * 100)
        pdeep = intround(len(deepis) / npairs * 100)
        pother = intround(len(otheris) / npairs * 100)
        supcorrs = corrs[supis]
        midcorrs = corrs[midis]
        deepcorrs = corrs[deepis]
        othercorrs = corrs[otheris]

        # pairwise separations:
        seps = np.zeros(npairs)
        for i, pair in enumerate(pairs):
            nid0, nid1 = nids[pair[0]], nids[pair[1]]
            seps[i] = dist(n[nid0].pos, n[nid1].pos)
        supseps = seps[supis]
        midseps = seps[midis]
        deepseps = seps[deepis]
        otherseps = seps[otheris]

        if psup > 0:
            a.errorbar(supseps.mean(), supcorrs.mean(),
                       xerr=supseps.std(), yerr=supcorrs.std(), color='r', ls='--')
        if pmid > 0:
            a.errorbar(midseps.mean(), midcorrs.mean(),
                       xerr=midseps.std(), yerr=midcorrs.std(), color='g', ls='--')
        if pdeep > 0:
            a.errorbar(deepseps.mean(), deepcorrs.mean(),
                       xerr=deepseps.std(), yerr=deepcorrs.std(), color='b', ls='--')
        if pother > 0:
            a.errorbar(otherseps.mean(), othercorrs.mean(),
                       xerr=otherseps.std(), yerr=othercorrs.std(), color='y', ls='--')
        a.scatter(seps, corrs, marker='o', c=c, edgecolor='none', s=10, zorder=100)
        a.set_xlim(left=0)
        # underplot horizontal line at y=0:
        a.axhline(y=0, c='e', ls='--', marker=None, zorder=-100)
        a.set_xlabel("pair separation (um)")
        a.set_ylabel("spike correlation")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        # add stuff to top right of plot:
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][self.r.tr.absname]
        a.text(0.99, 0.99, '%s\n'
                           'tres = %d ms\n'
                           'phase = %d deg\n'
                           'R = %r um\n'
                           'minrate = %.2f Hz\n'
                           'nneurons = %d\n'
                           'npairs = %d\n'
                           'sup = %r um\n'
                           'mid = %r um\n'
                           'deep = %r um\n'
                           'dt = %d min'
                           % (self.r.name, uns['CODETRES']//1000, uns['CODEPHASE'], self.R,
                              uns['MINRATE'], len(self.nids), npairs,
                              sup, mid, deep, intround(self.r.dtmin)),
                           transform=a.transAxes,
                           horizontalalignment='right',
                           verticalalignment='top')
        # make proxy artists for legend:
        sl = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='r', mec='r')
        ml = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='g', mec='g')
        dl = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='b', mec='b')
        ol = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='y', mec='y')
        # add legend:
        a.legend([sl, ml, dl, ol],
                 ['superficial: %d%%' % psup, 'middle: %d%%' % pmid, 'deep: %d%%' % pdeep,
                  'other: %d%%' % pother],
                 numpoints=1, loc='upper center',
                 handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def sepbin(self, binwidth=100, figsize=(7.5, 6.5)):
        """Plot mean pairwise spike correlations as a f'n of binned pair separation"""
        ## TODO: update for multiple tracks
        self.calc()
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        nids = self.nids
        corrs = self.corrs
        pairs = self.pairs
        npairs = len(pairs)
        n = self.r.n

        # identify pairs as superficial, middle, deep, or other:
        c, supis, midis, deepis, otheris = self.pair_laminarity(self.nids, pairs)
        # get percentages of each:
        pall = intround(npairs / npairs * 100)
        psup = intround(len(supis) / npairs * 100)
        pmid = intround(len(midis) / npairs * 100)
        pdeep = intround(len(deepis) / npairs * 100)
        pother = intround(len(otheris) / npairs * 100)
        allcorrs = corrs
        supcorrs = allcorrs[supis]
        midcorrs = allcorrs[midis]
        deepcorrs = allcorrs[deepis]
        othercorrs = allcorrs[otheris]

        # pairwise separations:
        allseps = np.zeros(npairs)
        for i, pair in enumerate(pairs):
            nid0, nid1 = nids[pair[0]], nids[pair[1]]
            allseps[i] = dist(n[nid0].pos, n[nid1].pos)
        supseps = allseps[supis]
        midseps = allseps[midis]
        deepseps = allseps[deepis]
        otherseps = allseps[otheris]

        #a.scatter(allseps, corrs, marker='o', c=c, edgecolor='none', s=10, zorder=100)
        depthinfo = {'other': (otherseps, othercorrs, 'y', 'other: %d%%' % pother),
                     'sup': (supseps, supcorrs, 'r', 'superficial: %d%%' % psup),
                     'mid': (midseps, midcorrs, 'g', 'middle: %d%%' % pmid),
                     'deep': (deepseps, deepcorrs, 'b', 'deep: %d%%' % pdeep),
                     'all': (allseps, allcorrs, 'e', 'all: %d%%' % pall)}

        for depth in ('other', 'sup', 'mid', 'deep', 'all'):
            seps, corrs, c, label = depthinfo[depth]
            if len(seps) == 0:
                continue # skip this depth type
            edges = np.arange(0, seps.max()+binwidth, binwidth) # incl left and right edges
            midedges = edges[1:-1] # exclude left and right edges
            nbins = len(edges) - 1
            means = np.zeros(nbins)
            stdevs = np.zeros(nbins)
            #sems = np.zeros(nbins)
            # bini == 0: left of leftmost midedge; bini == nbins: right of rightmost midedge
            binis = np.digitize(seps, midedges)
            for bini in range(nbins):
                bincorrs = corrs[binis == bini]
                n = len(bincorrs)
                if n == 0:
                    continue
                means[bini] = bincorrs.mean()
                stdevs[bini] = bincorrs.std()
                #sems[bini] = stdevs[bini] / np.sqrt(n)
            x = (edges[1:] + edges[:-1]) / 2 # bin midpoints
            # don't plot empty bins, for safety, only exclude if both mean and stdev are 0:
            keepis = (means != 0.0) + (stdevs != 0.0)
            x, means, stdevs = x[keepis], means[keepis], stdevs[keepis]
            a.errorbar(x, means, yerr=stdevs, marker='o', c=c, label=label, mec='none')
            
        a.set_xlim(left=0)
        # underplot horizontal line at y=0:
        a.axhline(y=0, c='e', ls='--', marker=None, zorder=-100)
        a.set_xlabel("pair separation (um)")
        a.set_ylabel("mean spike correlation")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        # add stuff to top right of plot:
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][self.r.tr.absname]
        a.text(0.99, 0.99, '%s\n'
                           'tres = %d ms\n'
                           'phase = %d deg\n'
                           'R = %r um\n'
                           'minrate = %.2f Hz\n'
                           'nneurons = %d\n'
                           'npairs = %d\n'
                           'sup = %r um\n'
                           'mid = %r um\n'
                           'deep = %r um\n'
                           'dt = %d min\n'
                           'binwidth = %d um'
                           % (self.r.name, uns['CODETRES']//1000, uns['CODEPHASE'], self.R,
                              uns['MINRATE'], len(self.nids), npairs,
                              sup, mid, deep, intround(self.r.dtmin), binwidth),
                           transform=a.transAxes,
                           horizontalalignment='right',
                           verticalalignment='top')
        # add legend:
        a.legend(numpoints=1, loc='upper center',
                 handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self

    def pos(self, maxsep=150, figsize=(7.5, 6.5)):
        """Plot spike corrs between cell pairs that fall within a threshold separation
        distance maxsep of each other, as a function of position down length of probe"""
        ## TODO: update for multiple tracks
        self.calc()
        assert len(np.unique(self.pairs)) == len(self.nids) # sanity check
        n = self.r.n
        # use just y position of each neuron:
        pos = np.array([n[nid].pos[1] for nid in self.nids]) # nneurons position array (um)
        pos.shape = -1, 1 # make it 2D for pdist
        sep = pdist(pos) # 1D vector form, requires pairwise combinatorial indexing
        #sep = squareform(sep) # matrix form, with simple 2D indexing
        
        # keep only those pairs within maxsep separation:
        pairis = sep <= maxsep # boolean array
        pairs = self.pairs[pairis]
        sep = sep[pairis] # keep only the relevant pair separation values
        npairs = len(pairs)
        pairpos = np.zeros((npairs, 2))
        for pairi, pair in enumerate(pairs):
            nid0, nid1 = pair
            pairpos[pairi] = np.mean([pos[nid0], pos[nid1]], axis=0) # mean position of pair
        corrs = self.corrs[pairis]
        
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        ypos = pairpos[:, 1]
        # this is only useful when connecting the dots using plot():
        '''
        sortis = ypos.argsort()
        ypos = ypos[sortis]
        corrs = corrs[sortis]
        sep = sep[sortis]
        '''
        # underplot horizontal line at y=0:
        a.axhline(y=0, c='e', ls='--', marker=None)

        # scatter plot corrs vs ypos, black=0 um separation, white=maxsep um separation,
        # could also use cmap=mpl.cm.jet_r instead:
        a.scatter(ypos, corrs, c=sep, vmin=0, vmax=maxsep, cmap=mpl.cm.gray,
                  marker='.', s=100, lw=0.5)

        minpos = 0 #min(self.r.chanpos[:, 1])
        maxpos = max(self.r.chanpos[:, 1])
        a.set_xlim((minpos, maxpos))
        ymin, ymax = a.get_ylim()
        ymin = min(ymin, -0.05) # set to no more than this
        ymax = max(ymax, 0.3) # set to no less than this
        a.set_ylim((ymin, ymax))

        a.set_xlabel("vertical pair position (um)")
        a.set_ylabel("spike correlation")
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)

        a.text(0.995, 0.99, '%s\n'
                            'white = %d um maxsep\n'
                            'npairs = %d'
               % (self.r.name, maxsep, npairs), color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        #return ypos, corrs
        return self

    def sct(self, method='mean', inclusive=False):
        """Calculate pairwise spike correlations for each type of laminarity as a function of
        time. method can be 'weighted mean', 'mean', 'median', 'max', 'min' or 'all'"""
        ## can this work over multiple tracks (ie multiple code arrays), or do I need to
        ## enforce only a single track?
        uns = get_ipython().user_ns
        if self.width == None:
            self.width = intround(uns['SCWIDTH'] * 1000000) # convert from sec to us
        if self.tres == None:
            self.tres = intround(uns['SCTRES'] * 1000000) # convert from sec to us
        self.calc()
        allis = np.arange(self.npairs) # all indices into self.pairs
        c, supis, midis, deepis, otheris = self.pair_laminarity(self.nids, self.pairs,
                                                                inclusive=inclusive)
        laminarcorrs = []
        laminarnpairs = []
        for pairis in (allis, supis, midis, deepis, otheris):
            npairs = len(pairis)
            if npairs == 0:
                laminarcorrs.append(np.zeros(len(self.tranges)))
                laminarnpairs.append(npairs)
                continue
            corrs = self.corrs[pairis] # npairs x ntranges
            if method.startswith('weighted'):
                # weight each pair by its normalized ON count per trange
                counts = self.counts[pairis] # npairs * ntranges
                totalcounts = counts.sum(axis=0) # len(ntranges)
                # avoid div by 0, counts at such timepoints will be uniformly 0 anyway:
                zcountis = totalcounts == 0 # trange indices where totalcounts are 0
                totalcounts[zcountis] = 1
                weights = counts / totalcounts # npairs x ntranges
            if method == 'weighted mean':
                # sum over all weighted pairs:
                corrs = (corrs * weights).sum(axis=0) # len(ntranges)
            elif method == 'weighted median': # not entirely sure this is right:
                corrs = np.median(corrs * weights, axis=0) * self.npairs
            elif method == 'mean':
                corrs = corrs.mean(axis=0)
            elif method == 'median':
                corrs = np.median(corrs, axis=0)
            elif method == 'max':
                corrs = corrs.max(axis=0)
            elif method == 'min':
                corrs = corrs.min(axis=0)
            elif method == 'all':
                corrs = corrs.T # need transpose for some reason when plotting multiple traces
            else:
                raise ValueError("unknown method %r" % method)
            laminarcorrs.append(corrs)
            laminarnpairs.append(npairs)
        laminarcorrs = np.vstack(laminarcorrs)
        laminarnpairs = np.array(laminarnpairs)
        # get midpoint of each trange, convert from us to sec:
        t = self.tranges.mean(axis=1) / 1000000
        ylabel = method + ' spike correlations'
        return laminarcorrs, laminarnpairs, t, ylabel

    def plot(self, method='mean', inclusive=False, figsize=(20, 6.5)):
        """Plot pairwise spike correlations as a function of time. method can be
        'weighted mean', 'mean', 'median', 'max' or 'min'"""
        ## can this work over multiple tracks (ie multiple code arrays), or do I need to
        ## enforce only a single track?

        corrs, npairs, t, ylabel = self.sct(method=method, inclusive=inclusive)
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        # underplot horizontal line at y=0:
        a.axhline(y=0, c='e', ls='--', marker=None)
        #if corrs.ndim == 2:
        #    a.plot(t, corrs) # auto colours
        # plot according to laminarity:
        a.plot(t, corrs[0], 'e.-', label='all (%d)' % npairs[0])
        a.plot(t, corrs[1], 'r.-', label='superficial (%d)' % npairs[1])
        a.plot(t, corrs[2], 'g.-', label='middle (%d)' % npairs[2])
        a.plot(t, corrs[3], 'b.-', label='deep (%d)' % npairs[3])
        a.plot(t, corrs[4], 'y.-', label='other (%d)' % npairs[4], zorder=0)
        a.set_xlabel("time (sec)")
        ylabel = ylabel + " (%d pairs)" % self.npairs
        a.set_ylabel(ylabel)
        # limit plot to duration of acquistion, in sec:
        t0, t1 = np.asarray(self.r.trange) / 1000000
        ymax = max([0.1, corrs[[0,1,2,3]].max()])
        ymin = min([0.0, corrs[[0,1,2,3]].min()])
        a.set_ylim(ymin, ymax)
        a.set_xlim(t0, t1)
        #a.autoscale(axis='x', enable=True, tight=True)
        # turn off annoying "+2.41e3" type offset on x axis:
        formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        a.xaxis.set_major_formatter(formatter)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][self.r.tr.absname]
        a.text(0.998, 0.99,
               '%s\n'
               'sup = %r um\n'
               'mid = %r um\n'
               'deep = %r um'
               % (self.r.name, sup, mid, deep),
               color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents

    def si(self, method='mean', inclusive=False, sisource='lfp', kind=None, chani=-1,
           sirange=None, plot=True, layers=False, ms=5, figsize=(7.5, 6.5)):
        """Scatter plot spike correlations vs MUA or LFP synchrony index"""
        ## TODO: update for multiple recs
        rec = self.r
        uns = get_ipython().user_ns
        #t0 = time.time()

        if sisource not in ['lfp', 'mua']:
            raise ValueError('unknown sisource %r' % sisource)

        if kind == None:
            if sisource == 'lfp':
                kind = uns['LFPSIKIND']
            else:
                kind = uns['MUASIKIND']

        if layers == False:
            layers = ['all']
        elif layers == True:
            layers = ['sup', 'deep']
        LAYER2I = {'all':0, 'sup':1, 'mid':2, 'deep':3, 'other':4}
        layeris = [ LAYER2I[layer] for layer in layers ]

        # ct are center timepoints:
        corrs, npairs, ct, ylabel = self.sct(method=method, inclusive=inclusive)
        #print('sct(t) calc took %.3f sec' % (time.time()-t0))
        # sit are also center timepoints:
        if sisource == 'lfp':
            si, sit = rec.lfp.si(chani=chani, kind=kind, plot=False)
            si = np.vstack([si, si, si, si]) # make 4 x nt, just like for mua si
        else:
            si, sit, n = rec.mua_si(kind=kind, plot=False)
        #print('SI(t) calc took %.3f sec' % (time.time()-t0))
        #print('len(sit, ct) = %d, %d' % (len(sit), len(ct)))

        # get common time resolution:
        t, corrs, si = commontres(ct, corrs, sit, si)

        if not plot:
            return corrs, si, ylabel # return corrs and si of all laminarities

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)

        #ylim = corrs[layeris].min(), corrs[layeris].max()
        #yrange = ylim[1] - ylim[0]
        #extra = yrange*0.03 # 3 %
        #ylim = ylim[0]-extra, ylim[1]+extra
        ylim = uns['SCLIMITS']

        # keep only those points whose synchrony index falls within sirange:
        if sirange == None:
            finitesi = si[np.isfinite(si)]
            sirange = finitesi.min(), finitesi.max()
        sirange = np.asarray(sirange)
        keepis = (sirange[0] <= si[0]) * (si[0] <= sirange[1]) # boolean index array
        si = si[:, keepis]
        corrs = corrs[:, keepis]

        # plot linear regressions of corrs vs si[0]:
        if 'all' in layers:
            m0, b0, r0, p0, stderr0 = linregress(si[0], corrs[0])
            a.plot(sirange, m0*sirange+b0, 'e--')
        if 'sup' in layers:
            m1, b1, r1, p1, stderr1 = linregress(si[0], corrs[1])
            a.plot(sirange, m1*sirange+b1, 'r--')
        if 'mid' in layers:
            m2, b2, r2, p2, stderr2 = linregress(si[0], corrs[2])
            a.plot(sirange, m2*sirange+b2, 'g--')
        if 'deep' in layers:
            m3, b3, r3, p3, stderr3 = linregress(si[0], corrs[3])
            a.plot(sirange, m3*sirange+b3, 'b--')
        if 'other' in layers:
            m4, b4, r4, p4, stderr4 = linregress(si[0], corrs[4])
            a.plot(sirange, m4*sirange+b4, 'y--', zorder=0)

        # scatter plot corrs vs si[0], one colour per laminarity:
        if 'all' in layers:
            a.plot(si[0], corrs[0], 'e.', ms=ms, label='all (%d), m=%.3f, r=%.3f'
                                                       % (npairs[0], m0, r0))
        if 'sup' in layers:
            a.plot(si[0], corrs[1], 'r.', ms=ms, label='superficial (%d), m=%.3f, r=%.3f'
                                                       % (npairs[1], m1, r1))
        if 'mid' in layers:
            a.plot(si[0], corrs[2], 'g.', ms=ms, label='middle (%d), m=%.3f, r=%.3f'
                                                       % (npairs[2], m2, r2))
        if 'deep' in layers:
            a.plot(si[0], corrs[3], 'b.', ms=ms, label='deep (%d), m=%.3f, r=%.3f'
                                                       % (npairs[3], m3, r3))
        if 'other' in layers:
            a.plot(si[0], corrs[4], 'y.', ms=ms, label='other (%d), m=%.3f, r=%.3f'
                                                       % (npairs[4], m4, r4), zorder=0)
        #a.set_xlim(sirange)
        if kind[0] == 'n':
            a.set_xlim(-1, 1)
        a.set_ylim(ylim)
        #a.autoscale(enable=True, axis='y', tight=True)
        a.set_xlabel('%s SI (%s)' % (sisource.upper(), kind))
        ylabel = ylabel + " (%d pairs)" % self.npairs
        a.set_ylabel(ylabel)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        sup, mid, deep = uns['LAYERS'][rec.tr.absname]
        a.text(0.998, 0.99,
               '%s\n'
               'sup = %r um\n'
               'mid = %r um\n'
               'deep = %r um'
               % (rec.name, sup, mid, deep),
               color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents

    def mua(self, method='mean', inclusive=False, smooth=False, figsize=(7.5, 6.5)):
        """Scatter plot spike correlations vs multiunit activity"""
        ## TODO: update for multiple recs
        corrs, npairs, ct, ylabel = self.sct(method=method, inclusive=inclusive)
        mua, muat, n = self.r.mua(smooth=smooth, plot=False)
        # keep only MUA of all neurons, throw away laminar MUA information (for now at least):
        mua = mua[0] # 1D array

        # get common time resolution:
        t, corrs, mua = commontres(ct, corrs, muat, mua)

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        ylim = corrs[:5].min(), corrs[:5].max()
        yrange = ylim[1] - ylim[0]
        extra = yrange*0.03 # 3 %
        ylim = ylim[0]-extra, ylim[1]+extra
        muarange = np.array([mua.min(), mua.max()])

        # plot linear regressions:
        m0, b0, r0, p0, stderr0 = linregress(mua, corrs[0])
        m1, b1, r1, p1, stderr1 = linregress(mua, corrs[1])
        m2, b2, r2, p2, stderr2 = linregress(mua, corrs[2])
        m3, b3, r3, p3, stderr3 = linregress(mua, corrs[3])
        m4, b4, r4, p4, stderr4 = linregress(mua, corrs[4])
        a.plot(muarange, m0*muarange+b0, 'e--')
        a.plot(muarange, m1*muarange+b1, 'r--')
        a.plot(muarange, m2*muarange+b2, 'g--')
        a.plot(muarange, m3*muarange+b3, 'b--')
        a.plot(muarange, m4*muarange+b4, 'y--', zorder=0)

        # scatter plot corrs vs mua, one colour per laminarity:
        a.plot(mua, corrs[0], 'e.', label='all (%d), m=%.3f, r=%.3f' % (npairs[0], m0, r0))
        a.plot(mua, corrs[1], 'r.', label='superficial (%d), m=%.3f, r=%.3f'
                                          % (npairs[1], m1, r1))
        a.plot(mua, corrs[2], 'g.', label='middle (%d), m=%.3f, r=%.3f' % (npairs[2], m2, r2))
        a.plot(mua, corrs[3], 'b.', label='deep (%d), m=%.3f, r=%.3f' % (npairs[3], m3, r3))
        a.plot(mua, corrs[4], 'y.', label='other (%d), m=%.3f, r=%.3f'
                                          % (npairs[4], m4, r4), zorder=0)
        a.set_ylim(ylim)
        #a.autoscale(enable=True, axis='y', tight=True)
        a.set_xlabel("mean MUA (Hz), %d neurons" % n[0])
        ylabel = ylabel + " (%d pairs)" % self.npairs
        a.set_ylabel(ylabel)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        uns = get_ipython().user_ns
        sup, mid, deep = uns['LAYERS'][self.r.tr.absname]
        a.text(0.998, 0.99,
               '%s\n'
               'sup = %r um\n'
               'middle = %r um\n'
               'deep = %r um'
               % (self.r.name, sup, mid, deep),
               color='k', transform=a.transAxes,
               horizontalalignment='right', verticalalignment='top')
        a.legend(loc='upper left', handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents


class NeuropyWindow(QtGui.QMainWindow):
    """Base class for all of neuropy's tool windows"""
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.maximized = False

    def keyPressEvent(self, event):
        key = event.key()
        if key == Qt.Key_F11:
            self.toggleMaximized()
        else:
            QtGui.QMainWindow.keyPressEvent(self, event) # pass it on

    def mouseDoubleClickEvent(self, event):
        """Doesn't catch window titlebar doubleclicks for some reason (window manager
        catches them?). Have to doubleclick on a part of the window with no widgets in it"""
        self.toggleMaximized()

    def toggleMaximized(self):
        if not self.maximized:
            self.normalPos, self.normalSize = self.pos(), self.size()
            dw = QtGui.QDesktopWidget()
            rect = dw.availableGeometry(self)
            self.setGeometry(rect)
            self.maximized = True
        else: # restore
            self.resize(self.normalSize)
            self.move(self.normalPos)
            self.maximized = False


class RevCorrWindow(NeuropyWindow):
    def __init__(self, parent=None, title='RevCorrWindow', rfs=None,
                 nids=None, ts=None, scale=2.0, blanksize=(32, 32)):
        NeuropyWindow.__init__(self, parent)
        self.title = title
        self.rfs = rfs
        self.nids = nids
        self.ts = ts
        self.scale = scale # setting to non-integer will give uneven sized pixels

        cmap = mpl.cm.jet(np.arange(256), alpha=None, bytes=True) # 8 bit RGBA colormap
        # from Qt docs, need to use ARGB format:
        # http://qt-project.org/doc/qt-4.8/qimage.html#image-formats
        # convert to 8 bit ARGB colormap, but due to little-endianness, need to arrange
        # array columns in reverse BGRA order:
        cmap[:, [0, 1, 2, 3]] = cmap[:, [2, 1, 0, 3]]
        colortable = cmap.view(dtype=np.uint32).ravel().tolist() # QVector<QRgb> colors
        layout = QtGui.QGridLayout() # can set vert and horiz spacing
        #layout.setContentsMargins(0, 0, 0, 0) # doesn't seem to do anything

        # place time labels along top
        for ti, t in enumerate(ts):
            label = QtGui.QLabel(str(t))
            layout.addWidget(label, 0, ti+1)
        # plot each row, with its nid label
        for ni, nid in enumerate(nids):
            label = QtGui.QLabel('n'+str(nid)) # nid label on left
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            layout.addWidget(label, ni+1, 0)
            rf = rfs[ni]
            for ti, t in enumerate(ts):
                #data = np.uint8(np.random.randint(0, 255, size=(height, width)))
                if rf == None:
                    data = np.zeros(blanksize, dtype=np.uint8) # blank placeholder rf
                else:
                    data = rf[ti]
                width, height = data.shape
                image = QImage(data.data, width, height, QImage.Format_Indexed8)
                image.ndarray = data # hold a ref, prevent gc
                image.setColorTable(colortable)
                image = image.scaled(QSize(scale*width, scale*height)) # scale it
                pixmap = QPixmap.fromImage(image)
                label = QtGui.QLabel()
                label.setPixmap(pixmap)
                layout.addWidget(label, ni+1, ti+1) # can also control alignment

        mainwidget = QtGui.QWidget(self)
        mainwidget.setLayout(layout)

        scrollarea = QtGui.QScrollArea()
        scrollarea.setWidget(mainwidget)

        self.setCentralWidget(scrollarea)
        self.setWindowTitle(title)
        #palette = QPalette(QColor(255, 255, 255))
        #self.setPalette(palette) # set white background, or perhaps more

def mplrevcorr(title='RevCorrWindow', rfs=None, nids=None, ts=None, scale=2, dpi=100,
               margins=True):
    """MPL version of RevCorrWindow, good for saving RFs to a file. This one uses figimage to
    prevent image resampling. It seems that for this to save to file correctly, you need to
    explicitly set the output file to have the same DPI as the figure:

    gcf().savefig('filename.png', dpi=gcf().dpi)
    You might also want to set transparency=False kwarg to get rid of transparency, or set
    that as the default in rcParams['savefig.transparent'].
    
    Saving to PDF seems buggy. Use PNG instead.

    Do bbox_inches='tight', pad_inches=0 do anything? Don't seem to...
    """
    # spacing (inches):
    lm, rm = 0.5*margins, 0 # left and right margins
    tm, bm = 0.2*margins, 0 # top and bottom margins
    vs, hs = 0.04, 0.04 # vertical and horizontal spacing
    w, h = np.asarray(rfs[0][0].shape) * scale # in pixels
    #w, h = 32*scale, 32*scale # data width and height, assumed for now
    #assert w == h
    # number of inches in to get one screen pixel per image pixel:
    if dpi == True:
        dpi = 100
    winch = w / dpi # width in inches
    hinch = h / dpi # height in inches
    rfw, rfh = winch, hinch # rf width and height, in inches

    # size the figure:
    nn = len(nids)
    nt = len(ts)
    fw = lm + nt*rfw + (nt-1)*hs + rm # inches
    fh = tm + nn*rfh + (nn-1)*vs + bm # inches
    #print(fw, fh)
    f = plt.figure(figsize=(fw, fh), dpi=dpi)

    maxni = nn - 1
    # place time labels along top:
    y = (bm + nn*(rfh + vs)) / fh
    if margins:
        for ti, t in enumerate(ts):
            x = (lm + ti*(rfw + hs)) / fw
            plt.figtext(x, y, str(t)+' ms')
    # plot each row, with its nid label
    x = (lm - hs) / fw # right edge of text measured from left, normalized
    ims = []
    for ni, nid in enumerate(nids):
        rf = rfs[ni]
        # center of text measured from bottom, normalized:
        y = (bm + (maxni-ni)*(rfh+vs) + rfh/2) / fh
        if margins:
            plt.figtext(x, y, 'n'+str(nid),
                        verticalalignment='center', horizontalalignment='right')
        if rf == None:
            continue # skip blank rfs
        for ti, t in enumerate(ts):
            # x and y pos of rf, in pixels:
            x0 = (lm + ti*(rfw+hs)) * dpi
            y0 = (bm + (maxni-ni)*(rfh+vs)) * dpi
            #data = np.uint8(np.random.randint(0, 255, size=(32, 32)))
            data = rf[ti]
            if scale != 1:
                data = data.repeat(scale, axis=0).repeat(scale, axis=1)
            im = plt.figimage(data, x0, y0, origin='upper', cmap=mpl.cm.jet,
                              vmin=0, vmax=255)
            ims.append(im)
            #a.set_axis_off() # disable axes lines, ticks and labels

    gcfm().window.setWindowTitle(title)
    return ims

def mplrevcorraxes(title='RevCorrWindow', rfs=None, nids=None, ts=None, scale=2.0):
    """MPL version of RevCorrWindow, good for saving RFs to a file. This one uses one axes per
    RF, which is a bit slow. If the figure is sized exactly right, then there should be no
    resampling effects of the images. It seems that for this to save to file correctly, you
    need to explicitly set the output file to have the same DPI as the figure:

    gcf().savefig('filename.pdf', dpi=gcf().dpi)
    """
    # spacing (inches):
    lm, rm = 0, 0 # left and right margins
    tm, bm = 0, 0 # top and bottom margins
    vs, hs = 0.05, 0.05 # vertical and horizontal spacing
    # default is 80 dpi, and default RF image size is 32x32, so 32/80 = 0.4 inches needed
    # for 1 screen pixel per image pixel. So size rfw and rfh by multiple of 0.4?
    w, h = rfs[0][0].shape
    #w, h = 32, 32 # data width and height, assumed for now
    assert w == h
    # number of inches in to get one screen pixel per image pixel:
    basesize = w / mpl.rcParams['figure.dpi']
    rfw, rfh = basesize*scale, basesize*scale # rf width and height

    # size the figure:
    nn = len(nids)
    nt = len(ts)
    fw = lm + nt*rfw + (nt-1)*hs + rm # inches
    fh = tm + nn*rfh + (nn-1)*vs + bm # inches
    f = plt.figure(figsize=(fw, fh))

    # plot each row:
    maxni = nn - 1
    for ni, nid in enumerate(nids):
        rf = rfs[ni]
        for ti, t in enumerate(ts):
            # left bottom width and height of rf, in fractional figure units:
            lbwh = (lm + ti*(rfw+hs))/fw, (bm + (maxni-ni)*(rfh+vs))/fh, rfw/fw, rfh/fh
            a = f.add_axes(lbwh, axisbg='w')
            a.set_axis_off() # disable axes lines, ticks and labels
            data = rf[ti]
            #data = np.uint8(np.random.randint(0, 255, size=(32, 32)))
            im = a.imshow(data, interpolation='nearest', cmap=mpl.cm.jet,
                          vmin=0, vmax=255, aspect='equal', origin='upper',
                          extent=[0, data.shape[1], 0, data.shape[0]])

    gcfm().window.setWindowTitle(title)

'''
class NetstateReceptiveFieldFrame(ReceptiveFieldFrame):
    """A wx.Frame for plotting a scrollable 2D grid of netstate receptive fields, with
    netstate and time labels. rfs is a list of (nt, width, height) sized receptive fields of
    uint8 RGB data, one per netstate"""
    def __init__(self, parent=None, id=-1, title='NetstateReceptiveFieldFrame',
                 rfs=None, intcodes=None, t=None, scale=2.0):
        self.rfs = rfs
        self.intcodes = tolist(intcodes)
        self.t = t
        self.title = title
        wx.Frame.__init__(self, parent=parent, id=id, title=title,
                          style=wx.DEFAULT_FRAME_STYLE)
        self.panel = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)
        self.bitmaps = {}
        for ii, i in enumerate(self.intcodes):
            self.bitmaps[ii] = {}
            for ti, t in enumerate(self.t):
                rf = self.rfs[ii][ti]
                # expose rf as databuffer:
                im = wx.ImageFromData(width=rf.shape[0], height=rf.shape[1], data=rf.data)
                im = im.Scale(width=im.GetWidth()*scale, height=im.GetHeight()*scale)
                self.bitmaps[ii][t] = wx.StaticBitmap(parent=self.panel,
                                                      bitmap=im.ConvertToBitmap())
        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle(self.title)
        self.panel.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.panel.SetScrollRate(10, 10)

    def __do_layout(self):
        sizer_1 = wx.GridSizer(1, 1, 0, 0)
        # add an extra row and column for the text labels:
        grid_sizer_1 = wx.FlexGridSizer(rows=len(self.intcodes)+1, cols=len(self.t)+1,
                                        vgap=2, hgap=2)
        grid_sizer_1.Add((1, 1), 0, wx.ADJUST_MINSIZE, 0) # spacer in top left corner
        for t in self.t:
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "%sms" % t), 0,
                             wx.ADJUST_MINSIZE, 0) # text row along top
        for ii, i in enumerate(self.intcodes):
            # text down left side:
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "ns%d" % i), 0,
                wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ADJUST_MINSIZE, 0)
            for t in self.t:
                grid_sizer_1.Add(self.bitmaps[ii][t], 1,
                                 wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL, 0)
        self.panel.SetAutoLayout(True)
        self.panel.SetSizer(grid_sizer_1)
        grid_sizer_1.Fit(self.panel)
        sizer_1.Add(self.panel, 1, wx.ADJUST_MINSIZE|wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
'''

class Ising(object):
    """Maximum entropy Ising model"""
    def __init__(self, means, pairmeans, algorithm='CG'):
        """means is a list of mean activity values [-1 to 1] for each neuron code.
        pairmeans is list of products of activity values for all pairs of neuron codes.
        'Returns a maximum-entropy (exponential-form) model on a discrete sample space'
            -- scipy.maxent.model
        """
        from scipy import maxentropy

        # why are the abs values of these so big, shouldn't they be near 0?:
        #print('means:\n', means)
        #print('pairmeans:\n', pairmeans)

        nbits = len(means)
        npairs = len(pairmeans)
        assert npairs == nCr(nbits, 2) # sanity check
        self.intsamplespace = range(0, 2**nbits)
        table = getbinarytable(nbits=nbits) # words are in the columns, MSB at bottom row
        # all possible binary words (each is MSB to LSB), as arrays of 0s and 1s:
        self.binsamplespace = np.array([ table[::-1, wordi] for wordi in range(0, 2**nbits) ])
        self.samplespace = self.binsamplespace * 2 - 1 # convert 0s to -1s
        # return the i'th bit (LSB to MSB) of binary word x,
        # have to do i=i to statically assign value (gets around scope closure problem):
        f1s = [ lambda x, i=i: x[-1-i] for i in range(0, nbits) ]
        # return product of the i'th and j'th bit (LSB to MSB) of binary word
        f2s = []
        pairmeansi = 0
        for i in range(0, nbits):
            for j in range(i+1, nbits):
                if pairmeans[pairmeansi] != None: # None indicates we should ignore this pair
                    f2s.append(lambda x, i=i, j=j: x[-1-i] * x[-1-j])
                pairmeansi += 1
        #f2s = [ lambda x, i=i, j=j: x[-1-i] * x[-1-j] for i in range(0, nbits)
        #        for j in range(i+1, nbits) if pairmeans[i*nbits+j-1] != None ]
        f = np.concatenate((f1s, f2s))
        self.model = maxentropy.model(f, self.samplespace)
        #self.model.mindual = -10000
        #self.model.log = None # needed to make LBFGSB algorithm work
        # Now set the desired feature expectations
        means = np.asarray(means)
        pairmeans = np.asarray(pairmeans) # if it has Nones, it's an object array
        # remove the Nones, convert to list to get rid of object array, then convert
        # back to array to get a normal, non-object array (probably a float64 array):
        pairmeans = np.asarray(list(pairmeans[pairmeans != [None]]))
        npairs = len(pairmeans) # update npairs
        # add the one half in front of each coefficient, not too sure if
        # this should go here! causes convergence problems:
        #pairmeans /= 2.0
        K = np.concatenate((means, pairmeans))
        self.model.verbose = False

        # Fit the model
        """this has problems in scipy 0.9.0 and 0.10.1 but not 0.5.2. Raises error:
        
        "ValueError: operands could not be broadcast together with shapes (1024) (55)"

        This has to do with the call in model.logpmf() in maxentropy.py of maxentutils.py's
        innerprodtranspose(), which has this code in 0.5.2 around line 361:

            elif sparse.isspmatrix(A):
        return A.rmatvec(v).transpose()

        but was changed to this code in 0.9.0 and later:

        elif sparse.isspmatrix(A):
            return (A.conj().transpose() * v).transpose()

        which I guess returns an array of 1024 instead of the expected 55 (10 means plus 10
        choose 2 == 45 pairmeans). So maybe there's something wrong with this change that
        was made in scipy 0.9.0. Don't know about scipy 0.6.0 through 0.8.0 because I can't
        get any of those to compile without errors.
        """
        self.model.fit(K, algorithm=algorithm)

        self.hi = self.model.params[0:nbits]
        self.Jij = self.model.params[nbits:nbits+npairs]
        self.p = self.model.probdist()
        # sanity checks:
        assert (len(self.hi), len(self.Jij), len(self.p)) == (nbits, npairs, 2**nbits)
        #print('means:', means)
        #print('pairmeans:', pairmeans)
        #print('%d iters,' % self.model.iters)
        #print('hi:', self.hi.__repr__())
        #print('Jij:', self.Jij.__repr__())

        '''
        # Output the distribution
        print("\nFitted model parameters are:\n" + str(self.model.params))
        print("\nFitted distribution is:")
        for j in range(len(self.model.samplespace)):
            x = np.array(self.model.samplespace[j])
            x = (x+1)/2 # convert from -1s and 1s back to 0s and 1s
            print('\tx:%s, p(x):%s' % (x, p[j]))
        '''
        '''
        # Now show how well the constraints are satisfied:
        print()
        print("Desired constraints:")
        print("\tp['dans'] + p['en'] = 0.3")
        print("\tp['dans'] + p['" + a_grave + "']  = 0.5".encode('utf-8'))
        print()
        print("Actual expectations under the fitted model:")
        print("\tp['dans'] + p['en'] =", p[0] + p[1])
        print(("\tp['dans'] + p['" + a_grave + "']  = " + str(p[0]+p[2])).encode('utf-8'))
        # (Or substitute "x.encode('latin-1')" if you have a primitive terminal.)
        '''

class NeuropyScalarFormatter(mpl.ticker.ScalarFormatter):
    """Overloaded from mpl.ticker.ScalarFormatter for 4 reasons:
    1) turn off stupid offset
    2) increase maximum possible number of sigfigs
    3) increase +ve and -ve order of magnitude thresholds before switching to scientific
       notation
    4) keep exponents in engineering notation, ie multiples of 3
    """
    def __init__(self, useOffset=False, useMathText=False):
        # useOffset allows plotting small data ranges with large offsets:
        # for example: [1+1e-9,1+2e-9,1+3e-9]
        # useMathText will render the offset an scientific notation in mathtext
        # can't use this, because derived from an old-style class:
        #super(NeuropyScalarFormatter, self).__init__(useOffset=useOffset,
        #                                             useMathText=useMathText)
        mpl.ticker.ScalarFormatter.__init__(self, useOffset=useOffset, useMathText=useMathText)
        self.thousandsSep = '' # default to not using a thousands separator

    def _set_orderOfMagnitude(self, range):
        # if scientific notation is to be used, find the appropriate exponent
        # if using an numerical offset, find the exponent after applying the offset
        locs = np.absolute(self.locs)
        if self.offset: oom = math.floor(math.log10(range))
        else:
            if locs[0] > locs[-1]: val = locs[0]
            else: val = locs[-1]
            if val == 0: oom = 0
            else: oom = math.floor(math.log10(val))
        if oom < -3: # decreased -ve threshold for sci notation
            # stick to engineering notation, multiples of 3:
            self.orderOfMagnitude = (oom // 3)*3
        elif oom > 6: # increased +ve threshold for sci notation
            # stick to engineering notation, multiples of 3:
            self.orderOfMagnitude = (oom // 3)*3
        else:
            self.orderOfMagnitude = 0

    def _set_format(self):
        # set the format string to format all the ticklabels
        locs = (np.array(self.locs)-self.offset) / 10**self.orderOfMagnitude+1e-15
        # '%1.3f' changed to '%1.10f' to increase maximum number of possible sigfigs:
        sigfigs = [len(str('%1.10f'% loc).split('.')[1].rstrip('0')) for loc in locs]
        sigfigs.sort()
        self.format = '%1.' + str(sigfigs[-1]) + 'f'
        if self._usetex or self._useMathText: self.format = '$%s$'%self.format

    def pprint_val(self, x):
        xp = (x-self.offset)/10**self.orderOfMagnitude
        if np.absolute(xp) < 1e-8: xp = 0
        s = self.format % xp
        if self.thousandsSep: # add thousands-separating characters
            if s.count('.'): # it's got a decimal in there
                # use the regexp for floats:
                s = re.sub(r'(?<=\d)(?=(\d\d\d)+\.)', self.thousandsSep, s)
            else: # it's an int
                # use the regexp for ints:
                s = re.sub(r'(?<=\d)(?=(\d\d\d)+$)', self.thousandsSep, s)
        return s


class NeuropyAutoLocator(mpl.ticker.MaxNLocator):
    """A tick autolocator that generates more ticks than the standard mpl autolocator"""
    def __init__(self):
        # standard autolocator:
        #mpl.ticker.MaxNLocator.__init__(self, nbins=9, steps=[1, 2, 5, 10])
        mpl.ticker.MaxNLocator.__init__(self) # use MaxNLocator's defaults instead


def getargstr(obj):
    """Returns object's argument list as a string. Stolen from wx.py package?"""
    import inspect
    argstr = apply(inspect.formatargspec, inspect.getargspec(obj))
    if inspect.isfunction(obj):
        pass
    elif inspect.ismethod(obj):
        # stolen from wx.py.introspect.getCallTip:
        temp = argstr.split(',')
        if len(temp) == 1:  # No other arguments.
            argstr = '()'
        elif temp[0][:2] == '(*': # first param is like *args, not self
            pass
        else:  # Drop the first argument.
            argstr = '(' + ','.join(temp[1:]).lstrip()
    else:
        argstr = '()'
    return argstr
'''
def frame(**kwargs):
    """Returns a CanvasFrame object"""
    frame = CanvasFrame(**kwargs)
    frame.Show(True)
    return frame
frame.__doc__ += '\n' + CanvasFrame.__doc__
frame.__doc__ += '\n\n**kwargs:\n' + getargstr(CanvasFrame.__init__)

def barefigure(*args, **kwargs):
    """Creates a bare figure with no toolbar or statusbar"""
    figure(*args, **kwargs)
    gcfm().frame.GetStatusBar().Hide()
    gcfm().frame.GetToolBar().Hide()
barefigure.__doc__ += '\n' + figure.__doc__
'''
def lastcmd():
    """Returns a string containing the last command entered by the user in the
    IPython shell"""
    ip = get_ipython()
    return ip._last_input_line

def innerclass(cls):
    """Class decorator for making a class behave as a Java (non-static) inner
    class.

    Each instance of the decorated class is associated with an instance of its
    enclosing class. The outer instance is referenced implicitly when an
    attribute lookup fails in the inner object's namespace. It can also be
    referenced explicitly through the property '__outer__' of the inner
    instance.

    Title: Implementing Java inner classes using descriptors
    Submitter: George Sakkis - gsakkis at rutgers.edu
    Last Updated: 2005/07/08
    Version no: 1.1
    Category: OOP
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/409366
    """
    if hasattr(cls, '__outer__'):
        raise TypeError('Cannot set attribute "__outer__" in inner class')
    class InnerDescriptor(object):
        def __get__(self, outer, outercls):
            if outer is None:
                raise AttributeError('An enclosing instance that contains '
                           '%s.%s is required' % (cls.__name__, cls.__name__))
            clsdict = cls.__dict__.copy()
            # explicit read-only reference to the outer instance
            clsdict['__outer__'] = property(lambda self: outer)
            # implicit lookup in the outer instance
            clsdict['__getattr__'] = lambda self,attr: getattr(outer,attr)
            def __setattr__(this, attr, value):
                # setting an attribute in the inner instance sets the
                # respective attribute in the outer instance if and only if
                # the attribute is already defined in the outer instance
                if hasattr(outer, attr): setattr(outer,attr,value)
                else: super(this.__class__,this).__setattr__(attr,value)
            clsdict['__setattr__'] = __setattr__
            return type(cls.__name__, cls.__bases__, clsdict)
    return InnerDescriptor()

def intround(n):
    """Round to the nearest integer, return an integer. Works on arrays.
    Saves on parentheses, nothing more"""
    if iterable(n): # it's a sequence, return as an int64 array
        return np.int64(np.round(n))
    else: # it's a scalar, return as normal Python int
        return int(round(n))

def roundto(val, nearest):
    """Round val to nearest nearest, always rounding away from 0"""
    if val >= 0:
        return np.ceil(val / nearest) * nearest
    else:
        return np.floor(val / nearest) * nearest

def sigfig(x, n=1):
    """Return x rounded to n significant figures. Modified from
    http://code.activestate.com/lists/python-tutor/70739/"""
    return round(x, int(n - np.ceil(np.log10(abs(x)))))

def e10(x):
    """Return exponent of power of 10. Example: if x is 3.5e-10, return -10"""
    if x == 0:
        return 0
    return int(np.floor(np.log10(abs(x))))

def ceilsigfig(x, n=1):
    """Return x rounded up to n significant figures. This is useful when wanting to
    write p < value for significance tests. Works for negative numbers too.
    Example:  2.1e-10 -->  3e-10
    Example:  2.9e-10 -->  3e-10
    Example: -2.1e-10 --> -2e-10"""
    sfx = sigfig(x, n)
    if sfx < x: # it was rounded down
        sfx = sfx + 10**(e10(x)) # add one at the same decimal place
        # filter through sigfig() again to try and fix any float inaccuracy,
        # example: -2.6e-10 --> -3e-10 --> -1.9999999999999998e-10 --> -2e-10
        sfx = sigfig(sfx, n)
    return sfx

def floorsigfig(x, n=1):
    """Return x rounded down to n significant figures. This is useful when wanting to
    write p > value for significance tests. Works for negative numbers too.
    Example:  2.1e-10 -->  2e-10
    Example:  2.9e-10 -->  2e-10
    Example: -2.1e-10 --> -3e-10"""
    sfx = sigfig(x, n)
    if sfx > x: # it was rounded up
        sfx = sfx - 10**(e10(x)) # subtract one at the same decimal place
        # filter through sigfig() again to try and fix any float inaccuracy,
        # example: -2.4e-10 --> -2e-10 --> -2.9999999999999998e-10 --> -3e-10
        sfx = sigfig(sfx, n)
    return sfx

def pad0s(val, ndigits=2):
    """Returns a string rep of val, padded with enough leading 0s
    to give you a string rep with ndigits in it"""
    val = str(int(val))
    nzerostoadd = ndigits - len(val) # -ve values add no zeros to val
    val = '0'*nzerostoadd + val
    return val

def txtdin2binarydin(fin, fout):
    """Converts a csv text .din file to an int64 binary .din file"""
    fi = file(fin, 'r') # open the din file for reading in text mode
    fo = file(fout, 'wb') # for writing in binary mode
    for line in fi:
        line = line.split(',')
        '''
        # for old NVS display, converts from NVS condition numbers (which increment with
        # repeats) to dimstim sweepis (which don't)
        nruns = 18
        line[1] = int(line[1]) % nruns
        '''
        # write both values out as a C long longs, using the system's native ('@') byte order:
        fo.write( struct.pack('@qq', int(line[0]), int(line[1])) )
    fi.close()
    fo.close()
    print('Converted ascii din: %r to binary din: %r' % (fin, fout))

def convertalltxtdin2binarydin(path=None):
    """Converts all text .csv din files in path (or cwd) to 64 bit binary .din files of the
    same name"""
    if path == None:
        path = os.getcwd()

    listing = os.listdir(path)
    dinfnames = []

    for fname in listing:
        if fname.endswith('.csv'):
            # text din filenames without the .csv extension
            dinfnames.append(fname[:-len('.csv')])

    for dinfname in dinfnames:
        fin = os.path.join(path, dinfname) + '.csv'
        fout = os.path.join(path, dinfname) + '.din'
        #os.rename(fout, fin) # rename the csv .din file to .din.txt extension
        txtdin2binarydin(fin, fout) # convert the text .csv file to binary .din file
        #os.remove(fin) # delete the .din.txt file

def renameSpikeFiles(path, newname):
    """Renames all .spk files in path to newname, retaining their '_t##.spk' ending"""
    for fname in os.listdir(path):
        if fname.endswith('.spk'):
            i = fname.find('_t')
            if i != -1:
                newfname = newname+fname[i::]
                print(newfname)
                os.rename(os.path.join(path, fname), os.path.join(path, newfname))

def csv2binary(fin, multiplier=1e6, skipfirstline=True):
    """Converts spike data in a csv file, with cells in the columns and times down the rows,
    into int64 binary files, one for each neuron. Takes csv values and multiplies them by
    multiplier before saving"""
    fin = os.path.normpath(fin)
    fi = file(fin, 'r') # open csv file for reading in text mode
    print('Exporting %s to:' % fi.name)
    firstline = fi.next()
    nneurons = len(firstline.split(','))
    if not skipfirstline: # ie first line isn't just column headers
        fi.seek(0)
    data = [] # nested list, one entry per neuron
    for ni in range(nneurons):
        data.append([]) # init each neuron's list
    for line in fi:
        line = line.replace('\n', '') # strip the newline character
        line = line.split(',')
        for ni, strval in enumerate(line): # going horizontally across the line
            try:
                data[ni].append(intround(float(strval)*multiplier))
            except ValueError: # strval is empty string
                pass
    fi.close()
    #return data
    path = os.path.splitext(fi.name)[0] # extensionless path + filename
    try:
        os.mkdir(path) # make a dir with that name
    except OSError: # dir already exists
        pass
    # just the extensionless filename, replace spaces with underscores:
    tail = os.path.split(path)[-1].replace(' ', '_')
    for ni, neuron in enumerate(data):
        fname = (os.path.join(path, tail) + '_t' +
                 pad0s(ni, ndigits=len(str(nneurons))) + '.spk')
        fo = file(fname, 'wb') # for writing in binary mode
        print(fo.name)
        for spiketime in neuron:
            # write each spiketime to the file, there should be a more streamlined way
            # to do this. Write the value out as a C long long, using the system's native
            # ('@') byte order
            fo.write( struct.pack('@q', spiketime) )
        fo.close()
'''
def warn(msg, level=2, exit_val=1):
    """Standard warning printer. Gives formatting consistency. Stolen from IPython.genutils"""
    if level > 0:
        header = ['', '', 'WARNING: ', 'ERROR: ', 'FATAL ERROR: ']
        print >> sys.stderr, '%s%s' % (header[level],msg)
        if level == 4:
            print >> sys.stderr, 'Exiting.\n'
            sys.exit(exit_val)

def warn(msg):
    import warnings
    warnings.warn(msg, category=RuntimeWarning, stacklevel=2)

def unique(seq):
    """Return unique items from a 1-dimensional sequence. Stolen from numpy.unique().
    Dictionary setting is quite fast"""
    result = {}
    for item in seq:
        result[item] = None
    return result.keys()

def unique(objlist):
    """Returns the input list minus any repeated objects it may have had.
    Also defined in dimstim"""
    return list(set(objlist)) # this requires Python >= 2.4

def unique(objlist):
    """Does in-place removal of non-unique objects in a list of objects"""
    for (i,obj1) in enumerate(objlist):
        for (j,obj2) in enumerate(objlist):
            if i != j and obj1 == obj2:
                del objlist[j]
'''
def iterable(x):
    """Check if the input is iterable, stolen from numpy.iterable()"""
    try:
        iter(x)
        return True
    except:
        return False

def toiter(x):
    """Convert to iterable. If input is iterable, returns it. Otherwise returns it in a list.
    Useful when you want to iterate over an object (like in a for loop),
    and you don't want to have to do type checking or handle exceptions
    when the object isn't a sequence"""
    if iterable(x):
        return x
    else:
        return [x]

def tolist(x):
    """Convert to list. If input is a dict, return its values. If it's an array, convert it to
    a list (NOTE: this new feature, as of 2014-04-07, hasn't been widely tested in all places
    where tolist() is called). If it's already a list, return it. Otherwise, return input in a
    list. A little different from np.atleast_1d."""
    tx = type(x)
    if tx == dict:
        return list(x.values())
    elif tx == np.ndarray:
        return list(x)
    elif tx == list:
        return x
    else:
        return [x] # stick it in a list

def to2d(arr):
    """Convert a 1D array to a 2D array with just a singleton row. If arr is already
    2D, just return it. If it's anything more than 2D, raise an error"""
    nd = arr.ndim
    assert nd in [1, 2], 'array rank > 2'
    if nd == 1:
        arr = arr.reshape(1, -1)
    return arr

def joinpath(pathlist):
    """Unlike os.path.join(), take a list of path segments, return them joined in a string
    with local separators"""
    path = ''
    for segment in pathlist:
        path = os.path.join(path, segment)
    return path

def dist(a, b):
    """Return the Euclidean distance between two N-dimensional coordinates"""
    a = np.asarray(a)
    b = np.asarray(b)
    return np.sqrt(((a-b)**2).sum())

def approx(a, b, rtol=1.e-14, atol=1.e-14):
    """Return a boolean array describing which components of a and b are equal
    subject to given tolerances. The relative error rtol must be positive and << 1.0
    The absolute error atol comes into play for those elements of y that are very
    small or zero; it says how small x must be also. Copied and modified from
    numpy.allclose()"""
    x = np.array(a, copy=False)
    y = np.array(b, copy=False)
    #print(x.shape)
    #print(y.shape)
    return np.less(np.absolute(x-y), atol + rtol * np.absolute(y))

def pmf(a, bins=10, range=None, weights=None):
    """Return probability mass function of a, where sum of all bins is 1. If bins is
    iterable, make sure to include the rightmost bin edge. Unlike np.histogram, returned
    bins exclude the rightmost bin edge"""
    n, bins = np.histogram(a, bins=bins, range=range, weights=weights, density=False)
    n = n / float(sum(n)) # normalize by sum of bins to get pmf
    return n, bins[:-1]

def pmf2d(a, bins=10, range=None, weights=None):
    """Return 2D probability mass function of a, where sum of all bins is 1. If bins is
    iterable, make sure to include the rightmost bin edge. Unlike np.histogram, returned
    bins exclude the rightmost bin edge"""
    H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=range, normed=False,
                                       weights=weights)
    H = H / float(sum(H)) # normalize by sum of bins to get pmf
    return H, xedges[:-1], yedges[:-1]

def sah(t, y, ts, keep=False):
    """Resample using sample and hold. Returns resampled values at ts given the original
    points (t,y) such that the resampled values are just the most recent value in y (think
    of a staircase with non-uniform steps). Assumes that t is sorted. t and ts arrays should
    be of the same data type. Contributed by Robert Kern."""
    # find where ts falls in t, dec so you get indices that point to the most
    # recent value in y:
    i = np.searchsorted(t, ts) - 1
    # handle the cases where ts is smaller than the first point.
    '''this has an issue of not keeping the original data point where ts == t'''
    i = np.where(i < 0, 0, i)

    ### NOTE: can probably get around having to do this by using searchsorted's
    ### 'side' keyword

    if keep:
        # The following ensures that the original data point is kept when ts == t,
        # doesn't really work if the shortest ISI is less than tres in ts.
        # find changes in i, nonzero() method returns a tuple, pick the result for the
        # first dim with [0] index
        di = np.diff(i).nonzero()[0]
        # check at those change indices if t ~= ts (ignoring potential floating point
        # representational inaccuracies). If so, inc i at that point so you keep y at that
        # point.
        si = approx(t[1::], ts[di])
        #print(i)
        i[di[si]] += 1
        #print(i)
    return y[i]

def corrcoef(x, y):
    """Returns the correlation coefficient of signals x and y. This just uses np.corrcoef(),
    but converts to floats first, cuz np.corrcoef() seems to have issues with integer signals,
    especially those with zeros in them."""
    #assert len(x) == len(y), 'arrays need to be of equal length'
    x = np.float64(x)
    y = np.float64(y)
    # pick one of the 2 entries in the correlation coefficient matrix, on the -ve diagonal
    # (er, the one that goes from bottom left to top right, that's what I mean):
    return np.corrcoef(x, y)[0, 1]
    # this works just fine as well, easier to understand too:
    #return ((x * y).mean() - x.mean() * y.mean()) / (x.std() * y.std())

def bin(i, minbits=8):
    """Return a string with the binary representation of an integer, or sequence of integers.
    If necessary, will append leading zeros if result is less than minbits long.
    Uses np.binary_repr()"""
    ints = toiter(i) # ensure it's iterable
    sints = []
    for i in ints:
        s = np.binary_repr(i)
        nzerostoadd = minbits - len(s) # OK if this is -ve
        s = '0'*nzerostoadd + s # add enough leading zeros to get requested minbits
        sints.append(s)
    if len(sints) == 1:
        sints = sints[0] # pull it out of the list
    return sints

def binslow(i, minbits=8):
    """Return a string with the binary representation of an integer. If necessary, will
    append leading zeros if result is less than minbits long. Seems like np.binary_repr() is
    a somewhat faster alternative. First 2 lines stolen from Andrew Gaul <andrew@gaul.org>
    off the web"""
    l = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111',
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    s = ''.join(map(lambda x, l=l: l[int(x, 16)], hex(i)[2:]))
    s = s.lstrip('0') # strip s of leading zeros
    nzerostoadd = minbits - len(s)
    s = '0'*nzerostoadd + s # add enough leading zeros to get requested minbits
    return s

# an alternative would be to use int('10110', base=2) for each column, probably slower though
def binarray2int(bin):
    """Takes a 2D binary array (only 1s and 0s, with rows LSB to MSB from top to bottom)
    and returns the base 10 integer representations of the columns"""
    #assert type(bin) == type(np.array)
    bin = to2d(bin) # ensure it's 2D. If it's 1D, force it into having a singleton row
    nbits = bin.shape[0] # length of the first dimension, ie the number of rows
    multiplier = []
    for i in range(nbits):
        multiplier.append(2**i)
    # convert from list and transpose to a column vector (have to make it 2D to transpose):
    multiplier = np.array(multiplier, ndmin=2).transpose()
    #print(multiplier)
    x = bin*multiplier
    #print(x)
    # sum over the first dimension (the rows), that way, you're left with only columns in
    # a row vector:
    return x.sum(axis=0)

def getbinarytable(nbits=8):
    """Generate a 2D binary table containing all possible words for nbits, with bits in the
    rows and words in the columns (LSB to MSB from top to bottom)"""
    rowlength = 2**nbits
    '''
    x = np.zeros((nbits, 2**nbits)) # init an array
    for bit in range(nbits):
        pattern = [0]*2**bit
        pattern.extend([1]*2**bit)
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits /
                                             # 2**(bit+1) == 2**(nbits-bit-1)
        row = pattern*npatterns
        x[bit]=row
    return x
    '''
    '''
    x = np.zeros((nbits, 2**nbits), dtype=np.int8) # init an array
    for bit in range(nbits): # one row at a time
        pattern = np.array(0, dtype=np.int8).repeat(2**bit)
        pattern = np.concatenate((pattern, np.array(1, dtype=np.int8).repeat(2**bit)))
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits /
                                             # 2**(bit+1) == 2**(nbits-bit-1)
        row = np.tile(pattern, [1, npatterns])
        x[bit::,::] = row
    return x
    '''
    # this seems to be the fastest method:
    x = []
    for bit in range(nbits): # one row at a time
        pattern = np.array(0, dtype=np.int8).repeat(2**bit)
        pattern = np.concatenate((pattern, np.array(1, dtype=np.int8).repeat(2**bit)))
        # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        npatterns = rowlength / len(pattern)
        row = np.tile(pattern, [1, npatterns])
        x.append(row)
    return np.concatenate(x)

def enlarge(a, x=2, y=None):
    """Enlarges 2D image array a using simple pixel repetition in both dimensions.
    Enlarges by factor x horizontally and factor y vertically.
    If y is left as None, uses factor x for both dimensions."""
    a = np.asarray(a)
    assert a.ndim == 2
    if y == None:
        y = x
    for factor in (x, y):
        assert factor.__class__ == int
        assert factor > 0
    return a.repeat(y, axis=0).repeat(x, axis=1)

def charfind(string, char):
    """Finds char in string, returns matching indices. There's gotta be a built-in way to do
    this somewhere..."""
    assert len(char) == 1
    i = []
    # maybe more efficient to use .find() method on successively smaller slices of string
    for si, s in enumerate(string):
        if s == char:
            i.append(si)
    return i
'''
def shuffle(x):
    """Takes an input list x and returns a shuffled (without replacement) copy. Its only
    benefit over and above random.sample() is that you don't have to pass a second argument
    len(x) every time you use it. In NumPy, it's better (and faster) to use
    np.random.shuffle()"""
    return random.sample(x, len(x))
'''
def shuffle(seq):
    """Takes a sequence and returns a shuffled (without replacement) copy. Its only benefit
    over np.random.shuffle is that it returns a copy instead of shuffling in-place"""
    result = copy(seq)
    np.random.shuffle(result) # shuffles in-place, doesn't convert to an array
    return result
'''
def randomize(x):
    """Takes an input list x and returns a randomized (with replacement) output list of
    the same length, sampled from the input sequence"""
    y = [] # init output list
    for i in range(0, len(x)):
        y.append(random.choice(x))
    return y
'''
def randomize(seq):
    """Returns a randomized (with replacement) output sequence sampled from
    (and of the same length as) the input sequence"""
    n = len(seq)
    i = np.random.randint(n, size=n) # returns random ints from 0 to len(seq)-1
    if seq.__class__ == np.ndarray:
        return np.asarray(seq)[i] # use i as random indices into seq, return as an array
    else:
        return list(np.asarray(seq)[i]) # return as a list

def randsign(size=1):
    """Return random array of -1 and 1 integers"""
    rand = np.random.random(size)
    signs = np.ones(size, dtype=np.int)
    signs[rand < 0.5] = -1
    return signs

def fact(n):
    """Factorial!"""
    assert type(n) == int
    assert n >= 0
    if n == 0:
        n = 1 # 0! == 1!
    #return int(np.prod(np.arange(n, 1, -1))) # this easily overflows given int64
    # Python int math doesn't overflow:
    result = n
    for i in range(1, n):
        result *= i
    return result

def nPr(n, r):
    """n Pick r"""
    assert n >= r
    return fact(n) // fact(n-r)

def nCr(n, r):
    """n Choose r"""
    assert n >= r
    return fact(n) // (fact(n-r) * fact(r))

ncr = nCr # convenience f'ns
npr = nPr

def combgen(objects, r=2, i=None, level=0):
    """Generator that yields, without replacement, all length r possible combinations of
    objects from a length n sequence. Eg, if objects=[0,1,2] and r=2, this yields [0,1],
    [0,2], and [1,2], one at a time. A recursive generator is used in order to create the
    necessary r number of nested for loops. This is cool (my first generator!), but deep
    recursion is slow"""
    objects = np.asarray(objects)
    assert r <= len(objects)
    try: # recursive case
        if i == None:
            i = [0]*r # stores all the current index values for all r nested for loops
        if level == 0: # handles special case for starting index of top level for loop
            starti = 0
        else:
            # start this level's loop index at one greater than the previous level's
            # current loop index:
            starti = i[level-1] + 1
        # not too sure why this is n+1, but it works:
        for i[level] in range(starti, len(objects)+1):
            # iterate over next level's generator:
            for comb in combgen(objects, r=r, i=i, level=level+1):
                # yield whatever the next level (level+1) yields, pass it on up to the
                # previous level (level-1)
                yield comb 
    except IndexError:
        # base case, we're at the deepest recursion level (innermost for loop). IndexError
        # comes from i[level] being out of range:
        #if len(i) == 1:
        #    yield objects[i[0]] # no need to yield them in a list
        #else:
            # use the current index state for all levels to yield a combination of objects:
            yield objects[i]

def combs(objects, r=2):
    """Returns all nCr possible combinations of items in objects, in a 1D array of arrays.
    Generates code with the right number of nested for loops, faster than combgen()"""
    objects = np.asarray(objects)
    dtype = objects.dtype
    n = len(objects)
    assert r <= n
    i = np.asarray([0]*r)
    # stores all combinations, will be a 1D array of arrays:
    combs = np.empty(nCr(n, r), dtype=np.object)
    combi = -1

    code = ''
    tabs = ''
    code += tabs+'for i[0] in range(0, n):\n' # this is the outermost for loop
    tabs += '\t'
    for level in range(1, r): # here come the inner nested for loops...
        code += tabs+'for i['+str(level)+'] in range(i['+str(level-1)+']+1, n):\n'
        tabs += '\t'

    # here's the innermost part of the nested for loops
    code += tabs + 'combi += 1\n'
    code += tabs + 'combs[combi] = objects[i]\n'
    #print(code)

    exec(code) # run the generated code
    return combs
    '''
    # example of what the generated code looks like for r==3:
    for i[0] in range(0, n):
        for i[1] in range(i[0]+1, n):
            for i[2] in range(i[1]+1, n):
                combi += 1
                combs[combi] = objects[i]
    '''

def argcombs(objects, r=2):
    """Returns all nCr possible combinations of indices into objects.
    You'd think this would be faster than combs(), but it doesn't seem to be"""
    n = len(objects)
    assert n < 2**8 # this way, we can use uint8's instead of int32's to save memory
    assert r <= n
    i = np.asarray([0]*r)
    argcombs = np.zeros((nCr(n, r), r), dtype=np.uint8)
    combi = -1

    code = ''
    tabs = ''
    code += tabs+'for i[0] in range(0, n):\n' # this is the outermost for loop
    tabs += '\t'
    for level in range(1, r): # here come the inner nested for loops...
        code += tabs+'for i['+str(level)+'] in range(i['+str(level-1)+']+1, n):\n'
        tabs += '\t'

    # here's the innermost part of the nested for loops
    code += tabs + 'combi += 1\n'
    code += tabs + 'argcombs[combi, :] = i\n'
    #print(code)

    exec(code) # run the generated code
    return argcombs
    '''
    # example of what the generated code looks like for r==3:
    for i[0] in range(0, n):
        for i[1] in range(i[0]+1, n):
            for i[2] in range(i[1]+1, n):
                combi += 1
                argcombs[combi, :] = i
    '''

def nCrsamples(objects, r, nsamples=None):
    """Returns a list of nsamples unique samples, each of length r, sampled from objects"""
    nobjects = len(objects)
    maxnsamples = nCr(nobjects, r)
    if nsamples == None:
        nsamples = maxnsamples # return all possible combinations
    if nsamples > maxnsamples:
        # make sure we're not being asked for more than the maximum possible number of
        # unique samples
        raise ValueError('requested unique nsamples (%d) is larger than nobjects choose '
                         'r (%d C %d == %d)' % (nsamples, nobjects, r, maxnsamples))
    # Don't use an exhaustive table, because generating it and then sampling it almost
    # always takes longer than just picking combs at random and making sure they're unique
    '''
    if maxnsamples < something:
        # generate a table of all possible combinations, and then just pick nsamples from
        # it without replacement
        table = combs(objects, r)
        samples = random.sample(table, nsamples)
    '''
    if r == 1: # we're just choosing one item from objects at a time
        samples = random.sample(objects, nsamples)
    else: # sample at random and make sure each sample is unique:
        samples = []
        for samplei in range(nsamples):
            sample = sortedsample(objects, r)
            while sample in samples: # repeat until a unique one is found
                sample = sortedsample(objects, r)
            samples.append(sample)
    return samples

def sortedsample(objects, r):
    """Randomly sample r things from objects, sorting the result to remove permutations"""
    sample = random.sample(objects, r)
    sample.sort()
    return sample

'''
# this f'n isn't really needed, just use objlist.sort(key=lambda obj: obj.attrib)
def sortby(objs, attrib, cmp=None, reverse=False):
    """Returns objects list sorted according to the specified object attribute.
    attrib should be passed as a string"""
    # sort in-place:
    objs.sort(key=lambda obj: obj.__getattribute__(attrib), cmp=cmp, reverse=reverse)
    return objs
'''
def intersect1d(arrays, assume_unique=False):
    """Find the intersection of any number of 1D arrays.
    Return the sorted, unique values that are in all of the input arrays.
    Adapted from numpy.lib.arraysetops.intersect1d"""
    N = len(arrays)
    if N == 0:
        return np.asarray(arrays)
    arrays = list(arrays) # allow assignment
    if not assume_unique:
        for i, arr in enumerate(arrays):
            arrays[i] = np.unique(arr)
    aux = np.concatenate(arrays) # one long 1D array
    aux.sort() # sorted
    if N == 1:
        return aux
    shift = N-1
    return aux[aux[shift:] == aux[:-shift]]

def mean_accum(data):
    """Takes mean by accumulating over 0th axis in data,
    much faster than np.mean() because it avoids making any copies of the data
    Suggested by Tim Hochberg"""
    result = np.zeros(data[0].shape, np.float64) # init output array
    for dataslice in data:
        # this loop isn't so bad since the massive add step is the limiting factor
        result += dataslice
    result /= len(data)
    return result

def mean_accum2(data, indices):
    """A variant of mean_accum(), where you provide all the data and the indices into it
    to average over. This was Tim Hochberg's version"""
    result = np.zeros(data[0].shape, np.float64)
    for i in indices:
        result += data[i]
    result /= len(indices)
    return result

def normalize_range(a):
    """Normalize a such that all of its values span the interval [0, 1]"""
    a = a - a.min()
    a = a / a.max()
    return a

def normalize(p):
    """Normalize distribution p such that sum(p) == 1. Return zeros if sum(p) == 0"""
    p = np.asarray(p)
    if p.sum() == 0:
        return np.zeros(p.shape) # just return zeros
    else:
        return p / float(p.sum()) # return it normalized

def ensurenormed(p, atol=1e-8):
    """Ensures p is normalized. Returns p unchanged if it's already normalized,
    otherwise, prints a warning and returns it normalized. atol is how close to 1.0
    p.sum() needs to be"""
    p = np.asarray(p)
    psum = p.sum()
    if not approx(psum, 1.0, atol=atol): # make sure the probs sum to 1
        print("ps don't sum to 1, they sum to %f instead, normalizing for you" % psum)
        p = p / float(psum) # in case p are integers, don't divide in place
    return p

def logn(x, base=10):
    """Performs log of x with specified base"""
    return np.log(x) / np.log(base)

def log_no_sing(x, subval=0.0, base=np.e):
    """Performs log on array x, ignoring any zeros in x to avoid singularities,
    and returning subval in their place in the result"""
    x = np.asarray(x)
    singi = x==0 # find the singularities
    x[singi] = 1 # replace 'em with 1s, or anything else that's safe to take the log of
    result = logn(x, base=base) # now it's safe to take the log
    # substitute the result where the singularities were with the substitution value:
    result[singi] = subval
    return result

def log10_no_sing(x, subval=0.0):
    """Performs log10 on x, ignoring singularities"""
    return log_no_sing(x, subval=subval, base=10)

def log2_no_sing(x, subval=0.0):
    """Performs log2 on x, ignoring singularities"""
    return log_no_sing(x, subval=subval, base=2)

def entropy(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob
    values in p"""
    p = ensurenormed(p)
    return -(p * np.log2(p)).sum()

def entropy_no_sing(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob values
    in p. Ignore singularities in p (assumes their contribution to entropy is zero)"""
    p = ensurenormed(p)
    return -(p * log2_no_sing(p, subval=0.0)).sum()

def MI(XY):
    """Given the joint PDF of two variables, return the mutual information (in bits)
    between the two.
    I = sum_X sum_Y P(x, y) * log2( P(x, y) / (P(x) * P(y)) )
    where P(x) and P(y) are the marginal distributions taken from the joint
    Is this slow? Needs optimization? Already implemented in scipy?"""
    XY = np.asarray(XY)
    assert XY.ndim == 2
    XY = ensurenormed(XY)
    # calculate the marginal probability distributions for X and Y from the joint
    X = XY.sum(axis=1) # sum over the rows of the joint, get a vector nrows long
    Y = XY.sum(axis=0) # sum over the cols of the joint, get a vector ncols long
    I = 0.0
    for xi, x in enumerate(X):
        for yi, y in enumerate(Y):
            if XY[xi, yi] == 0 or (x * y) == 0: # avoid singularities
                pass # just skip it, assume info contributed is 0 (?????????????????)
            else:
                I += XY[xi, yi] * np.log2( XY[xi, yi] / (x * y) )
    return I

def MIbinarrays(Nbinarray=None, Mbinarray=None, verbose=False):
    """Calculates information that N cells provide about M cells (ie,
    their mutual information), as a fraction of the M cells' marginal entropy.
    Takes cell activities as binary arrays (on or off), with
    cells in the rows and time bins in the columns."""
    Nbinarray = to2d(Nbinarray) # make it 2D if it's 1D
    N = len(Nbinarray) # gets the number of rows
    Nintcodes = binarray2int(Nbinarray)
    Mbinarray = to2d(Mbinarray) # make it 2D if it's 1D
    M = len(Mbinarray) # gets the number of rows
    Mintcodes = binarray2int(Mbinarray)
    # build up joint pdf of all the possible N words, and the two possible N+1th values
    # (0 and 1)
    # values 0 to 2**N - 1, plus 2**N which is needed as the rightmost bin edge for
    # histogram2d:
    xedges = np.arange(2**N+1) # include rightmost edge
    yedges = np.arange(2**M+1) # include rightmost edge
    bins = [xedges, yedges]
    # generate joint pdf, *edgesout exclude rightmost edge:
    jpdf, xedgesout, yedgesout = pmf2d(Nintcodes, Mintcodes, bins)
    #print('jpdf\n', jpdf.__repr__())
    #print('jpdf.sum()', jpdf.sum())
    assert (np.float64(xedges)[:-1] == xedgesout).all() # sanity check
    assert (np.float64(yedges)[:-1] == yedgesout).all()
    # pdf of N cells
    #Npdf, Nedges = pmf(Nintcodes, bins=range(2**N + 1))
    #print('first 100 Npdf\n', Npdf[:100].__repr__())
    # pdf of M cells
    #Mpdf, Medges = pmf(Mintcodes, bins=range(2**M + 1))
    #print('first 100 Mpdf\n', Mpdf[:100].__repr__())
    marginalMpdf = jpdf.sum(axis=0)
    # make sure what you get from the joint is what you get when just building up the
    # pdf straight up on its own:
    #assert approx(Mpdf, marginalMpdf).all()
    I = MI(jpdf)
    # mutual info as fraction of entropy in M group of cells:
    IdivS = I / entropy(marginalMpdf)
    if verbose:
        print('nids', nids)
        print('mids', mids)
        #print('Mpdf', Mpdf)
        #print('entropy(Mpdf)', entropy(Mpdf))
        print('marginal Mpdf', marginalMpdf)
        print('entropy(marginal Mpdf)', entropy(marginalMpdf))
        print('I', I)
        print('I/entropy', IdivS)
    if not 0.0 <= IdivS <= 1.0+1e-10:
        import pdb; pdb.set_trace()
        print('IdivS is out of range')
        print('IdivS is %.16f' % IdivS)

    return dictattr(I=I, IdivS=IdivS)

def DKL(p, q):
    """Kullback-Leibler divergence from true probability distribution p
    to arbitrary distribution q"""
    assert len(p) == len(q)
    p = ensurenormed(p)
    q = ensurenormed(q)
    # avoid singularities:
    return sum([ pi * np.log2(pi/float(qi)) for pi, qi in zip(p, q) if pi != 0 and qi != 0 ] ) 

def DJS(p, q):
    """Jensen-Shannon divergence, a symmetric measure of divergence between
    distributions p and q"""
    p = np.asarray(p) # required for adding p and q
    q = np.asarray(q)
    m = 1 / 2.0 * (p + q)
    return 1 / 2.0 * ( DKL(p, m) + DKL(q, m) )

def lstrip(s, strip):
    """What I think str.lstrip should really do"""
    if s.startswith(strip):
        return s[len(strip):] # strip it
    else:
        return s

def rstrip(s, strip):
    """What I think str.rstrip should really do"""
    if s.endswith(strip):
        return s[:-len(strip)] # strip it
    else:
        return s

def strip(s, strip):
    """What I think str.strip should really do"""
    return rstrip(lstrip(s, strip), strip)

def lrstrip(s, lstr, rstr):
    """Strip lstr from start of s and rstr from end of s"""
    return rstrip(lstrip(s, lstr), rstr)

def pathdecomp(path):
    """Decompose (fully split) all components of path into a list of strings
    If the first string is empty, that means the second string was relative to
    filesystem root"""
    return path.split(os.path.sep)

def eof(f):
    """Return whether file pointer is a end of file"""
    orig = f.tell()
    f.seek(0, 2) # seek 0 bytes from end
    return f.tell() == orig

def td2usec(td):
    """Convert datetime.timedelta to microseconds"""
    sec = td.total_seconds() # float
    usec = intround(sec * 1000000) # round to nearest us
    return usec

def issorted(x):
    """Check if x is sorted"""
    try:
        if x.dtype.kind == 'u':
            # x is unsigned int array, risk of int underflow in np.diff
            x = np.int64(x)
    except AttributeError:
        pass # no dtype, not an array
    return (np.diff(x) >= 0).all() # is difference between consecutive entries >= 0?
    # or, you could compare the array to an explicitly sorted version of itself,
    # and see if they're identical

def inverse_uquadratic_cdf(y, a=0, b=1):
    assert b > a
    alpha = 12 / ((b - a)**3)
    beta = (b + a) / 2
    print(y * 3 / alpha - (beta - a)**3)
    return cbrt(y * 3 / alpha - (beta - a)**3) + beta

def sample_uquadratic(a=0, b=1, size=None):
    """Randomly sample the U-quadratic distribution. Good for modelling
    bimodal distributions. a and b specify upper and lower bounds.
    See:
    http://en.wikipedia.org/wiki/UQuadratic_distribution
    http://distributome.org/js/exp/UQuadraticExperiment.html
    """
    assert b > a
    x = np.random.random(size=size) # sample uniform distrib
    x = (b - a) * x + a # scale so that min(x) == a and max(x) == b
    return inverse_uquadratic_cdf(x, a, b)

def split_tranges(tranges, width, tres):
    """Split up tranges into lots of smaller (typically overlapping) tranges, with width and
    tres. Usually, tres < width, but this also works for width < tres.
    Test with:

    print(split_tranges([(0,100)], 1, 10))
    print(split_tranges([(0,100)], 10, 1))
    print(split_tranges([(0,100)], 10, 10))
    print(split_tranges([(0,100)], 10, 8))
    print(split_tranges([(0,100)], 3, 10))
    print(split_tranges([(0,100)], 10, 3))
    print(split_tranges([(0,100)], 3, 8))
    print(split_tranges([(0,100)], 8, 3))
    """
    newtranges = []
    for trange in tranges:
        t0, t1 = trange
        assert width < (t1 - t0)
        # calculate left and right edges of subtranges that fall within trange:
        # This is tricky: find maximum left edge such that the corresponding maximum right
        # edge goes as close as possible to t1 without exceeding it:
        tend = (t1-width+tres) // tres*tres # there might be a nicer way, but this works
        ledges = np.arange(t0, tend, tres)
        redges = ledges + width
        subtranges = [ (le, re) for le, re in zip(ledges, redges) ]
        newtranges.append(subtranges)
    return np.vstack(newtranges)

def laminarity(ypos, trackabsname):
    """Return boolean arrays indicating whether depths ypos are superficial, middle,
    or deep layer (or none of the above)"""
    uns = get_ipython().user_ns
    try:
        (sup0, sup1), (mid0, mid1), (deep0, deep1) = uns['LAYERS'][trackabsname]
    except KeyError: # trackabsname doesn't exist as key in LAYERS global
        # set layers such that all cells are considered 'mid'
        (sup0, sup1), (mid0, mid1), (deep0, deep1) = (0,0), (0, np.inf), (0,0)
    # boolean neuron indices:
    supis = (sup0 <= ypos) * (ypos < sup1) # True values are superficial
    midis = (mid0 <= ypos) * (ypos < mid1) # True values are middle
    deepis = (deep0 <= ypos) * (ypos < deep1) # True values are deep
    #otheris = not(supis + midis + deepis) # True values are other, not needed
    return supis, midis, deepis

def pair_laminarity(nids, ypos, trackabsname, pairs, inclusive=False):
    """Color cell pairs according to whether they're superficial, deep, or other. If
    inclusive, label a pair as superficial if *either* of the cells are superficial. Ditto
    for deep. This means many pairs will be counted as both superficial and deep.
    Return RGB colours and indices into self.pairs"""
    # y positions of all nids:
    supis, midis, deepis = laminarity(ypos, trackabsname)
    npairs = len(pairs)
    c = np.empty((npairs, 3), dtype=float) # color RGB array
    cc = mpl.colors.colorConverter
    REDRGB = cc.to_rgb('r')
    GREENRGB = cc.to_rgb('g')
    BLUERGB = cc.to_rgb('b')
    YELLOWRGB = cc.to_rgb('y')
    c[:] = YELLOWRGB # init to yellow, other pairs remain yellow
    if inclusive:
        for i, (ni0, ni1) in enumerate(pairs):
            if supis[ni0] or supis[ni1]:
                c[i] = REDRGB # at least one cell is superficial
            if midis[ni0] or midis[ni1]:
                c[i] = GREENRGB # at least one cell is middle
            if deepis[ni0] or deepis[ni1]:
                c[i] = BLUERGB # at least one cell is deep
    else:
        for i, (ni0, ni1) in enumerate(pairs):
            if supis[ni0] and supis[ni1]:
                c[i] = REDRGB # both cells are superficial
            if midis[ni0] and midis[ni1]:
                c[i] = GREENRGB # both cells are middle
            if deepis[ni0] and deepis[ni1]:
                c[i] = BLUERGB # both cells are deep
    # overwrite boolean neuron indices with boolean pair indices:
    supis, = np.where((c == REDRGB).all(axis=1))
    midis, = np.where((c == GREENRGB).all(axis=1))
    deepis, = np.where((c == BLUERGB).all(axis=1))
    otheris, = np.where((c == YELLOWRGB).all(axis=1))
    return c, supis, midis, deepis, otheris

def rainbow_text(a, x, y, words, colors, **kwargs):
    """
    Take a list of ``words`` and ``colors`` and place them next to each other, with
    words[i] being shown in colors[i]. All keyword arguments are passed to plt.text, so you
    can set the font size, family, etc. Note that horizontal and vertical alignment
    kwargs don't seem to work very well. Also note that although it looks pretty good in the
    QtAgg backend, in some of the important backends, like PNG and PDF, this doesn't space the
    words properly.

    Adapted from Paul Ivanov:
    https://github.com/matplotlib/matplotlib/issues/697#issuecomment-3859591
    """
    f = a.figure
    t = a.transData
    #t = a.transAxes

    # draw horizontal text:
    for w, c in zip(words, colors):
        text = a.text(x, y, " "+w+" ", color=c, transform=t, **kwargs)
        text.draw(f.canvas.get_renderer())
        ex = text.get_window_extent()
        t = mpl.transforms.offset_copy(text._transform, x=ex.width, units='dots')

def mergeuniquedictvals(dicts):
    """Merge a collection of dicts into a single dict, concatenating the list values
    of same-named keys into a single sorted list of unique values for each key"""
    keys = set()
    for d in dicts:
        keys.update(list(d))
    md = {} # merged dict
    for key in keys:
        md[key] = set()
        for d in dicts:
            if key in d: # don't require all dicts to have the same set of keys
                md[key].update(d[key])
        md[key] = sorted(md[key])
    return md

def commontres(t0, y0, t1, y1):
    """Return signals y0 and y1 with common time resolution of overlap periods.
    y0 and y1 can be of arbitrary dimension, but their last dimension must correspond to
    their t0 and t1 timepoints respectively. Assume timepoints t0 and t1 are sorted"""
    assert y0.shape[-1] == len(t0)
    assert y1.shape[-1] == len(t1)

    # find time overlap of both signals:
    mint = max(t0[0], t1[0])
    maxt = min(t0[-1], t1[-1])
    if mint > maxt:
        raise RuntimeError("time series don't overlap")
    trange = mint, maxt # overlapping time range
    i, j = t0.searchsorted(trange)
    t0 = t0[i:j]
    y0 = y0[..., i:j]
    i, j = t1.searchsorted(trange)
    t1 = t1[i:j]
    y1 = y1[..., i:j]
    # choose the lower resolution signal as the reference over the overlapping interval:
    if len(t1) < len(t0):
        # t1 is lower res, find where t1 fits into t0:
        t0i = t0.searchsorted(t1)
        if t0i[-1] >= len(t0):
            t0i[-1] = len(t0) - 1 # prevent right edge out of bounds into t0 and y0
        t = t0[t0i]
        y0 = y0[..., t0i]
    else:
        # t0 is lower res, find where t0 fits into t1:
        t1i = t1.searchsorted(t0)
        if t1i[-1] >= len(t1):
            t1i[-1] = len(t1) - 1 # prevent right edge out of bounds into t1 and y1
        t = t1[t1i]
        y1 = y1[..., t1i]

    return t, y0, y1

def parse_source(source):
    """Collect recordings from source, whether a track, a list of recordings, or a dict of
    track:rid keyvals. Return all collected recordings in a list, as well as the set of all
    associated tracks. Both are sorted in order of their absnames"""
    from track import Track # do this here to prevent circular import
    uns = get_ipython().user_ns
    if type(source) == Track:
        rids = uns['BSMSNSDBRIDS'][source.absname]
        recs = [ source.r[rid] for rid in rids ]
    elif type(source) == list: # assume it's a list of recordings
        recs = source
    elif type(source) == dict: # assume it's a dict of {animal.tr: [rids]} key-value pairs
        recs = []
        for animal_track, rids in source.items():
            animalname, trackname = animal_track.split('.')
            tr = uns[animalname].__getattribute__(trackname)
            recs.extend([ tr.r[rid] for rid in rids ])
    else:
        raise ValueError('unknown source type %r' % type(source))
    # set of all associated tracks:
    tracks = list(set([ rec.tr for rec in recs ]))
    # sort recordings by their absnames:
    recis = np.argsort([ rec.absname for rec in recs ])
    recs = [ recs[reci] for reci in recis ]
    # sort tracks by their absnames:
    try:
        trackis = np.argsort([ track.absname for track in tracks ])
    except AttributeError: # track is None and has no .absname?
        trackis = np.arange(len(tracks))
    tracks = [ tracks[tracki] for tracki in trackis ]
    return recs, tracks

def get_ssnids(recs, stranges=None, kind='active'):
    """Return superset (and sets) of IDs of active or all (at least 1 spike) neurons of all
    recordings in recs (all from the same track). Potentially constrict to spike tranges
    only"""
    recsecnids = [] # holds arrays of  nids of each recording section
    track = recs[0].tr
    # make sure they're all from the same track:
    for rec in recs:
        assert rec.tr == track
    if stranges == None:
        stranges = [None]*len(recs)
    # collect nids for each spike trange:
    for rec, strange in zip(recs, stranges):
        if strange == None:
            tranges = None
        else:
            tranges = [np.asarray(strange)]
        recsecnids.append(rec.get_nids(tranges=tranges, kind=kind))
    ssnids = np.unique(np.hstack(recsecnids)) # superset of nids from rec sections
    return ssnids, recsecnids

def get_nids(recs, tracks, kind=None):
    """Given sorted lists of recordings and unique tracks, return lists of sorted nids, one
    per track, according to value of kind ('active', 'quiet' or 'all'). If kind is 'active' or
    'quiet', degree of activity in each track is measured across only the recordings specified
    for that track, not across all recordings in that track"""
    if kind == None:
        kind = 'active'
    validkinds = ['active', 'quiet', 'all']
    if kind not in validkinds:
        raise ValueError('kind must be one of %r' % validkinds)

    if kind == 'all':
        if len(tracks) > 1: # recordings from multiple tracks, collect nids per track
            nidss = [ sorted(track.alln) for track in tracks ]
        else: # all recordings from same track
            if len(recs) == 1:
                nids = sorted(recs[0].alln) # all of single recording's nids
            elif len(recs) > 1:
                nids = sorted(tracks[0].alln) # all of single parent track's nids
            nidss = [ nids ]
        return nidss # list of nids lists

    # kind is 'active' or 'quiet'
    nidss = []
    uns = get_ipython().user_ns
    MINRATE = uns['MINRATE']
    totaldtsec = 0.0
    for track in tracks:
        dtsec = 0.0
        if track == None and len(recs) == 1:
            allnids = np.asarray(sorted(recs[0].alln))
        else:
            allnids = np.asarray(sorted(track.alln)) # all nids in this track
        nspikes = {}.fromkeys(allnids, 0) # nid to nspikes mapping, init each entry to 0
        for rec in recs:
            if rec.tr != track: # recording doesn't belong to this track
                continue
            for n in rec.alln.values():
                nspikes[n.id] += n.nspikes
            dtsec += rec.dtsec
        rates = np.asarray([ nspikes[nid] / dtsec for nid in allnids ])
        if kind == 'active':
            nids = allnids[rates >= MINRATE]
        else: # kind == 'quiet'
            nids = allnids[rates < MINRATE]
        nidss.append(nids)
        totaldtsec += dtsec
    print('totaldtsec =', totaldtsec)
    return nidss # list of nids lists

def g(x0, sx, x):
    """1-D Gaussian"""
    return np.exp( -(x-x0)**2 / (2*sx**2) )

def g2(x0, y0, sx, sy, x, y):
    """2-D Gaussian"""
    arg = -(x-x0)**2 / (2*sx**2) - (y-y0)**2 / (2*sy**2)
    return np.exp(arg)

def g3(x0, y0, z0, sx, sy, sz, x, y, z):
    """3-D Gaussian"""
    return np.exp( -(x-x0)**2 / (2*sx**2) - (y-y0)**2 / (2*sy**2) - (z-z0)**2 / (2*sz**2) )

def eucd(coords):
    """Generates Euclidean distance matrix from a
    sequence of n m-dimensional coordinates. Nice and fast.
    Written by Willi Richert
    Taken from:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/498246
    on 2006/11/11
    """
    coords = np.asarray(coords)
    n, m = coords.shape
    delta = np.zeros((n, n), dtype=np.float64)
    for d in xrange(m):
        data = coords[:, d]
        delta += (data - data[:, np.newaxis]) ** 2
    return np.sqrt(delta)

def maxabs(a, axis=None):
    """Return slice of a, keeping only those values that are furthest away from 0 along axis"""
    maxa = a.max(axis=axis)
    mina = a.min(axis=axis)
    p = abs(maxa) > abs(mina) # bool, or indices where +ve values win
    n = abs(mina) > abs(maxa) # bool, or indices where -ve values win
    if axis == None:
        if p: return maxa
        else: return mina
    shape = list(a.shape)
    shape.pop(axis)
    out = np.zeros(shape, dtype=a.dtype)
    out[p] = maxa[p]
    out[n] = mina[n]
    return out

def argfwhm(a, exti, fraction=0.5, method='inner'):
    """Find timepoints of full width half max (or whatever fraction is) around extremum
    at index exti in 1D array a. If method is 'inner', search from peak outwards. If method
    is 'outer', search from outer edges of a inwards."""
    #a = abs(a)
    fm = a[exti] * fraction # fraction of max
    d = a - fm
    lis = np.diff(np.sign(d[:exti])).nonzero()[0]
    ris = np.diff(np.sign(d[exti:])).nonzero()[0] + exti + 1
    if len(lis) == 0 or len(ris) == 0:
        # signal doesn't change enough on either side of exti to calculate FWHM
        raise ValueError("exti %d has no FWHM" % exti)
    if method == 'inner':
        # return rightmost of left indices, and leftmost of right indices:
        return lis[-1], ris[0]
    elif method == 'outer':
        # return leftmost of left indices, and rightmost of right indices:
        return lis[0], ris[-1]
    else:
        raise RuntimeError('unknown method %r' % method)

def sparseness(x):
    """Sparseness measure, from Vinje and Gallant, 2000. This is basically 1 minus the ratio
    of the square of the sums over the sum of the squares of the values in signal x"""
    n = len(x)
    return (1 - (x.sum()/n)**2 / np.sum((x**2)/n)) / (1 - 1/n)
