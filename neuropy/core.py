"""Miscellaneous functions and classes"""

from __future__ import division

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

import matplotlib as mpl
import matplotlib.cm
import pylab as pl
from pylab import get_current_fig_manager as gcfm

CODEKIND = 'binary'
CODETRES = 20000 # us
CODEPHASE = 0 # deg
CODEWORDLEN = 10 # in bits

TAB = '    ' # 4 spaces
EPOCH = datetime.datetime(1899, 12, 30, 0, 0, 0) # epoch for datetime stamps in .ptcs


class dictattr(dict):
    """Dictionary with attribute access. Copied from dimstim.Core"""
    def __init__(self, *args, **kwargs):
        super(dictattr, self).__init__(*args, **kwargs)
        for k, v in kwargs.iteritems():
            self.__setitem__(k, v) # call our own __setitem__ so we get keys as attribs even on kwarg init
    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError, '%r object has no attribute %r' % ('dictattr', key)
    
    def __setattr__(self, key, val):
        self[key] = val
    
    def __setitem__(self, key, val):
        super(dictattr, self).__setitem__(key, val)
        if key.__class__ == str and not key[0].isdigit(): # key isn't a number or a string starting with a number
            key = key.replace(' ', '_') # get rid of any spaces
            self.__dict__[key] = val # make the key show up as an attrib upon dir()


class PTCSHeader(object):
    """Polytrode clustered spikes file header"""
    def __init__(self):
        self.VER2FUNC = {1: self.read_ver_1, 2: self.read_ver_2} # call the appropriate method

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
            (source file name, probably .srf, padded with null bytes if needed for 8 byte alignment)
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
        """Same as version 1. NVS created some version 1 files incorrectly, and
        incremented to version 2 for the correctly exported ones"""
        return self.read_ver_1(f)


class PTCSNeuronRecord(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, header):
        self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2} # call the appropriate method
        self.header = header
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[self.header.nsamplebytes]

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
        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes) # spike timestamps (us)

    def read_wave(self, f):
        """Read wavedata/wavestd bytes"""
        nbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nwavedata/nwavestd bytes, padded
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
    

class PopulationRaster(object):
    """A population spike raster plot. nids are indices of neurons to
    plot in the raster, in order from bottom to top.
    jumpts are a sequence of timepoints (in us) that can then be quickly cycled
    through in the plot using keyboard controls.
    Defaults to absolute time origin (when acquisition began)"""
    def __init__(self, recording=None, experiments=None, nids=None,
                 jumpts=None, binwidth=None, relativet0=False, units='msec',
                 publication=False):
        self.r = recording
        if experiments == None:
            self.e = recording.e # dictionary
        else:
            self.e = experiments # should also be a dict
        if binwidth == None:
            self.binwidth = CODETRES
        else:
            self.binwidth = binwidth
        assert self.binwidth >= 10000
        self.plotbinedges = False # keyboard controlled
        if relativet0: # set time origin to start of first experiment
            firstexp = min(self.e.keys())
            self.t0 = self.e[firstexp].trange[0] # start time of the first experiment
        else: # use absolute time origin
            self.t0 = 0 # leave the time origin at when acquisition began
        experimentmarkers = [] # a flat list of all experiment start and stop times, in sorted order
        for e in self.e.values():
            experimentmarkers.extend(e.trange)
        self.experimentmarkers = np.asarray(experimentmarkers) - self.t0 # make 'em relative to t0
        self.experimentmarkers.sort() # just in case exps weren't in sorted order for some reason

        self.jumpts = np.asarray(jumpts)
        self.jumpts.sort() # just in case jumpts weren't in sorted order for some reason

        units2tconv = {'usec': 1e0, 'msec': 1e3, 'sec': 1e6}
        self.units = units
        self.tconv = units2tconv[units]

        self.publication = publication

        self.neurons = self.r.n # still a dictattr
        if nids != None:
            self.nids = nids
        else:
            self.nids = self.r.n.keys()
            self.nids.sort() # keep it tidy

    def plot(self, left=None, width=200000):
        """Plots the raster, units are us wrt self.t0"""
        if left == None:
            try:
                left = self.experimentmarkers[0] # init left window edge to first exp marker, ie start of first experiment
            except IndexError: # ain't no experiments, no associated markers
                left = min([ neuron.spikes[0] for neuron in self.r.n.values() ]) # set it to the earliest spike in the population
        try:
            self.f
        except AttributeError: # prepare the fig if it hasn't been done already
            figheight = 1.25+0.2*len(self.nids)
            self.f = pl.figure(figsize=(14, figheight))
            self.a = self.f.add_subplot(111)
            self.formatter = NeuropyScalarFormatter() # better behaved tick label formatter
            self.formatter.thousandsSep = ',' # use a thousands separator
            if not self.publication:
                self.a.xaxis.set_major_locator(NeuropyAutoLocator()) # better behaved tick locator
                self.a.xaxis.set_major_formatter(self.formatter)
                self.a.set_yticks([]) # turn off y axis
            gcfm().window.setWindowTitle(lastcmd())
            #self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100)) # create a long tooltip with newline to get around bug where newlines aren't recognized on subsequent self.tooltip.SetTip() calls
            #self.tooltip.Enable(False) # leave disabled for now
            #self.tooltip.SetDelay(0) # set popup delay in ms
            #gcfm().canvas.SetToolTip(self.tooltip) # connect the tooltip to the canvas
            self.a.set_xlabel('time (%s)' % self.units)
            if not self.publication:
                self.yrange = (0, len(self.nids))
            else:
                self.a.set_ylabel('cell index') # not really cell id, it's the ii (the index into the id)
                self.yrange = (-0.5, len(self.nids)-0.5)
            self.a.set_ylim(self.yrange)
            #aheight = min(0.025*len(self.nids), 1.0)
            bottominches = 0.75
            heightinches = 0.15+0.2*len(self.nids)
            bottom = bottominches / figheight
            height = heightinches / figheight
            if not self.publication:
                self.a.set_position([0.02, bottom, 0.96, height])
            else:
                self.a.set_position([0.05, bottom, 0.96, height])
            self.f.canvas.mpl_connect('motion_notify_event', self._onmotion)
            self.f.canvas.mpl_connect('key_press_event', self._onkeypress)

        self.left = left
        self.width = width
        # plot experiment start and endpoints
        for etrange in self.experimentmarkers.reshape(-1, 2): # reshape the flat array into a new nx2, each row is a trange
            estart = etrange[0]-self.t0
            eend = etrange[1]-self.t0
            if left <=  estart and estart <= left+width: # experiment start point is within view
                startlines = self.a.vlines(x=estart/self.tconv, ymin=self.yrange[0], ymax=self.yrange[1], color='k', linestyle='-') # marks exp start, convert to ms
                #startlines[0].set_color((0, 1, 0)) # set to bright green, can't do this for line coll
            if left <= eend and eend <= left+width: # experiment end point is within view
                endlines = self.a.vlines(x=eend/self.tconv, ymin=self.yrange[0], ymax=self.yrange[1], color='k', linestyle='-') # marks exp end, convert t
                #endlines[0].set_color((1, 0, 0)) # set to bright red, can't do this for line coll
        # plot the bin edges. Not taking into account self.t0 for now, assuming it's 0
        if self.plotbinedges:
            leftbinedge = (left // self.binwidth + 1)*self.binwidth
            binedges = np.arange(leftbinedge, left+width, self.binwidth)
            binlines = self.a.vlines(x=binedges/self.tconv, ymin=self.yrange[0], ymax=self.yrange[1], color='b', linestyle=':') # convert t
        # plot the rasters
        for nii, ni in enumerate(self.nids):
            neuron = self.neurons[ni]
            x = (neuron.cut((self.t0+left, self.t0+left+width)) - self.t0) / self.tconv # make spike times always relative to t0, convert t
            if not self.publication:
                self.a.vlines(x=x, ymin=nii, ymax=nii+1, color='k', linestyle='-')
            else:
                self.a.vlines(x=x, ymin=nii-0.5, ymax=nii+0.5, color='k', linestyle='-')
        self.a.set_xlim(left/self.tconv, (left+width)/self.tconv) # convert t

    def _panx(self, npages=None, left=None):
        """Pans the raster along the x axis by npages, or to position left"""
        self.a.lines = [] # first, clear all the vlines, this is easy but a bit innefficient, since we'll probably be redrawing most of the ones we just cleared
        if left != None: # use left
            self.plot(left=left, width=self.width)
        else: # use npages instead
            self.plot(left=self.left+self.width*npages, width=self.width)
        self.f.canvas.draw() # redraw the figure

    def _zoomx(self, factor):
        """Zooms the raster along the x axis by factor"""
        self.a.lines = [] # first, clear all the vlines, this is easy but a bit innefficient, since we'll probably be redrawing most of the ones we just cleared
        centre = (self.left + self.left+self.width) / 2.0
        width = self.width / float(factor)
        left = centre - width / 2.0
        self.plot(left=left, width=width)
        self.f.canvas.draw() # redraw the figure

    def _goto(self):
        """Bring up a dialog box to jump to timepoint, mark it with a dotted line"""
        ted = wx.TextEntryDialog(parent=None, message='Go to timepoint (ms):', caption='Goto',
                                 defaultValue=str(intround(self.left / self.tconv)), #wx.EmptyString,
                                 style=wx.TextEntryDialogStyle, pos=wx.DefaultPosition)
        if ted.ShowModal() == wx.ID_OK: # if OK button has been clicked
            response = ted.GetValue()
            try:
                left = float(response)
                self.plot(left=left*self.tconv, width=self.width)
                self.f.canvas.draw() # redraw the figure
            except ValueError: # response wasn't a valid number
                pass

    def _cyclethousandssep(self):
        """Cycles the tick formatter through thousands separators"""
        if self.formatter.thousandsSep == ',':
            self.formatter.thousandsSep = ' '
        elif self.formatter.thousandsSep == ' ':
            self.formatter.thousandsSep = None
        else:
            self.formatter.thousandsSep = ','
        self.f.canvas.draw() # redraw the figure

    def _togglebinedges(self):
        """Toggles plotting of bin edges"""
        self.plotbinedges = not self.plotbinedges
        self._panx(npages=0) # replot and redraw by panning by 0

    def _onmotion(self, event):
        """Called during mouse motion over figure. Pops up neuron and
        experiment info in a tooltip when hovering over a neuron row."""
        if event.inaxes: # if mouse is inside the axes
            nii = int(math.floor(event.ydata)) # use ydata to get index into sorted list of neurons
            ni = self.nids[nii]
            neuron = self.neurons[ni]
            currentexp = None
            for e in self.e.values(): # for all experiments
                estart = (e.trange[0]-self.t0)/self.tconv
                eend = (e.trange[1]-self.t0)/self.tconv
                if estart < event.xdata  < eend:
                    currentexp = e
                    break # don't need to check any of the other experiments
            tip = 't: %.3f ms\n' % event.xdata # print timepoint down to nearest us, in units of ms
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
        #print key
        # you can also just use the backend-neutral event.key, but that doesn't recognize as many keypresses, like pgup, pgdn, etc.
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
                i = self.jumpts.searchsorted(self.left, side='left') # current position of left edge of the window in jumpts list
                i = max(0, i-1) # decrement by 1, do bounds checking
                self._panx(left=self.jumpts[i])
            elif key == ord(']'): # skip forwards to next jump point
                i = self.jumpts.searchsorted(self.left, side='right') # current position of left edge of the window in jumpts list
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
                i = self.experimentmarkers.searchsorted(self.left, side='left') # current position of left edge of the window in experimentmarkers list
                i = max(0, i-1) # decrement by 1, do bounds checking
                self._panx(left=self.experimentmarkers[i])
            elif key == wx.WXK_RIGHT: # skip forwards to next experiment marker
                i = self.experimentmarkers.searchsorted(self.left, side='right') # current position of left edge of the window in experimentmarkers list
                i = min(i, len(self.experimentmarkers)-1) # bounds checking
                self._panx(left=self.experimentmarkers[i])
            elif key == wx.WXK_UP: # zoom in faster
                self._zoomx(3.0)
            elif key == wx.WXK_DOWN: # zoom out faster
                self._zoomx(1/3.0)


class Codes(object):
    """A 2D array where each row is a neuron code, and each column
    is a binary population word for that time bin, sorted LSB to MSB from top to bottom.
    neurons is a list of Neurons, also from LSB to MSB. Order in neurons is preserved."""
    def __init__(self, neurons=None, tranges=None, kind=CODEKIND, tres=CODETRES,
                 phase=CODEPHASE, shufflecodes=False):
        self.neurons = neurons
        self.tranges = tolist(tranges)
        self.kind = kind
        self.tres = tres
        self.phase = phase
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
        self.s = [] # stores the corresponding spike times for each neuron, for reference
        self.c = [] # stores the 2D code array
        # append neurons in their order in self.neurons, store them LSB to MSB from top to
        # bottom
        for neuron in self.neurons:
            codeo = neuron.code(tranges=self.tranges, kind=self.kind, tres=self.tres,
                                phase=self.phase)
            # build up nested list (ie, 2D) of spike times, each row will have different
            # length:
            self.s.append(codeo.s)
            if self.shufflecodes:
                c = codeo.c.copy() # make a copy (wanna leave the codeo's codetrain untouched)
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
        self.c = np.concatenate(self.c).reshape(nneurons, nbins)

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
    
class CodeCorrPDF(object):
    """A PDF of the correlations of the codes of all cell pairs (or of all cell pairs within
    some torus of radii R=(R0, R1) in um) in this Recording. See Schneidman2006 fig 1d"""
    def __init__(self, recording=None, experiments=None, nids=None, kind=CODEKIND,
                 tres=CODETRES, phase=CODEPHASE):
        self.r = recording
        if experiments != None:
            try:
                assert experiments.__class__ == dictattr
            except AssertionError: # maybe it's a seq of exp ids?
                eids = toiter(experiments)
                experiments = dictattr()
                for eid in eids:
                    experiments[eid] = self.r.e[eid]
        self.e = experiments # save it, should be a dictattr if not None
        if self.e != None: # specific experiments were specified
            self.tranges = [ e.trange for e in self.e.values() ]
        else:
            self.tranges = [ self.r.trange ] # use the Recording's trange
        if nids == None:
            nids = self.r.n.keys()
        self.nids = nids
        self.kind = kind
        self.tres = tres
        self.phase = phase
    '''
    # this was used to save on calc time by seeing if a CCPDF object with the same attribs
    # had already been calc'd, seems dumb and unsafe, commented out
    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their n and c attribs, if they exist, to prevent comparing them below,
        # since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd]
          for key in ['corrs', 'n', 'c', 'crange', 'nbins', 'normed']
          if d.has_key(key) ]
        if self.__class__ == other.__class__ and selfd == otherd:
            return True
        else:
            return False
    '''
    def calc(self, R=None, shuffleids=False):
        """Calculates self constrained to self.tranges and torus described by R"""
        if R != None:
            assert len(R) == 2 and R[0] < R[1]  # should be R = (R0, R1) torus
        self.R = R
        self.shuffleids = shuffleids
        nids = self.nids
        nneurons = len(nids)
        # it's more efficient to precalculate the means and stds of each cell's codetrain,
        # and then reuse them in calculating the correlation coefficients
        # store each code mean in a dict:
        means = dict( ( nid, self.r.n[nid].code(tranges=self.tranges,
                                                kind=self.kind,
                                                tres=self.tres,
                                                phase=self.phase).c.mean() ) for nid in nids )
        # store each code std in a dict:
        stds  = dict( ( nid, self.r.n[nid].code(tranges=self.tranges,
                                                kind=self.kind,
                                                tres=self.tres,
                                                phase=self.phase).c.std() ) for nid in nids ) 
        if self.shuffleids:
        # shuffled neuron ids, this is a control to see if it's the locality of neurons
        # included in the analysis, or the number of neurons included that's important. It
        # seems that both are.
            snids = shuffle(nids)
        else:
            snids = nids
        self.corrs = []
        for nii1 in range(nneurons):
            for nii2 in range(nii1+1, nneurons):
                ni1 = nids[nii1]
                sni1 = snids[nii1]
                ni2 = nids[nii2]
                sni2 = snids[nii2]
                # calc the pair's code correlation if there's no torus specified, or if
                # the pair's separation falls with bounds of specified torus:
                if R == None or (self.R[0] < dist(self.r.n[sni1].pos, self.r.n[sni2].pos)
                                 < self.R[1]):
                    code1 = self.r.n[ni1].code(tranges=self.tranges, kind=self.kind,
                                               tres=self.tres, phase=self.phase).c
                    code2 = self.r.n[ni2].code(tranges=self.tranges, kind=self.kind,
                                               tres=self.tres, phase=self.phase).c
                    # (mean of product - product of means) / by product of stds
                    cc = (((code1 * code2).mean() - means[ni1] * means[ni2])
                          / (stds[ni1] * stds[ni2]))
                    self.corrs.append(cc)
        self.corrs = np.array(self.corrs)
        self.npairs = len(self.corrs)
        '''
        # simpler, but slower way:
        self.corrs = [ self.r.codecorr(nids[nii1], nids[nii2], tranges=self.tranges,
                                       kind=self.kind, self.tres, self.phase)
                       for nii1 in range(0,nneurons) for nii2 in range(nii1+1,nneurons) ]
        '''
    def plot(self, figsize=(7.5, 6.5), crange=[-0.1, 0.2], limitstats=True,
             nbins=30, normed='pdf'):
        """Plots the corrs. If limitstats, the stats displayed exclude any corr values that
        fall outside of crange"""
        self.crange = crange
        nbins = max(nbins, 2*intround(np.sqrt(self.npairs)))
        self.nbins = nbins
        self.normed = normed
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        try: # figure out the bin edges
            bins = np.linspace(start=self.crange[0], stop=self.crange[1], num=self.nbins,
                               endpoint=True)
        except TypeError: # self.crange is None, let histogram() figure out the bin edges
            bins = self.nbins
        self.n, self.c = histogram(self.corrs, bins=bins, normed=self.normed)

        if limitstats:
            corrs = self.corrs[(self.corrs >= crange[0]) * (self.corrs <= crange[1])]
            n, c = histogram(corrs, bins=bins, normed=self.normed)
        else:
            corrs = self.corrs
            n = self.n
            c = self.c
        self.mean = np.mean(corrs)
        self.median = np.median(corrs)
        argmode = n.argmax()
        self.mode = np.mean([c[argmode], c[argmode + 1]]) # find middle of tallest bin
        self.stdev = np.std(corrs)

        try:
            barwidth = (self.crange[1] - self.crange[0]) / float(self.nbins)
        except TypeError: # self.crange is None, take width of first bin in self.c
            barwidth = self.c[1] - self.c[0]
        a.bar(left=c, height=n, width=barwidth, bottom=0, color='k',
              yerr=None, xerr=None, ecolor='k', capsize=3)
        try:
            a.set_xlim(self.crange)
        except TypeError: # self.crange is None
            pass
        gcfm().window.setWindowTitle(lastcmd())
        #titlestr = 'neuron pair code correlation pdf'
        #titlestr += '\n%s' % lastcmd()
        titlestr = '%s' % lastcmd()
        a.set_title(titlestr)
        
        if self.normed:
            if self.normed == 'pmf':
                a.set_ylabel('probability mass')
            else: # self.normed = 'pdf'
                a.set_ylabel('probability density')
        else:
            a.set_ylabel('count')
        a.set_xlabel('correlation coefficient')
        
        # add stuff to top right of plot:
        a.text(0.99, 0.99, '%s\n'
                           'mean = %.3f\n'
                           'median = %.3f\n'
                           'mode = %.3f\n'
                           'stdev = %.3f\n'
                           'R = %r\n'
                           'ratethresh = %.2f Hz\n'
                           'nneurons = %d\n'
                           'npairs = %d\n'
                           'dt = %d min'
                           % (self.r.name, self.mean, self.median, self.mode, self.stdev,
                              self.R, get_ipython().user_ns['QUIETMEANRATETHRESH'],
                              len(self.nids), self.npairs, intround(self.r.dtmin)),
                           transform = a.transAxes,
                           horizontalalignment='right',
                           verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self
'''
class CanvasFrame(wx.Frame):
    """A minimal wx.Frame containing a matplotlib figure"""
    def __init__(self, title='frame', size=(550, 350)):
        wx.Frame.__init__(self, None, -1, title=title, size=size)
        self.SetBackgroundColour(wx.NamedColor("WHITE"))
        self.figure = pl.figure.Figure(figsize=(5, 4), dpi=100)
        #self.axes = self.figure.add_subplot(111)
        #t = np.arange(0.0, 3.0, 0.01)
        #s = sin(2*pi*t)
        #self.axes.plot(t,s)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)

        # Capture the paint message, slows frame down a little, can be commented out
        #wx.EVT_PAINT(self, self.OnPaint)
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        self.SetSizer(self.sizer)
        self.Fit()
    def OnPaint(self, event):
        self.canvas.draw()
'''

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
                 nids=None, ts=None, scale=2.0):
        NeuropyWindow.__init__(self, parent)
        self.title = title
        self.rfs = rfs
        self.nids = nids
        self.ts = ts
        self.scale = scale # setting to a float will give uneven sized pixels

        cmap = mpl.cm.jet(np.arange(256), bytes=True) # 8 bit RGBA colormap
        #cmap[:, [0, 1, 2, 3]] = cmap[:, [3, 0, 1, 2]] # 8 bit ARGB colormap
        # from Qt docs, sounds like I should be using ARGB format, but seems like
        # RGBA is the format that works in PyQt4
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

'''
class ReceptiveFieldFrame(wx.Frame):
    """A wx.Frame for plotting a scrollable 2D grid of receptive fields, with neuron and time labels.
    rfs is a list of (nt, width, height) sized receptive fields of uint8 RGB data, one per neuron"""
    def __init__(self, parent=None, id=-1, title='ReceptiveFieldFrame', rfs=None,
                 neurons=None, t=None, scale=2.0):
        self.rfs = rfs
        self.neurons = neurons
        self.t = t
        self.title = title
        wx.Frame.__init__(self, parent=parent, id=id, title=title, style=wx.DEFAULT_FRAME_STYLE)
        self.panel = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)

        self.bitmaps = {}
        for ni, n in enumerate(self.neurons):
            self.bitmaps[ni] = {}
            for ti, t in enumerate(self.t):
                rf = self.rfs[ni][ti]
                im = wx.ImageFromData(width=rf.shape[0], height=rf.shape[1], data=rf.data) # expose rf as databuffer
                im = im.Scale(width=im.GetWidth()*scale, height=im.GetHeight()*scale)
                self.bitmaps[ni][t] = wx.StaticBitmap(parent=self.panel, bitmap=im.ConvertToBitmap())

        #self.Bind(wx.EVT_PAINT, self.OnPaint)
        #self.Bind(wx.EVT_MOUSEWHEEL, self.OnMouseWheel)
        self.__set_properties()
        self.__do_layout()
    """
    def OnPaint(self, event):
        #self.canvas.draw()
        event.Skip()
    """
    def OnMouseWheel(self, event):
        """This could be useful..."""
        pass

    def __set_properties(self):
        self.SetTitle(self.title)
        self.panel.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.panel.SetScrollRate(10, 10)

    def __do_layout(self):
        sizer_1 = wx.GridSizer(1, 1, 0, 0)
        grid_sizer_1 = wx.FlexGridSizer(rows=len(self.neurons)+1, cols=len(self.t)+1, vgap=2, hgap=2) # add an extra row and column for the text labels
        grid_sizer_1.Add((1, 1), 0, wx.ADJUST_MINSIZE, 0) # spacer in top left corner
        for t in self.t:
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "%sms" % t), 0, wx.ADJUST_MINSIZE, 0) # text row along top
        for ni, n in enumerate(self.neurons):
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "n%d" % n.id), 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ADJUST_MINSIZE, 0) # text down left side
            for t in self.t:
                grid_sizer_1.Add(self.bitmaps[ni][t], 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL, 0)
        self.panel.SetAutoLayout(True)
        self.panel.SetSizer(grid_sizer_1)
        grid_sizer_1.Fit(self.panel)
        #grid_sizer_1.SetSizeHints(self.panel) # prevents the panel from being resized to something smaller than the above fit size
        """
        # might be a more direct way to set these:
        for rowi in range(1, len(self.ns)+1):
            print 'rowi:', rowi
            grid_sizer_1.AddGrowableRow(rowi)
        for coli in range(1, len(self.ts)+1):
            print 'coli:', coli
            grid_sizer_1.AddGrowableCol(coli)
        """
        sizer_1.Add(self.panel, 1, wx.ADJUST_MINSIZE|wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        #sizer_1.SetSizeHints(self) # prevents the frame from being resized to something smaller than the above fit size
        self.Layout()


class NetstateReceptiveFieldFrame(ReceptiveFieldFrame):
    """A wx.Frame for plotting a scrollable 2D grid of netstate receptive fields, with netstate and time labels.
    rfs is a list of (nt, width, height) sized receptive fields of uint8 RGB data, one per netstate"""
    def __init__(self, parent=None, id=-1, title='NetstateReceptiveFieldFrame',
                 rfs=None, intcodes=None, t=None, scale=2.0):
        self.rfs = rfs
        self.intcodes = tolist(intcodes)
        self.t = t
        self.title = title
        wx.Frame.__init__(self, parent=parent, id=id, title=title, style=wx.DEFAULT_FRAME_STYLE)
        self.panel = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)
        self.bitmaps = {}
        for ii, i in enumerate(self.intcodes):
            self.bitmaps[ii] = {}
            for ti, t in enumerate(self.t):
                rf = self.rfs[ii][ti]
                im = wx.ImageFromData(width=rf.shape[0], height=rf.shape[1], data=rf.data) # expose rf as databuffer
                im = im.Scale(width=im.GetWidth()*scale, height=im.GetHeight()*scale)
                self.bitmaps[ii][t] = wx.StaticBitmap(parent=self.panel, bitmap=im.ConvertToBitmap())
        self.__set_properties()
        self.__do_layout()

    def __set_properties(self):
        self.SetTitle(self.title)
        self.panel.SetBackgroundColour(wx.Colour(255, 255, 255))
        self.panel.SetScrollRate(10, 10)

    def __do_layout(self):
        sizer_1 = wx.GridSizer(1, 1, 0, 0)
        grid_sizer_1 = wx.FlexGridSizer(rows=len(self.intcodes)+1, cols=len(self.t)+1, vgap=2, hgap=2) # add an extra row and column for the text labels
        grid_sizer_1.Add((1, 1), 0, wx.ADJUST_MINSIZE, 0) # spacer in top left corner
        for t in self.t:
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "%sms" % t), 0, wx.ADJUST_MINSIZE, 0) # text row along top
        for ii, i in enumerate(self.intcodes):
            grid_sizer_1.Add(wx.StaticText(self.panel, -1, "ns%d" % i), 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ADJUST_MINSIZE, 0) # text down left side
            for t in self.t:
                grid_sizer_1.Add(self.bitmaps[ii][t], 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL, 0)
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
        #print 'means:\n', means
        #print 'pairmeans:\n', pairmeans

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
        #print 'means:', means
        #print 'pairmeans:', pairmeans
        #print '%d iters,' % self.model.iters
        #print 'hi:', self.hi.__repr__()
        #print 'Jij:', self.Jij.__repr__()

        '''
        # Output the distribution
        print "\nFitted model parameters are:\n" + str(self.model.params)
        print "\nFitted distribution is:"
        for j in range(len(self.model.samplespace)):
            x = np.array(self.model.samplespace[j])
            x = (x+1)/2 # convert from -1s and 1s back to 0s and 1s
            print '\tx:%s, p(x):%s' % (x, p[j])
        '''
        '''
        # Now show how well the constraints are satisfied:
        print
        print "Desired constraints:"
        print "\tp['dans'] + p['en'] = 0.3"
        print ("\tp['dans'] + p['" + a_grave + "']  = 0.5").encode('utf-8')
        print
        print "Actual expectations under the fitted model:"
        print "\tp['dans'] + p['en'] =", p[0] + p[1]
        print ("\tp['dans'] + p['" + a_grave + "']  = " + str(p[0]+p[2])).encode('utf-8')
        # (Or substitute "x.encode('latin-1')" if you have a primitive terminal.)
        '''
'''
class Cat15Movie(object):
    """dimstim >= 0.16 Experiments use the dimstim Experiment (subclassed by say, Movie) object directly"""
    def __init__(self, fname=None, name=None, parent=None):
        """Movies don't need parents, they can just exist on their own and be used by anyone"""
        self.level = 5 # level in the hierarchy
        self.fname = fname
        self.name = name
        self.parent = parent # save parent object, might be an Experiment, might not
        if self.name == None and self.fname != None:
            self.path, self.fname = os.path.split(self.fname) # separate path from fname
            self.name = os.path.splitext(self.fname)[0] # extentionless fname
            if self.name not in _data.movies:
                _data.movies[self.name] = self # add self to _data.movies dictattr
        else:
            pass # both self.name and self.fname are None, this happens when executing Cat 15 textheaders, where you init a movie with m = Movie(), and only later assign its fname field. In this case, the .loadCat15exp() method handles adding movies init'd from textheader to the _data.movies dictattr

    def load(self, asarray=True, flip=False):
        """Load movie frames"""
        try:
            self.frames # movie's already been loaded, don't do anything
            return
        except AttributeError:
            pass
        try:
            self.frames = _data.movies[self.name].frames # if a Movie init'd with the same name already has its data loaded, use it
            return
        except AttributeError:
            pass

        self.f = file(self.fname, 'rb') # open the movie file for reading in binary format
        headerstring = self.f.read(5)
        if headerstring == 'movie': # a header has been added to the start of the file
            self.ncellswide, = struct.unpack('H', self.f.read(2)) # 'H'== unsigned short int
            self.ncellshigh, = struct.unpack('H', self.f.read(2))
            self.nframes, = struct.unpack('H', self.f.read(2))
            if self.nframes == 0: # this was used in Cat 15 mseq movies to indicate 2**16 frames, shouldn't really worry about this, cuz we're using slightly modified mseq movies now that we don't have the extra frame at the end that the Cat 15 movies had (see comment in Experiment module), and therefore never have a need to indicate 2**16 frames
                self.nframes = 2**16
            self.offset = self.f.tell() # header is 11 bytes long
        else: # there's no header at the start of the file, set the file pointer back to the beginning and use these hard coded values:
            self.f.seek(0)
            self.ncellswide = self.ncellshigh = 64
            self.nframes = 6000
            self.offset = self.f.tell() # header is 0 bytes long
        self.framesize = self.ncellshigh*self.ncellswide

        # read in all of the frames
        # maybe check first to see if file is > 1GB, if so, _loadaslist() to prevent trying to allocate one huge piece of contiguous memory and raising a MemoryError, or worse, segfaulting
        if asarray:
            self._loadasarray(flip=flip)
        else:
            self._loadaslist(flip=flip)
        leftover = self.f.read() # check if there are any leftover bytes in the file
        if leftover != '':
            pprint(leftover)
            print self.ncellswide, self.ncellshigh, self.nframes
            raise RuntimeError, 'There are unread bytes in movie file %r. Width, height, or nframes is incorrect in the movie file header.' % self.fname
        self.f.close() # close the movie file
        treestr = self.level*TAB + self.fname
        print treestr

    def _loadasarray(self, flip=False):
        self.frames = np.fromfile(self.f, np.uint8, count=self.nframes*self.framesize)
        self.frames.shape = (self.nframes, self.ncellshigh, self.ncellswide)
        self.f.seek(self.offset + self.nframes*self.framesize) # seek to what should be EOF
        if flip:
            self.frames = self.frames[::, ::-1, ::] # flip all frames vertically for OpenGL's bottom left origin

    def _loadaslist(self, flip=False):
        self.frames = []
        for framei in xrange(self.nframes): # one frame at a time...
            frame = np.fromfile(self.f, np.uint8, count=self.framesize) # load the next frame
            frame.shape = (self.ncellshigh, self.ncellswide)
            if flip:
                frame = frame[::-1, ::] # flip all frames vertically for OpenGL's bottom left origin
            self.frames.append(frame)
'''

class NeuropyScalarFormatter(mpl.ticker.ScalarFormatter):
    """Overloaded from mpl.ticker.ScalarFormatter for 4 reasons:
    1) turn off stupid offset
    2) increase maximum possible number of sigfigs
    3) increase +ve and -ve order of magnitude thresholds before switching to scientific notation
    4) keep exponents in engineering notation, ie multiples of 3
    """
    def __init__(self, useOffset=False, useMathText=False):
        # useOffset allows plotting small data ranges with large offsets:
        # for example: [1+1e-9,1+2e-9,1+3e-9]
        # useMathText will render the offset an scientific notation in mathtext
        #super(NeuropyScalarFormatter, self).__init__(useOffset=useOffset, useMathText=useMathText) # can't use this, cuz derived from an old-style class
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
            self.orderOfMagnitude = (oom // 3)*3 # stick to engineering notation, multiples of 3
        elif oom > 6: # increased +ve threshold for sci notation
            self.orderOfMagnitude = (oom // 3)*3 # stick to engineering notation, multiples of 3
        else:
            self.orderOfMagnitude = 0

    def _set_format(self):
        # set the format string to format all the ticklabels
        locs = (np.array(self.locs)-self.offset) / 10**self.orderOfMagnitude+1e-15
        sigfigs = [len(str('%1.10f'% loc).split('.')[1].rstrip('0')) for loc in locs] # '%1.3f' changed to '%1.10f' to increase maximum number of possible sigfigs
        sigfigs.sort()
        self.format = '%1.' + str(sigfigs[-1]) + 'f'
        if self._usetex or self._useMathText: self.format = '$%s$'%self.format

    def pprint_val(self, x):
        xp = (x-self.offset)/10**self.orderOfMagnitude
        if np.absolute(xp) < 1e-8: xp = 0
        s = self.format % xp
        if self.thousandsSep: # add thousands-separating characters
            if s.count('.'): # it's got a decimal in there
                s = re.sub(r'(?<=\d)(?=(\d\d\d)+\.)', self.thousandsSep, s) # use the regexp for floats
            else: # it's an int
                s = re.sub(r'(?<=\d)(?=(\d\d\d)+$)', self.thousandsSep, s) # use the regexp for ints
        return s


class NeuropyAutoLocator(mpl.ticker.MaxNLocator):
    """A tick autolocator that generates more ticks than the standard mpl autolocator"""
    def __init__(self):
        #mpl.ticker.MaxNLocator.__init__(self, nbins=9, steps=[1, 2, 5, 10]) # standard autolocator
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
    '''Class decorator for making a class behave as a Java (non-static) inner
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
    '''
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
    """Round to the nearest integer, return an integer. Works on arrays,
    saves on parentheses, nothing more"""
    if iterable(n): # it's a sequence, return as an int64 array
        return np.int64(np.round(n))
    else: # it's a scalar, return as normal Python int
        return int(round(n))

def pad0s(val, ndigits=2):
    """Returns a string rep of val, padded with enough leading 0s
    to give you a string rep with ndigits in it"""
    if val.__class__ != int:
        raise ValueError, '%r isn\'t an int' % val
    val = str(val)
    nzerostoadd = ndigits - len(val)
    val = '0'*nzerostoadd + val
    return val

def txtdin2binarydin(fin, fout):
    """Converts a csv text .din file to an int64 binary .din file"""
    fi = file(fin, 'r') # open the din file for reading in text mode
    fo = file(fout, 'wb') # for writing in binary mode
    for line in fi:
        line = line.split(',')
        '''
        # for old NVS display, converts from NVS condition numbers (which increment with repeats) to dimstim sweepis (which don't)
        nruns = 18
        line[1] = int(line[1]) % nruns
        '''
        fo.write( struct.pack('@qq', int(line[0]), int(line[1])) ) # write both values out as a C long longs, using the system's native ('@') byte order
    fi.close()
    fo.close()
    print 'Converted ascii din: %r to binary din: %r' % (fin, fout)

def convertalltxtdin2binarydin(path=None):
    """Converts all text .csv din files in path (or cwd) to 64 bit binary .din files of the same name"""
    if path == None:
        path = os.getcwd()

    listing = os.listdir(path)
    dinfnames = []

    for fname in listing:
        if fname.endswith('.csv'):
            dinfnames.append(fname[:-len('.csv')]) # text din filenames without the .csv extension

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
                print newfname
                os.rename(os.path.join(path, fname), os.path.join(path, newfname))

def csv2binary(fin, multiplier=1e6, skipfirstline=True):
    """Converts spike data in a csv file, with cells in the columns and times down the rows,
    into int64 binary files, one for each neuron. Takes csv values and multiplies them by
    multiplier before saving"""
    fin = os.path.normpath(fin)
    fi = file(fin, 'r') # open csv file for reading in text mode
    print 'Exporting %s to:' % fi.name
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
    tail = os.path.split(path)[-1].replace(' ', '_') # just the extensionless filename, replace spaces with underscores
    for ni, neuron in enumerate(data):
        fname = os.path.join(path, tail) + '_t' + pad0s(ni, ndigits=len(str(nneurons))) + '.spk'
        fo = file(fname, 'wb') # for writing in binary mode
        print fo.name
        for spiketime in neuron: # write each spiketime to the file, there should be a more streamlined way to do this
            fo.write( struct.pack('@q', spiketime) ) # write the value out as a C long long, using the system's native ('@') byte order
        fo.close()

def warn(msg, level=2, exit_val=1):
    """Standard warning printer. Gives formatting consistency. Stolen from IPython.genutils"""
    if level > 0:
        header = ['', '', 'WARNING: ', 'ERROR: ', 'FATAL ERROR: ']
        print >> sys.stderr, '%s%s' % (header[level],msg)
        if level == 4:
            print >> sys.stderr,'Exiting.\n'
            sys.exit(exit_val)
'''
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
    """Returns the input list minus any repeated objects it may have had. Also defined in dimstim"""
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
    """Convert to list. If input is a dict, returns its values. If it's already a list,
    returns it. Otherwise, input is returned in a list."""
    if type(x) == dict:
        return list(x.values())
    elif type(x) == list:
        return x
    else:
        return [x] # stick it in a list

def to2d(arr):
    """Converts a 1D array to a 2D array with just a singleton row
    If arr is already 2D, just returns it. If it's anything more than 2D, raises an error"""
    nd = arr.ndim
    assert nd in [1, 2], 'array rank > 2'
    if nd == 1:
        arr = arr.reshape(1, -1)
    return arr

def joinpath(pathlist):
    """Unlike os.path.join(), takes a list of path segments, returns them joined in a string
    with local separators"""
    path = ''
    for segment in pathlist:
        path = os.path.join(path, segment)
    return path

def dist(a, b):
    """Returns the Euclidean distance between two coordinates.
    Both a and b must be (x, y) tuples"""
    return math.sqrt( (b[0]-a[0])**2 + (b[1]-a[1])**2 )

def approx(a, b, rtol=1.e-14, atol=1.e-14):
    """Returns a boolean array describing which components of a and b are equal
    subject to given tolerances. The relative error rtol must be positive and << 1.0
    The absolute error atol comes into play for those elements of y that are very
    small or zero; it says how small x must be also. Copied and modified from
    numpy.allclose()"""
    x = np.array(a, copy=False)
    y = np.array(b, copy=False)
    #print x.shape
    #print y.shape
    return np.less(np.absolute(x-y), atol + rtol * np.absolute(y))

def histogram(a, bins=10, range=None, normed=False):
    """Builds a histogram, stolen from numpy.histogram(), modified to allow
    normed='pdf' or normed='pmf' (prob mass function)"""
    a = np.asarray(a).ravel()
    if not iterable(bins):
        if range is None:
            range = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins, endpoint=False)
    n = np.sort(a).searchsorted(bins)
    n = np.concatenate([n, [len(a)]])
    n = n[1:]-n[:-1]
    if normed:
        if normed == 'pdf':
            db = bins[1] - bins[0]
            return 1.0/(a.size*db) * n, bins
        elif normed == 'pmf':
            return n/float(sum(n)), bins
    else:
        return n, bins

def histogramSorted(sorteda, bins=10, range=None, normed=False):
    """Builds a histogram, stolen from numpy.histogram(), modified to assume
    sorted input and to allow normed='pdf' or normed='pmf' (prob mass function)"""
    a = np.asarray(sorteda).ravel()
    if not iterable(bins):
        if range is None:
            range = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins, endpoint=False)
    #n = np.sort(a).searchsorted(bins)
    n = a.searchsorted(bins)
    n = np.concatenate([n, [len(a)]]) # this adds a bin that includes overflow points
    n = n[1:]-n[:-1] # subtracts a shifted version of itself
    if normed:
        if normed == 'pdf':
            db = bins[1] - bins[0]
            return 1.0/(a.size*db) * n, bins
        elif normed == 'pmf':
            return n/float(sum(n)), bins
    else:
        return n, bins

def histogram2d(x, y, bins=10, range=None, normed=False):
    """Compute the 2D histogram for a dataset (x,y) given the edges or
    the number of bins.

    Stolen from np.histogram2d() in numpy 1.0
    Modified by mspacek to allow normed='pdf' or normed='pmf' (prob mass function)

    histogram2d(x, y, bins=10, range=None, normed=False) -> H, xedges, yedges

    Compute the 2D histogram from samples x,y.

    Parameters
    ----------
    x,y: 1D data series. Both arrays must have the same length.
    bins: Number of bins -or- [nbin x, nbin y] -or-
         [bin edges] -or- [x bin edges, y bin edges].
    range:  A sequence of lower and upper bin edges (default: [min, max]).
    normed: True or False.

    The histogram array is a count of the number of samples in each
    two dimensional bin.
    Setting normed to 'pdf' returns a density rather than a bin count.
    """
    try:
        N = len(bins)
    except TypeError:
        N = 1
        bins = [bins]
    x = np.asarray(x)
    y = np.asarray(y)
    if range is None:
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
    else:
        xmin, xmax = range[0]
        ymin, ymax = range[1]
    if N == 2:
        if np.isscalar(bins[0]):
            xnbin = bins[0]
            xedges = np.linspace(xmin, xmax, xnbin+1)
        else:
            xedges = np.asarray(bins[0], float)
            xnbin = len(xedges)-1
        if np.isscalar(bins[1]):
            ynbin = bins[1]
            yedges = np.linspace(ymin, ymax, ynbin+1)
        else:
            yedges = np.asarray(bins[1], float)
            ynbin = len(yedges)-1
    elif N == 1:
        ynbin = xnbin = bins[0]
        xedges = np.linspace(xmin, xmax, xnbin+1)
        yedges = np.linspace(ymin, ymax, ynbin+1)
    else:
        yedges = np.asarray(bins, float)
        xedges = yedges.copy()
        ynbin = len(yedges)-1
        xnbin = len(xedges)-1

    dxedges = np.diff(xedges)
    dyedges = np.diff(yedges)

    # Flattened histogram matrix (1D)
    hist = np.zeros((xnbin)*(ynbin), int)

    # Count the number of sample in each bin (1D)
    xbin = np.digitize(x, xedges)
    ybin = np.digitize(y, yedges)

    # Values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    xdecimal = int(-np.log10(dxedges.min()))+6
    ydecimal = int(-np.log10(dyedges.min()))+6
    on_edge_x = np.where(np.around(x, xdecimal) == np.around(xedges[-1], xdecimal))[0]
    on_edge_y = np.where(np.around(y, ydecimal) == np.around(yedges[-1], ydecimal))[0]
    xbin[on_edge_x] -= 1
    ybin[on_edge_y] -= 1
    # Remove the true outliers
    outliers = (xbin==0) | (xbin==xnbin+1) | (ybin==0) | (ybin==ynbin+1)
    xbin = xbin[outliers==False] - 1
    ybin = ybin[outliers==False] - 1

    # Compute the sample indices in the flattened histogram matrix.
    if xnbin >= ynbin:
        xy = ybin*(xnbin) + xbin
    else:
        xy = xbin*(ynbin) + ybin

    # Compute the number of repetitions in xy and assign it to the flattened
    # histogram matrix.
    flatcount = np.bincount(xy)
    indices = np.arange(len(flatcount))
    hist[indices] = flatcount

    # Shape into a proper matrix
    shape = np.sort([xnbin, ynbin])
    histmat = hist.reshape(shape)
    if (shape == (ynbin, xnbin)).all():
        histmat = histmat.T

    if normed:
        if normed == 'pdf':
            diff2 = np.outer(dxedges, dyedges)
            histmat = histmat / diff2 / histmat.sum()
        elif normed == 'pmf':
            histmat = histmat / float(histmat.sum())
        else:
            raise ValueError, 'unknown normed value %s' % normed
    return histmat, xedges, yedges
'''
def histogram2dold(x, y, bins, normed=False):
    """Compute the 2D histogram for a dataset (x,y) given the edges or
    the number of bins. Stolen from np.histogram2d() (numpy 1.0b5), modified to allow
    normed='pdf' or normed='pmf' (prob mass function)

    NOTE: THIS HAS A SERIOUS BUG, IT FAILS TO TAKE THE TRANSPOSE AT ONE POINT, OR SOMETHING,
    DON'T USE!!!!!!!!

    Returns histogram, xedges, yedges.
    The histogram array is a count of the number of samples in each bin.
    The array is oriented such that H[i,j] is the number of samples falling
        into binx[j] and biny[i].
    Data falling outside of the edges are not counted.
    """
    try:
        N = len(bins)
    except TypeError:
        N = 1
        bins = [bins]
    if N == 2:
        if np.isscalar(bins[0]):
            xnbin = bins[0]
            xedges = np.linspace(x.min(), x.max(), xnbin+1)
        else:
            xedges = asarray(bins[0], float)
            xnbin = len(xedges)-1
        if np.isscalar(bins[1]):
            ynbin = bins[1]
            yedges = np.linspace(y.min(), y.max(), ynbin+1)
        else:
            yedges = asarray(bins[1], float)
            ynbin = len(yedges)-1
    elif N == 1:
        ynbin = xnbin = bins[0]
        xedges = np.linspace(x.min(), x.max(), xnbin+1)
        yedges = np.linspace(y.max(), y.min(), ynbin+1)
        xedges[-1] *= 1.0001
        yedges[-1] *= 1.0001
    else:
        yedges = asarray(bins, float)
        xedges = yedges.copy()
        ynbin = len(yedges)-1
        xnbin = len(xedges)-1

    # Flattened histogram matrix (1D)
    hist = np.zeros((xnbin)*(ynbin), int)

    # Count the number of sample in each bin (1D)
    xbin = np.digitize(x,xedges)
    ybin = np.digitize(y,yedges)

    # Remove the outliers
    outliers = (xbin==0) | (xbin==xnbin+1) | (ybin==0) | (ybin == ynbin+1)

    xbin = xbin[outliers==False]
    ybin = ybin[outliers==False]

    # Compute the sample indices in the flattened histogram matrix.
    if xnbin >= ynbin:
        xy = ybin*(xnbin) + xbin
        shift = xnbin + 1
    else:
        xy = xbin*(ynbin) + ybin
        shift = ynbin + 1

    # Compute the number of repetitions in xy and assign it to the flattened
    #  histogram matrix.
    flatcount = np.bincount(xy)
    indices = np.arange(len(flatcount)-shift)
    hist[indices] = flatcount[shift:]

    # Shape into a proper matrix
    histmat = hist.reshape(xnbin, ynbin)

    if normed == 'pdf':
        diff2 = np.outer(np.diff(yedges), np.diff(xedges))
        histmat = histmat / diff2 / histmat.sum()
    elif normed == 'pmf':
        histmat = histmat / float(histmat.sum())
    return histmat, xedges, yedges
'''
def sah(t, y, ts, keep=False):
    """Resample using sample and hold. Returns resampled values at ts given the original points (t,y)
    such that the resampled values are just the most recent value in y (think of a staircase with non-uniform steps).
    Assumes that t is sorted. t and ts arrays should be of the same data type. Contributed by Robert Kern."""
    i = np.searchsorted(t, ts) - 1 # find where ts falls in t, dec so you get indices that point to the most recent value in y
    i = np.where(i < 0, 0, i) # handle the cases where ts is smaller than the first point.
    '''this has an issue of not keeping the original data point where ts == t'''

    ###NOTE: can probably get around having to do this by using searchsorted's new 'side' keyword

    if keep:
        # The following ensures that the original data point is kept when ts == t, doesn't really work if the shortest ISI is less than tres in ts
        di = np.diff(i).nonzero()[0] # find changes in i, nonzero() method returns a tuple, pick the result for the first dim with [0] index
        si = approx(t[1::], ts[di]) # check at those change indices if t ~= ts (ignoring potential floating point representational inaccuracies). If so, inc i at that point so you keep y at that point.
        #print i
        i[di[si]] += 1
        #print i
    return y[i]

def corrcoef(x, y):
    """Returns the correlation coefficient of signals x and y. This just uses np.corrcoef(),
    but converts to floats first, cuz np.corrcoef() seems to have issues with integer signals,
    especially those with zeros in them."""
    #assert len(x) == len(y), 'arrays need to be of equal length'
    x = np.float64(x)
    y = np.float64(y)
    return np.corrcoef(x, y)[0, 1] # pick one of the 2 entries in the correlation coefficient matrix, on the -ve diagonal (er, the one that goes from bottom left to top right, that's what I mean)
    #return ((x * y).mean() - x.mean() * y.mean()) / (x.std() * y.std()) # this works just fine as well, easier to understand too

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
    """Return a string with the binary representation of an integer. If necessary, will append leading zeros
    if result is less than minbits long. Seems like np.binary_repr() is a somewhat faster alternative.
    First 2 lines stolen from Andrew Gaul <andrew@gaul.org> off the web"""
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
    multiplier = np.array(multiplier, ndmin=2).transpose() # convert from list and transpose to a column vector (have to make it 2D to transpose)
    #print multiplier
    x = bin*multiplier
    #print x
    return x.sum(axis=0) # sum over the first dimension (the rows), that way, you're left with only columns in a row vector

def getbinarytable(nbits=8):
    """Generates a 2D binary table containing all possible words for nbits, with bits in the rows and words in the columns (LSB to MSB from top to bottom)"""
    rowlength = 2**nbits
    '''
    x = np.zeros((nbits, 2**nbits)) # init an array
    for bit in range(nbits):
        pattern = [0]*2**bit
        pattern.extend([1]*2**bit)
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        row = pattern*npatterns
        x[bit]=row
    return x
    '''
    '''
    x = np.zeros((nbits, 2**nbits), dtype=np.int8) # init an array
    for bit in range(nbits): # one row at a time
        pattern = np.array(0, dtype=np.int8).repeat(2**bit)
        pattern = np.concatenate((pattern, np.array(1, dtype=np.int8).repeat(2**bit)))
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
        row = np.tile(pattern, [1, npatterns])
        x[bit::,::] = row
    return x
    '''
    # this seems to be the fastest method:
    x = []
    for bit in range(nbits): # one row at a time
        pattern = np.array(0, dtype=np.int8).repeat(2**bit)
        pattern = np.concatenate((pattern, np.array(1, dtype=np.int8).repeat(2**bit)))
        npatterns = rowlength / len(pattern) # == 2**nbits / len(pattern) == 2**nbits / 2**(bit+1) == 2**(nbits-bit-1)
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
    """Finds char in string, returns matching indices. There's gotta be a built-in way to do this somewhere..."""
    assert len(char) == 1
    i = []
    # maybe more efficient to use .find() method on successively smaller slices of string
    for si, s in enumerate(string):
        if s == char:
            i.append(si)
    return i
'''
def shuffle(x):
    """Takes an input list x and returns a shuffled (without replacement) copy. Its only benefit
    over and above random.sample() is that you don't have to pass a second argument len(x)
    every time you use it.
    In NumPy, it's better (and faster) to use np.random.shuffle()"""
    return random.sample(x, len(x))
'''
def shuffle(seq):
    """Takes a sequence and returns a shuffled (without replacement) copy.
    Its only benefit over np.random.shuffle is that it returns a copy instead of shuffling in-place"""
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

def fact(n):
    """Factorial!"""
    assert n.__class__ == int
    assert n >= 0
    if n == 0:
        n = 1 # 0! == 1!
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
    return nPr(n, r) // fact(r)

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
    combs = np.empty(nCr(n, r), dtype=np.object) # stores all combinations, will be a 1D array of arrays
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
    #print code

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
    #print code

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
    maxnsamples = nCr(len(objects), r)
    if nsamples == None:
        nsamples = maxnsamples # return all possible combinations
    if nsamples > maxnsamples: # make sure we're not being asked for more than the maximum possible number of unique samples
        raise ValueError, 'requested unique nsamples (%d) is larger than len(objects) choose r (%d C %d == %d)' % (nsamples, len(objects), r, maxnsamples)
    # I've set the criteria for generating a table to be never, cuz generating the table and then sampling it almost always takes longer (at least for maxnsamples as high as 325, say) than just picking combs at random and making sure they're unique
    if maxnsamples < 0: # generate a table of all possible combinations, and then just pick nsamples from it without replacement
        table = combs(objects, r)
        samples = random.sample(table, nsamples)
    elif r == 1: # we're just choosing one item from objects at a time
        samples = random.sample(objects, nsamples)
    else: # the number of possible combs is inconveniently large to completely tabulate, pick some combinations at random and make sure each comb is unique
        samples = []
        samplei = 0
        while samplei < nsamples:
            sample = random.sample(objects, r) # choose r objects at random
            sample.sort() # sort for sake of comparison with other samples, important cuz this removes any differences due to permuatations (as opposed to combs)
            if sample not in samples: # make sure they're not the same set of objects as any previous set in samples
                samples.append(sample) # add it to the list of samples
                samplei += 1
    return samples

'''
# this f'n isn't really needed, just use objlist.sort(key=lambda obj: obj.attrib)
def sortby(objs, attrib, cmp=None, reverse=False):
    """Returns objects list sorted according to the specified object attribute.
    attrib should be passed as a string"""
    objs.sort(key=lambda obj: obj.__getattribute__(attrib), cmp=cmp, reverse=reverse) # sort in-place
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
    for dataslice in data: # this for loop isn't such a bad thing cuz the massive add step inside the loop is the limiting factor
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

def normalize(seq):
    """Normalizes a sequence, returning zeros if sum(seq) == 0"""
    a = np.asarray(seq)
    if a.sum() == 0: # numpy doesn't raise ZeroDivisionErrors for some reason
        return np.zeros(a.shape) # just return zeros
    else:
        return a / float(a.sum()) # return it normalized

def ensurenormed(p, atol=1e-8):
    """Ensures p is normalized. Returns p unchanged if it's already normalized,
    otherwise, prints a warning and returns it normalized. atol is how close to 1.0
    p.sum() needs to be"""
    p = np.asarray(p)
    psum = p.sum()
    if not approx(psum, 1.0, atol=atol): # make sure the probs sum to 1
        print 'ps don''t sum to 1, they sum to %f instead, normalizing for you' % psum
        p /= float(psum)
    return p

def logy(x, base=10):
    """Performs log of x with specified base"""
    return np.log(x) / np.log(base)

def log_no_sing(x, subval=0.0, base=np.e):
    """Performs log on array x, ignoring any zeros in x to avoid singularities,
    and returning subval in their place in the result"""
    x = np.asarray(x)
    singi = x==0 # find the singularities
    x[singi] = 1 # replace 'em with 1s, or anything else that's safe to take the log of
    result = logy(x, base=base) # now it's safe to take the log
    result[singi] = subval # substitute the result where the singularities were with the substitution value
    return result

def log10_no_sing(x, subval=0.0):
    """Performs log10 on x, ignoring singularities"""
    return log_no_sing(x, subval=subval, base=10)

def log2_no_sing(x, subval=0.0):
    """Performs log2 on x, ignoring singularities"""
    return log_no_sing(x, subval=subval, base=2)

def entropy(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob values in p"""
    p = ensurenormed(p)
    return -(p * np.log2(p)).sum()

def entropy_no_sing(p):
    """Returns the entropy (in bits) of the prob distribution described by the prob values in p
    Ignore singularities in p (assumes their contribution to entropy is zero)"""
    p = ensurenormed(p)
    return -(p * log2_no_sing(p, subval=0.0)).sum()

def MI(XY):
    """Given the joint PDF of two variables, returns the mutual information (in bits) between the two
    I = sum_X sum_Y P(x, y) * log2( P(x, y) / (P(x) * P(y)) )
    where P(x) and P(y) are the marginal distributions taken from the joint
    Is this slow? Needs optimization? Already implemented in scipy???????????????"""
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
    # build up joint pdf of all the possible N words, and the two possible N+1th values (0 and 1)
    xedges = np.arange(2**N+1) # values 0 to 2**N - 1, plus 2**N which is needed as the rightmost bin edge for histogram2d (annoying)
    yedges = np.arange(2**M+1)
    bins = [xedges, yedges]
    jpdf, xedgesout, yedgesout = histogram2d(Nintcodes, Mintcodes, bins, normed='pmf') # generate joint pdf
    #print 'jpdf\n', jpdf.__repr__()
    #print 'jpdf.sum()', jpdf.sum()
    assert (np.float64(xedges) == xedgesout).all() # sanity check
    assert (np.float64(yedges) == yedgesout).all()
    # pdf of N cells
    #Npdf, Nedges = histogram(Nintcodes, bins=range(2**N), normed='pmf')
    #print 'first 100 Npdf\n', Npdf[:100].__repr__()
    # pdf of M cells
    #Mpdf, Medges = histogram(Mintcodes, bins=range(2**M), normed='pmf')
    #print 'first 100 Mpdf\n', Mpdf[:100].__repr__()
    marginalMpdf = jpdf.sum(axis=0)
    #assert approx(Mpdf, marginalMpdf).all() # make sure what you get from the joint is what you get when just building up the pdf straight up on its own
    I = MI(jpdf)
    IdivS = I / entropy(marginalMpdf) # return mutual info as fraction of entropy in M group of cells
    if verbose:
        print 'nids', nids
        print 'mids', mids
        #print 'Mpdf', Mpdf
        #print 'entropy(Mpdf)', entropy(Mpdf)
        print 'marginal Mpdf', marginalMpdf
        print 'entropy(marginal Mpdf)', entropy(marginalMpdf)
        print 'I', I
        print 'I/entropy', IdivS
    if not 0.0 <= IdivS <= 1.0+1e-10:
        import pdb; pdb.set_trace()
        print 'IdivS is out of range'
        print 'IdivS is %.16f' % IdivS

    return dictattr(I=I, IdivS=IdivS)

def DKL(p, q):
    """Kullback-Leibler divergence from true probability distribution p to arbitrary distribution q"""
    assert len(p) == len(q)
    p = ensurenormed(p)
    q = ensurenormed(q)
    return sum([ pi * np.log2(pi/float(qi)) for pi, qi in zip(p, q) if pi != 0 and qi != 0 ] ) # avoid singularities

def DJS(p, q):
    """Jensen-Shannon divergence, a symmetric measure of divergence between distributions p and q"""
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

def eof(f):
    """Return whether file pointer is a end of file"""
    orig = f.tell()
    f.seek(0, 2) # seek 0 bytes from end
    return f.tell() == orig
