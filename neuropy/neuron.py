"""Defines the Neuron class and all of its support classes"""

from __future__ import division

import os
import StringIO
import time
import hashlib

import numpy as np
import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import util # .pyx file

import pylab as pl
import matplotlib as mpl
from pylab import get_current_fig_manager as gcfm

import core
from core import rstrip, getargstr, iterable, toiter, tolist, intround
from core import mean_accum, lastcmd, RevCorrWindow
from core import PTCSNeuronRecord, SPKNeuronRecord
from dimstimskeletal import Movie


class BaseNeuron(object):
    """A neuron's spike data spans all the experiments within a recording.
    If different recordings have sorts with the same name, you can assume they were
    extracted in the same spike sorting session, and that therefore the neuron ids
    are the same"""
    def __init__(self, path, sort=None):
        self.level = 5 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.path = path
        self.sort = sort

    def get_name(self):
        fname = os.path.split(self.path)[-1] # pathless
        return os.path.splitext(fname)[0] # extensionless

    name = property(get_name)

    id = property(lambda self: self.record.nid)
    descr = property(lambda self: self.record.descr)
    #clusterscore = property(lambda self: self.record.clusterscore)
    pos = property(lambda self: (self.record.xpos, self.record.ypos))
    nchans = property(lambda self: self.record.nchans)
    chans = property(lambda self: self.record.chans)
    maxchan = property(lambda self: self.record.maxchan)
    wavedata = property(lambda self: self.record.wavedata)
    wavestd = property(lambda self: self.record.wavestd)
    nspikes = property(lambda self: self.record.nspikes)
    spikes = property(lambda self: self.record.spikes)
    
    '''
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.sort.writetree(string)

    def load(self, f=None):
        #treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        #self.writetree(treestr+'\n')
        #print(treestr) # print string to tree hierarchy and screen
        ext = os.path.splitext(self.path)[1]
        if ext == '.ptcs':
            self.loadptcs(f)
        elif ext == '.spk':
            self.loadspk()
    '''
    def loadptcs(self, f, header):
        """Read in the next neuron record in an open .ptcs file"""
        nrec = PTCSNeuronRecord(header)
        nrec.read(f)
        self.record = nrec
        self.post_load()

    def loadspk(self):
        """Read in a .spk file containing purely spike times"""
        nrec = SPKNeuronRecord(self.path)
        nrec.read()
        self.record = nrec
        self.post_load()

    def post_load(self):
        if self.nspikes == 0:
            raise RuntimeError('neuron %d in %s has no spikes' % (self.id, self.path))
        self.trange = self.spikes[0], self.spikes[-1]
        self.dt = self.trange[1] - self.trange[0]
        self.dtsec = self.dt / 1e6
        self.dtmin = self.dtsec / 60
        self.dthour = self.dtmin / 60


class NeuronBasics(object):
    """Mix-in class that defines basic Neuron methods"""
    def cut(self, *args):
        """Returns a view of the Neuron's spike times where tstart <= spikes <= tend
        *args can be: nothing (returns all spikes), None, tstart, or (tstart, tend)
        None and 0 for tstart are shorthand for 'from first spike'
        None and -1 for tend are shorthand for 'to last spike'"""
        print("WARNING: Neuron.cut should be deprecated, use spikes.searchsorted "
              "directly instead")
        tstart, tend = self._parse_cutargs(*args)
        '''
        # this is what we're trying to do:
        return self.spikes[ (self.spikes >= tstart) & (self.spikes <= tend) ]
        self.searchsorted(values) method does it faster. It returns an index where
        the value would fit in self. The index is such that self[index-1] < value <= self[index].
        In this formula self[self.size]=inf and self[-1]= -inf
        Another possibility could be to use a masked array instead?
        '''
        lo, hi = self.spikes.searchsorted([tstart, tend]) # returns indices where tstart and tend would fit in spikes
        if tend == self.spikes[min(hi, self.nspikes-1)]: # if tend matches a spike (protect from going out of index bounds when checking)
            hi += 1 # inc to include a spike if it happens to exactly equal tend. This gives us end inclusion
            hi = min(hi, self.nspikes) # limit hi to max slice index (==max value index + 1)
        cutspikes = self.spikes[lo:hi] # slice it
        return cutspikes

    def cutrel(self, *args):
        """Cuts Neuron spike times and returns them relative to tstart"""
        print("WARNING: Neuron.cutrel should be deprecated, use spikes.searchsorted "
              "directly instead")
        tstart, tend = self._parse_cutargs(*args)
        cutspikes = self.cut(tstart, tend)
        # let's keep all the returned spikes as integers, us is more than accurate enough
        tstart = np.int64(round(tstart))
        return cutspikes - tstart

    def _parse_cutargs(self, *args):
        """Takes set of args as passed to cut or cutrel and returns (tstart, tend)"""
        print("WARNING: Neuron._parse_cutargs should be deprecated, use spikes.searchsorted "
              "directly instead")
        tstart = None
        tend = None
        if len(args) == 0: # passed nothing
            pass
        elif len(args) == 1: # passed None, or just tstart, or a (tstart, tend) sequence
            if not iterable(args[0]):
                if args[0] == None:
                    pass
                else: # just tstart was passed
                    tstart = args[0]
            elif len(args[0]) == 2: # it's a sequence
                tstart = args[0][0]
                tend = args[0][1]
            else:
                raise ValueError('sequence is too long')
        elif len(args) == 2: # passed tstart and tend as separate args
            tstart = args[0]
            tend = args[1]
        else:
            raise ValueError('too many arguments')
        if tstart in [None, 0]:
            # shorthand for "from first spike" - would be problematic if a spike existed at t=0
            tstart = self.spikes[0]
        if tend in [None, -1]:
            # shorthand for "to last spike" - would be problematic if a spike existed at t=-1
            tend = self.spikes[-1]
        return (tstart, tend)

    def copy(self):
        """Returns a copy of the Neuron"""
        return copy(self)

    def isi(self, trange=None):
        """Returns the inter-spike intervals of the Neuron's spike train in trange"""
        return diff(self.cut(trange))

    def iisii(self, trange=None):
        """Returns the inter-ISI intervals (the differences between consecutive ISIs)
        of the Neuron's spike train in trange.
        See delta def'n in: 2002 Segev, et al - Long term behavior of lithographically..."""
        return diff(self.isi(trange))


class XCorr(object):
    def __init__(self, n0, nid1, trange=50):
        """Cross-correlation object. n0 is a Neuron, nid1 is a nid, +/- trange is in ms"""
        self.n0 = n0
        self.n1 = n0.sort.n[nid1]
        self.autocorr = self.n0.id == self.n1.id
        trange = abs(trange) * 1000 # convert to us
        self.trange = np.array([-trange, trange]) # convert to a +/- array, in us

    def calc(self):
        t0 = time.time()
        dts = util.xcorr(self.n0.spikes, self.n1.spikes, trange=self.trange) # in us
        print('xcorr calc took %.3f sec' % (time.time()-t0))
        self.dts = np.array(dts)
        if self.autocorr:
            self.dts = self.dts[self.dts != 0] # remove 0s for autocorr
        return self

    def plot(self, nbins=None, rate=False, figsize=(7.5, 6.5)):
        """style can be 'rate', but defaults to count"""
        if nbins == None:
            nbins = intround(np.sqrt(len(self.dts))) # good heuristic
        dts = self.dts / 1000 # in ms, converts to float64 array
        trange = self.trange / 1000 # in ms, converts to float64 array
        nbins = max(20, nbins) # enforce min nbins
        nbins = min(200, nbins) # enforce max nbins
        t = np.linspace(start=trange[0], stop=trange[1], num=nbins, endpoint=True)
        n = np.histogram(dts, bins=t, density=False)[0]
        binwidth = t[1] - t[0] # all should be equal width
        if rate: # normalize by binwidth and convert to float:
            n = n / float(binwidth)
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        a.bar(left=t[:-1], height=n, width=binwidth) # omit last right edge in t
        a.set_xlim(t[0], t[-1])
        a.set_xlabel('ISI (ms)')
        if rate:
            a.set_ylabel('spike rate (Hz)')
        else:
            a.set_ylabel('count')
        #a.set_title('n%d spikes relative to n%d spikes' % (self.n1.id, self.n0.id))
        title = lastcmd() + ', binwidth: %.2f ms' % binwidth
        a.set_title(title)
        gcfm().window.setWindowTitle(title)
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        return self


class NeuronXCorr(object):
    """Mix-in class that defines the xcorr Neuron method"""
    def xcorr(self, nid, trange=50):
        """Return cross-correlation of self and other nid over +/- trange in ms"""
        xco = XCorr(self, nid, trange=trange)
        xco.calc()
        return xco


class BinaryCode(object):
    """Quantize a spike train, cut according to tranges in us, shifted by shift ms,
    into a binary signal with values CODEVALS and time resolution CODETRES in us.
    CODEPHASE specifies where to start the codetrain in time, relative to the nearest
    multiple of CODETRES before each trange. Phase is in degrees of a single bin period. -ve
    phase is leading (codetrain starts earlier in time), +ve is lagging (codetrain starts
    later in time)"""
    def __init__(self, spikes=None, tranges=None, shift=0):
        uns = get_ipython().user_ns
        self.kind = 'binary'
        self.spikes = spikes
        if tranges == None:
            tranges = [ self.neuron.trange ] # just the one full trange
        self.tranges = tranges
        self.shift = shift
        self.codevals = uns['CODEVALS']
        self.tres = uns['CODETRES']
        self.phase = uns['CODEPHASE']

    def get_hash(self):
        """Return characteristic hash of everything that affects calc output"""
        h = hashlib.md5()
        for thing in [self.kind, self.spikes, self.tranges, self.shift,
                      self.codevals, self.tres, self.phase]:
            if type(thing) in [int, list]:
                thing = np.asarray(thing) # ints and lists can't be hashed
            h.update(thing)
        return h.hexdigest()

    hash = property(get_hash)

    def calc(self):
        """.t and .s attribs are commented out to save substantial memory"""
        self.t = [] # bin times
        #self.s = [] # relevant spike times, potentially shifted by self.shift
        self.c = [] # code values for each bin
        shift = intround(self.shift * 1000) # convert self.shift in ms to int us
        for trange in self.tranges:
            # make the start of the timepoints be an even multiple of self.tres. Round down
            # to the nearest multiple. This way, timepoints will line up for different code
            # objects
            # left edge of first code bin:
            tstart = trange[0] - (trange[0] % self.tres)
            if self.phase: # add phase offset relative to tstart
                tstart += self.phase / 360.0 * self.tres
            tstart = intround(tstart) # keep it int
            tend = intround(trange[1]) # ditto
            # t sequence demarcates left bin edges, add extra tres to end to make t end
            # inclusive, keep 'em in us integers:
            t = np.arange(tstart, tend+self.tres, self.tres) # should come out as int64
            # get relevant spike times s, cut over originally specified trange, not from
            # start to end of newly generated code bin timepoints:
            lo, hi = self.spikes.searchsorted(trange)
            s = self.spikes[lo:hi] + shift
            c = np.empty(len(t), dtype=np.int8) # init binary code array
            c[:] = self.codevals[0] # init code to low value
            # searchsorted returns indices where s fits into t. Sometimes more than one
            # spike will fit into the same time bin, which means searchsorted will return
            # multiple occurences of the same index. You can set c at these indices to 1 a
            # multiple number of times, or more efficiently, do an np.unique on it to
            # only set each index to 1 once.
            # dec index by 1 so that you get indices that point to the most recent bin edge.
            # For each bin that has at least 1 spike in it, set its value to high:
            c[np.unique(t.searchsorted(s)) - 1] = self.codevals[1]
            self.t.append(t)
            #self.s.append(s)
            self.c.append(c)
        # horizontally concatenate results from each trange:
        self.t = np.hstack(self.t)
        #self.s = np.hstack(self.s)
        self.c = np.hstack(self.c)
        del self.spikes # no need for spikes any more, save memory

    def plot(self):
        """Plot some kind of long grid of white and black elements?"""
        raise NotImplementedError
        

class NeuronCode(object):
    """Mix-in class that defines the spike code related Neuron methods"""
    def code(self, tranges=None, shift=0):
        """Returns an existing Code object, or creates and calcs a new one if necessary"""
        try:
            self._codes
        except AttributeError: # self._codes doesn't exist yet
            self._codes = {} # create a dict that'll hold Code objects for this Neuron
        kind = get_ipython().user_ns['CODEKIND']
        if kind == 'binary': # init a new BinaryCode object
            co = BinaryCode(self.spikes, tranges, shift)
            co_hash = co.hash
        else:
            raise ValueError('Unknown kind: %r' % kind)
        try:
            co = self._codes[co_hash] # there's already one that's been calc'd
        except KeyError:
            co.calc()
            self._codes[co_hash] = co # add it to the Code dict
        return co
    code.__doc__ += '\n\nbinary:\n' + BinaryCode.__doc__


class BaseRate(object):
    """Abstract firing rate class"""
    def __init__(self, neuron=None, trange=None):
        self.neuron = neuron
        if trange == None:
            self.trange = self.neuron.trange
        else:
            self.trange = trange

    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their r and t attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd] for key in ['r', 't', 'rawr', 'rawt'] if d.has_key(key) ]
        if type(self) == type(other) and selfd == otherd:
            return True
        else:
            return False

    def plot(self):
        pl.figure()
        pl.plot(self.t, self.r)
        # diagnostic for comparing interpolated to raw nisi rate:
        #pl.plot(self.rawt, self.rawr, 'r+')
        pl.xlabel('t')
        pl.ylabel('spike rate')


class BinRate(BaseRate):
    """Uses simple binning to calculate firing rate"""
    def __init__(self, neuron=None, trange=None, tres=100000):
        super(BinRate, self).__init__(neuron=neuron, trange=trange)
        self.kind = 'bin'
        self.tres = tres

    def calc(self):
        # make the start of the timepoints be an even multiple of self.tres. Round down to
        # the nearest multiple. Do the same for the end of the timepoints. This way,
        # timepoints will line up for different code objects
        tstart = self.trange[0] - (self.trange[0] % self.tres)
        tend   = self.trange[1] - (self.trange[1] % self.tres)
        # t sequence demarcates left bin edges, add tres to trange[1] to make t end inclusive:
        t = np.arange(tstart, tend+self.tres, self.tres)
        s = self.neuron.cut(self.trange) # spike times
        self.r, self.t = np.histogram(self.neuron.spikes, bins=t)
        self.r = self.r / float(self.tres) * 1000000 # spikes/sec

    def plot(self):
        super(BinRate, self).plot()
        pl.title('neuron %d - binned spike rate' % self.neuron.id)


class nISIRate(BaseRate):
    """Uses nisi inter spike intervals to calculate firing rate, with a fixed number of spikes nisi+1 per bin. nisi == 1 is the ISI rate"""
    def __init__(self, neuron=None, trange=None, nisi=3, tres=50000, interp='sah'):
        super(nISIRate, self).__init__(neuron=neuron, trange=trange)
        self.kind = 'nisi'
        self.nisi = nisi
        self.tres = tres
        self.interp = interp

    def calc(self):
        s = self.neuron.cut(self.trange) # spike times
        #n0 = self.n-1 # 0-based n
        # compare s to a shifted version of itself:
        diff = s[self.nisi::] - s[:-self.nisi:] # (n0 to end, single steps) - (beginning to n0 from end, single steps)
        r = float(self.nisi+1) / diff * 1000000 # spikes/sec
        t = s[self.nisi::] # for the corresponding timepoints, pick the spike time at the end of the group of n spikes to keep it causal
        #self.rawr = r
        #self.rawt = t
        if not self.interp:
            self.t = t
            self.r = r
            return
        # make the start of our interpolated timepoints be an even multiple of self.tres. Round down to the nearest multiple. This way, the timepoints will line up, even if different Rates have different starting points, like neuron.rate() vs experiment.rate()
        tstart = t[0] - (t[0] % self.tres)
        # should we have tend = t[-1] + self.tres ?
        self.t = np.arange(tstart, t[-1], self.tres) # new set of timepoints to interpolate over
        if self.interp == 'sah':
            self.r = sah(t, r, self.t, keep=False)
        elif self.interp == 'linear':
            f = sig.interpolate.interp1d(t, r, kind='linear') # returns an interpolation f'n
            self.r = f(self.t) # interpolate over the new timepoints
            # or maybe try sig.resample() instead
        else:
            raise ValueError('unknown interpolation method: %s' % self.interp)

    def plot(self):
        super(nISIRate, self).plot()
        pl.title('neuron %d - %d-inter-spike-interval spike rate, %s interpolation'
                     % (self.neuron.id, self.nisi, self.interp))


class WnISIRate(BaseRate):
    """Uses a weighted sum of various n inter spike intervals to calculate rate"""
    pass
    #def plot(self)


class IDPRate(BaseRate):
    """Instantaneous discharge probability. See Pauluis and Baker, 2000"""
    def __init__(self, neuron=None, IDP_a=4, tres=50000):
        super(IDPRate, self).__init__(neuron=neuron, trange=trange)
        self.kind = 'idp'
        self.tres = tres

    def calc(self):

        # Step 1: find sudden changes in ISI distribution
        Is = self.neuron.isi() # intervals
        # y=gdtrc(a,b,x) returns the integral from x to infinity of the gamma probability density function.  SEE gdtr, gdtri
        b = IDP_a
        increaseOnsets = [] # stores the interval indices of significant sudden increases in ISI rate
        decreaseOnsets = [] # stores the interval indices of significant sudden decreases in ISI rate
        for (n,I) in enumerate(Is):
            try:
                IDP_mu = Is(n+1) # the next interval
                a = IDP_a / float(IDP_mu)
                p1 = gdtrc(a,b,I) # integrate from I to infinity
                IDP_mu = Is(n+2) # the interval after that
                a = IDP_a / float(IDP_mu)
                p2 = gdtrc(a,b,I) # integrate from I to infinity
                if p1 < pthresh and p2 < pthresh: # then this nth interval was a sudden increase in ISI firing rate
                    inreaseOnsets.append(n)
            except IndexError:
                pass
            try:
                IDP_mu = Is(n-1) # the previous interval
                a = IDP_a / float(IDP_mu)
                p1 = gdtrc(a,b,I) # integrate from I to infinity
                IDP_mu = Is(n-2) # the interval before that
                a = IDP_a / float(IDP_mu)
                p2 = gdtrc(a,b,I) # integrate from I to infinity
                if p1 < pthresh and p2 < pthresh: # then this nth interval was a sudden increase in ISI firing rate
                    dereaseOnsets.append(n)
            except IndexError:
                pass

        print 'IDPRATE IS INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        raise NotImplementedError

        # Step 2: extend back half an ISI from sudden rate increase, and forward half an ISI from sudden rate decrease

        # Step 3: sample our new ISI rate with extended intervals with sufficient resolution
        isi = self.neuron.rate(kind='nisi', nisi=1, interp=None) # ISI rate object with with no interpolation
        isi.r
        isi.t
        isi = sah(isi)

        # use np.piecewise() ?
    pass
    #def plot(self)


class GaussRate(BaseRate):
    """Uses a sliding Gaussian window to calculate firing rate"""
    def __init__(self, neuron=None, trange=None, width=200000):
        super(GaussRate, self).__init__(neuron=neuron, trange=trange)
        self.kind = 'gauss'
        self.width = width

    def calc(self):
        # use scipy.signal.gaussian, use just left half to keep it causal
        raise NotImplementedError

    def plot(self):
        super(GaussRate, self).plot()
        pl.title('Gaussian sliding window spike rate')


class RectRate(BaseRate):
    """Uses a sliding rectangular window to calculate firing rate"""
    def __init__(self, neuron=None, trange=None, width=200000):
        super(RectRate, self).__init__(neuron=neuron, trange=trange)
        self.kind = 'rect'
        self.width = width

    def calc(self):
        # use scipy.signal.boxcar
        raise NotImplementedError

    def plot(self):
        super(RectRate, self).plot()
        pl.title('rectangular sliding window spike rate')


class RatePDF(object):
    """Firing rate probability distribution function"""
    def __init__(self, neuron=None, rrange=(0, 200), nbins=100, scale='log', density=True,
                 **kwargs):
        # rrange == rate range, ie limits of pdf x axis; nbins == number of rate bins
        self.neuron = neuron
        self.rrange = rrange
        self.nbins = nbins
        self.scale = scale
        self.density = density
        self.rate = neuron.rate(**kwargs) # returns a Rate object of some kind

    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their n and r and logrrange attribs, if they exist, to prevent comparing
        # them below, since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd] for key in ['n', 'r', 'logrrange']
          if d.has_key(key) ]
        if type(self) == type(other) and selfd == otherd:
            return True
        else:
            return False

    def calc(self):
        if self.scale == 'log':
            safe = list(self.rrange) # convert from immutable tuple to mutable list
            # prevent from taking log(0):
            for i, r in enumerate(safe):
                if r == 0:
                    safe[i] = 0.1 # set to 0.1 Hz
            self.logrrange = log10(tuple(safe)) # convert back to tuple, is now safe to take log
            # r sequence demarcates left rate bin edges
            r = np.logspace(start=self.logrrange[0], stop=self.logrrange[1], num=self.nbins, endpoint=True, base=10.0)
        elif self.scale == 'linear':
            r = np.linspace(start=self.rrange[0], stop=self.rrange[1], num=self.nbins, endpoint=True)
        else:
            raise ValueError('unknown scale: %r' % scale)
        self.n, self.r = np.histogram(self.rate.r, bins=r, density=self.density)

    def plot(self):
        pl.figure()
        if self.scale == 'log':
            barwidth = list(diff(self.r)) # each bar will have a different width, convert to list so you can append
            # need to add one more entry to barwidth to the end to get nbins of them:
            #barwidth.append(barwidth[-1]) # not exactly correct
            logbinwidth = (self.logrrange[1]-self.logrrange[0]) / float(self.nbins)
            barwidth.append(10**(self.logrrange[1]+logbinwidth) - self.r[-1]) # should be exactly correct
        elif self.scale == 'linear':
            barwidth = (self.rrange[1]-self.rrange[0]) / float(self.nbins)
        else:
            raise ValueError('unknown scale: %r' % scale)
        pl.bar(left=self.r, height=self.n, width=barwidth)
        pl.axes().set_xscale(self.scale, basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
        title('neuron %d - %s spike rate PDF' % (self.neuron.id, self.rate.kind))
        if self.density:
            pl.ylabel('probability density')
        else:
            pl.ylabel('count')
        pl.xlabel('spike rate')


class NeuronRate(object):
    """Mix-in class that defines the spike rate related Neuron methods"""
    def rate(self, kind='nisi', **kwargs):
        """Returns an existing Rate object, or creates a new one if necessary"""
        try:
            self._rates
        except AttributeError: # self._rates doesn't exist yet
            self._rates = [] # create a list that'll hold Rate objects
        if kind == 'bin':
            ro = BinRate(neuron=self, **kwargs) # init a new BinRate object
        elif kind == 'nisi':
            ro = nISIRate(neuron=self, **kwargs) # init a new nISIRate object
        elif kind == 'wnisi':
            ro = WnISIRate(neuron=self, **kwargs) # init a new WnISIRate object
        elif kind == 'idp':
            ro = IDPRate(neuron=self, **kwargs) # init a new IDPRate object
        elif kind == 'gauss':
            ro = GaussRate(neuron=self, **kwargs) # init a new GaussRate object
        elif kind == 'rect':
            ro = RectRate(neuron=self, **kwargs) # init a new RectRate object
        else:
            raise ValueError('unknown kind: %r' % self.kind)
        for rate in self._rates:
            if ro == rate: # need to define special == method for class Rate()
                return rate # returns the first Rate object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._rates
        ro.calc() # no matching Rate was found, calculate it
        self._rates.append(ro) # add it to the Rate object list
        return ro
    rate.__doc__ += '\n\n**kwargs:'
    _rateargs = '\nbin: '+getargstr(BinRate.__init__)
    _rateargs += '\nnisi: '+getargstr(nISIRate.__init__)
    _rateargs += '\nwnisi: '+getargstr(WnISIRate.__init__)
    _rateargs += '\nidp: '+getargstr(IDPRate.__init__)
    _rateargs += '\ngauss: '+getargstr(GaussRate.__init__)
    _rateargs += '\nrect: '+getargstr(RectRate.__init__)
    rate.__doc__ += _rateargs

    def ratepdf(self, **kwargs):
        """Returns an existing RatePDF object, or creates a new one if necessary"""
        try:
            self._ratepdfs
        except AttributeError: # self._ratepdfs doesn't exist yet
            self._ratepdfs = [] # create a list that'll hold RatePDF objects
        rpdf = RatePDF(neuron=self, **kwargs)
        for ratepdf in self._ratepdfs:
            if rpdf == ratepdf: # need to define special == method for class RatePDF()
                return ratepdf # returns the first RatePDF object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._ratepdfs
        rpdf.calc() # no matching RatePDF object was found, calculate it
        self._ratepdfs.append(rpdf) # add it to the RatePDF object list
        return rpdf
    ratepdf.__doc__ += '\n\n**kwargs:'
    ratepdf.__doc__ += '\nRatePDF: '+getargstr(RatePDF.__init__)
    ratepdf.__doc__ += '\nrate: '+getargstr(rate)
    ratepdf.__doc__ += _rateargs


class RevCorr(object):
    """Base class for doing reverse correlation of spikes to stimulus"""
    def __init__(self, neuron=None, experiment=None, trange=None, nt=10):
        self.neuron = neuron
        self.experiment = experiment
        if trange == None:
            self.trange = self.experiment.trange
        else:
            self.trange = trange
        assert type(self.experiment.e) == Movie
        self.movie = self.experiment.e
        try:
            self.movie.frames # check if movie frames have been loaded from file
        except AttributeError:
            # Load as 3D array instead of as a list of 2D arrays, more convenient for
            # analysis, although will cause memory problems for really big (>1GB movies).
            # Don't flip the movie frames vertically for OpenGL's bottom left origin, since
            # we aren't using OpenGL for analysis:
            self.movie.load(asarray=True, flip=False)
        self.nt = nt # number of revcorr timepoints
        self.tis = range(0, nt, 1) # revcorr timepoint indices
        # revcorr timepoint values, stored in a list, not an array. Bad behaviour happens
        # during __eq__ below if attribs are numpy arrays cuz comparing numpy arrays returns
        # an array of booleans, not just a simple boolean:
        self.ts = [ intround(ti * self.movie.dynamic.sweepSec * 1000) for ti in self.tis ]
        self.ndinperframe = intround(self.movie.dynamic.sweepSec * 1000000 /
                                     self.experiment.REFRESHTIME)
        #self.movie.frames = np.asarray(self.movie.frames)
        self.width = self.movie.frames.shape[-1] # (nframes, height, width)
        self.height = self.movie.frames.shape[-2]
        self.done = False # hasn't yet successfully completed its calc() method

    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # delete their rcdini and rf attribs, if they exist, to prevent comparing them below,
        # since those attribs may not have yet been calculated:
        [ d.__delitem__(key) for d in [selfd, otherd]
            for key in ['rcdini', 'rf', 'done'] if d.has_key(key) ]
        if type(self) == type(other) and selfd == otherd:
            return True
        else:
            return False

    def calc(self):
        """General calc step that has to be performed for all kinds of reverse correlations"""
        spikes = self.neuron.cut(self.trange)
        # revcorr dini. Find where the spike times fall in the din, dec so you get indices
        # that point to the most recent din value for each spike:
        self.rcdini = self.experiment.din[:, 0].searchsorted(spikes) - 1
        #self.din = self.experiment.din[rcdini, 1] # get the din (frame indices) at the rcdini

    def plot(self, interp='nearest', normed=True, title='RevCorrWindow', scale=2.0):
        """Plots the spatiotemporal RF as bitmaps in a wx.Frame"""
        # create a copy to manipulate for display purposes, (nt, width, height):
        rf = self.rf.copy()
        if normed: # normalize across the timepoints for this RevCorr
            norm = mpl.colors.normalize(vmin=rf.min(), vmax=rf.max(), clip=True)
            rf = norm(rf) # normalize the rf the same way across all timepoints
        else: # don't normalize across timepoints, leave each one to autoscale
            for ti in range(self.nt):
                norm = mpl.colors.normalize(vmin=None, vmax=None, clip=True)
                rf[ti] = norm(rf[ti]) # normalize the rf separately at each timepoint
        rf *= 255 # scale up to 8 bit values
        rf = rf.round().astype(np.uint8) # downcast from float to uint8
        win = RevCorrWindow(title=title, rfs=[rf], nids=[self.neuron.id], 
                            ts=self.ts, scale=scale)
        win.show()
        return win # necessary in IPython
    '''
    def oldplot(self, interp='nearest', normed=True):
        """Plots the RFs as images, returns all the image objects"""
        mpl.rcParams['interactive'] = False # there's gotta be a figure or frame object attribute you can set instead of having to change the global rcParams
        #mpl.rcParams['toolbar'] = None # turn off toolbars for this figure. There's gotta be a more OO way...
        nt = self.nt
        figure(figsize=(nt*0.9, 1.5)) # in inches
        gcfm().frame.GetStatusBar().Hide()
        gcfm().frame.GetToolBar().Hide()
        gcfm().frame.Fit()
        #frame()
        # looks like the minimum frame width in Windows is 112 pixels, minimum height is 27 pix
        # with no toolbar, menubar, statusbar, a
        # (112, 112) frame has usable area (104, 85), which means (8, 27) are used by border and caption regions
        frameedgepix = (8, 27)
        framewidthpix = 100*nt+frameedgepix[0]
        frameheightpix = 20+100+frameedgepix[1]
        gcfm().frame.SetSize((framewidthpix, frameheightpix)) # hack
        # look up wxSizers to see how to get around stupid mpl problems and size stuff properly
        #gcfm().frame.Fit()
        #gcf().set_figsize_inches([8, 0.8])
        #as = [] # stores axes handles
        ias = [] # stores image axes handles
        hspace = 0.01
        hborder = 0.15
        vborder = 0.03
        width = (1.0 - hspace*(nt-1) - vborder*2) / nt
        height = 1 - hborder*2
        bottom = hborder
        if normed:
            vmin, vmax = self.rf.min(), self.rf.max()
        else:
            vmin, vmax = None, None
        for ti, t in zip(self.tis, self.t):
            left = (width + hspace)*ti + vborder
            a = axes([left, bottom, width, height])
            # should there be something in here to ensure the image uses up the whole axes space with no borders?
            ia = imshow(self.rf[ti], vmin=vmin, vmax=vmax, interpolation=interp)
            ia.axes.axison = False # turns off x and y axis
            ias.append(ia)
            gcf().text(left, 1-hborder/2.0, '%dms' % t, horizontalalignment='center', verticalalignment='center') # the -0.06 is a hack, verticalalignment is off, cuz of hidden statusbar/toolbar?
        mpl.rcParams['interactive'] = True
        pl.show()
        #mpl.rcParams['toolbar'] = 'toolbar2' # turn toolbars back on for subsequent figures
        #return ias # this prints a whole bunch to screen if not bound to a var, kinda annoying, not too useful anyway
        # use wx's gcf().canvas.Refresh() to update the window, if it doesn't do so automatically when you modify its contents. Or, you can use matplotlib's pl.draw() command instead, better yet, use fig.canvas.draw() explicitly
    '''

class STA(RevCorr):
    """Spike-triggered average revcorr object"""
    def calc(self):
        RevCorr.calc(self) # run the base calc() steps first
        #sys.stdout.write('n%d' % self.neuron.id) # prevents trailing space and newline
        # init a 3D matrix to store the STA at each timepoint. rf == 'receptive field'
        self.rf = np.zeros([self.nt, self.height, self.width], dtype=np.float64)
        # converting from uint8 to float64 seems to speed up mean() method a bit:
        #frames = np.float64(self.movie.frames)
        frames = self.movie.frames
        tstart = time.clock()
        for ti in self.tis:
            # this can unintentionally introduce -ve valued indices at the left boundary:
            rcdini = self.rcdini - ti*self.ndinperframe
            rcdini = rcdini[rcdini >= 0] # remove any -ve valued indices
            # get the din values (frame indices) at the rcdini for this timepoint:
            frameis = self.experiment.din[rcdini, 1]
            """
            In ptc15, we erroneously duplicated the first frame of the mseq movies at the
            end, giving us one more frame (0 to 65535 for mseq32) than we should have had (0
            to 65534 for mseq32). We're now using the correct movies, but the din for Cat 15
            mseq experiments still have those erroneous frame indices (65535 and 16383 for
            mseq32 and mseq16 respectively), so we'll just ignore them for revcorr purposes.
            """
            if 'mseq32' in self.movie.static.fname.lower():
                frameis = frameis[frameis != 65535] # remove all occurences of 65535
            elif 'mseq16' in self.movie.static.fname.lower():
                frameis = frameis[frameis != 16383] # remove all occurences of 16383
            # take the mean of all the frames at this timepoint:
            # slowest way:
            #self.rf[ti] = frames[frameis].mean(axis=0)
            # faster way:
            #self.rf[ti] = frames.mean(axis=0)
            # much faster way:
            # collect the relevant frames for this timepoint, take is much faster than
            # direct indexing
            pickedframes = frames.take(frameis, axis=0)
            self.rf[ti] = mean_accum(pickedframes)
            #self.rf[ti] = mean_accum2(frames, frameis)
        self.done = True

    def plot(self, interp='nearest', normed=True, scale=2.0):
        win = RevCorr.plot(self, interp=interp, normed=normed, title=lastcmd(),
                           scale=scale)
        return win # necessary in IPython
    plot.__doc__ = RevCorr.plot.__doc__


class STC(RevCorr):
    """Spike-triggered correlation revcorr object"""
    def calc(self):
        print 'INCOMPLETE'
        super(STC, self).calc() # run the general calc() steps

    def plot(self):
        print 'INCOMPLETE'
        vals = super(STC, self).plot()
        return vals
    plot.__doc__ = RevCorr.plot.__doc__


class NeuronRevCorr(object):
    """Mix-in class that defines the reverse correlation related Neuron methods"""
    def sta(self, experiment=None, **kwargs):
        """Returns an existing STA RevCorr object, or creates a new one if necessary"""
        try:
            self._stas
        except AttributeError: # self._stas doesn't exist yet
            self._stas = [] # create a list that'll hold STA objects
        if experiment == None: # no Experiment was passed, use the first experiment this Neuron was involved in
            experiment = self.sort.r.e[0]
        else:
            if type(experiment) == int:
                # assume experiment is an Experiment id, get the associated object
                experiment = self.sort.r.e[experiment]
            # else: experiment is probably an Experiment object
        stao = STA(neuron=self, experiment=experiment, **kwargs) # init a new STA object
        for sta in self._stas:
            if stao == sta: # need to define special == method for RevCorr class
                if sta.done:
                    return sta # returns the first STA object whose attributes match what's desired, and whose calculation done flag is set. This saves on calc() time and avoids wasting memory with unnecessary sta objects
                else:
                    sta.calc() # re-run its calc()
                    return sta
        stao.calc() # no matching STA was found, calculate it
        if stao.done: # if calc() was allowed to go to completion
            self._stas.append(stao) # add it to the STA object list
        return stao # return it, even if it isn't done
    sta.__doc__ += '\n\n**kwargs:\n'
    sta.__doc__ += getargstr(STA.__init__)

    def stc(self, experiment=None, **kwargs):
        """Returns an existing STC RevCorr object, or creates a new one if necessary"""
        try:
            self._stcs
        except AttributeError: # self._stcs doesn't exist yet
            self._stcs = [] # create a list that'll hold STC objects
        if experiment == None: # no Experiment was passed, use the first experiment this Neuron was involved in
            experiment = self.sort.r.e[0]
        else:
            try: # assume experiment is an Experiment id, get the associated object
                experiment = self.sort.r.e[experiment]
            except KeyError: # experiment is probably an Experiment object
                pass
        stco = STC(neuron=self, experiment=experiment, **kwargs) # init a new STC object
        for stc in self._stcs:
            if stco == stc: # need to define special == method for class STC
                return stc # returns the first STC object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._stcs
        stco.calc() # no matching STC was found, calculate it
        self._stcs.append(stco) # add it to the STC object list
        return stco
    stc.__doc__ += '\n\n**kwargs:\n'
    stc.__doc__ += getargstr(STC.__init__)


class Tune(object):
    """Stimulus tuning analysis object"""
    def __init__(self, neuron=None, experiment=None):
        self.neuron = neuron
        self.experiment = experiment
        self.done = False # hasn't yet successfully completed its calc() method
        
    def calc(self, tdelay=None):
        """tdelay: time delay in us to use between stimulus and response"""
        spikes = self.neuron.spikes
        din = self.experiment.din
        ndin = len(din)
        sweepis = np.unique(din[:, 1]) # all possible sweep indices
        if tdelay == None:
            if 'flash' in self.experiment.name: # flashgrating or flashbar
                tdelay = 40000 # akin to a revcorr timepoint for STA
            else:
                tdelay = 0
        self.tdelay = tdelay
        # find positions of each sweep index in the din, and generate array of tranges
        # during which that stimulus condition was on
        self.tranges = {} # index into using sweepi
        self.counts = {} # index into using sweepi
        for sweepi in sweepis:
            dinis = np.where(din[:, 1] == sweepi)[0] # screen refresh indices
            deltaiis = np.where(np.diff(dinis) != 1)[0] # look for non-consecutive values
            startiis = np.insert(deltaiis+1, 0, 0) # prepend with 0
            endiis = np.append(deltaiis, len(dinis)-1)
            rangeiis = np.vstack([startiis, endiis]).T
            rangeis = dinis[rangeiis]
            rangeis[:, 1] += 1 # end inclusive
            rangeis[-1, 1] = min(rangeis[-1, 1], ndin-1) # except for very end
            tranges = din[rangeis, 0] + tdelay
            self.tranges[sweepi] = tranges
            # find which spike indices the start and end of each trange would fall
            # between. Take difference between those two spike indices to get spike
            # count for that trange. Repeat for all tranges for this sweepi to get
            # array of spike counts, one for each trange.
            self.counts[sweepi] = np.diff(spikes.searchsorted(tranges), axis=1).flatten()
        self.done = True
        
    def plot(self, var='ori', fixed=None):
        """var: string name of variable you want to plot a tuning curve for
        fixed: dict with keys containing names of vars to keep fixed when building tuning
        curve, and values containing each var's value(s) to fix at
        
        Ex: r71.n[1].tune().plot('phase0', fixed={'ori':138, 'sfreqCycDeg':[0.4, 0.8]})
        """
        if not self.done:
            self.calc(tdelay=self.tdelay)
        if fixed != None:
            fixedsweepis = []
            for fixedvar, fixedvals in fixed.items():
                vals = self.experiment.sweeptable.data[fixedvar]
                if fixedvar == 'ori': # correct for orientation offset by adding
                    if (vals > 180).any():
                        maxori = 360
                    else:
                        maxori = 180
                    vals = vals.copy() # don't modify the sweeptable!
                    vals += self.experiment.s.orioff # static parameter
                    vals %= maxori
                sweepis = []
                for fixedval in toiter(fixedvals):
                    sweepis.append(np.where(vals == fixedval)[0])
                sweepis = np.concatenate(sweepis)
                fixedsweepis.append(sweepis)
            # intersect all fixedvar arrays in fixedsweepis:
            fixedsweepis = core.intersect1d(fixedsweepis)
            #print(fixedsweepis)
        # get values for var at all unique sweep indices:
        vals = self.experiment.sweeptable.data[var]
        if var == 'ori': # correct for orientation offset by adding
            if (vals > 180).any():
                maxori = 360
            else:
                maxori = 180
            vals = vals.copy() # don't modify the sweeptable!
            vals += self.experiment.s.orioff # static parameter
            vals %= maxori
        x = np.unique(vals) # x axis
        y = np.zeros(len(x), dtype=int) # spike counts for each variable value
        for vali, val in enumerate(x):
            sweepis = np.where(vals == val)[0]
            if fixed != None:
                sweepis = np.intersect1d(sweepis, fixedsweepis, assume_unique=True)
                print sweepis
            for sweepi in sweepis:
                y[vali] += self.counts[sweepi].sum()
        # create a new figure:
        f = pl.figure()
        a = f.add_subplot(111)
        a.plot(x, y, 'k.-')
        a.set_xlabel(var)
        a.set_ylabel('spike count')
        titlestr = lastcmd()
        titlestr += ' nid%d' % self.neuron.id
        a.set_title(titlestr)
        f.canvas.window().setWindowTitle(titlestr)
        a.text(0.99, 0.99, 'peak=(%s, %s)' % (x[y.argmax()], y.max()),
               transform=a.transAxes,
               horizontalalignment='right',
               verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        self.f = f
        self.x, self.y = x, y
        return self


class NeuronTune(object):
    """Mix-in class that defines stimulus tuning analysis method"""
    def tune(self, eid=0, tdelay=None):
        """Return stimulus tuning analysis object"""
        experiment = self.sort.r.e[eid] # get experiment from parent recording
        tuneo = Tune(self, experiment)
        tuneo.calc(tdelay)
        return tuneo


class Neuron(NeuronTune,
             NeuronRevCorr,
             NeuronRate,
             NeuronCode,
             NeuronXCorr,
             NeuronBasics,
             BaseNeuron):
    """Inherit all the Neuron classes into a single Neuron class"""
    pass


class TrackNeuron(NeuronRate,
                  NeuronCode,
                  NeuronXCorr,
                  NeuronBasics):
    """A neuron that spans all recordings in a track, and therefore can't have any
    experiment-specific analyses"""
    def __init__(self, sort):
        self.level = 4 # level in the hierarchy, just below TrackSort
        self.path = sort.path
        self.sort = sort # a TrackSort
