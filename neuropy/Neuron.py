"""Defines the Neuron class and all of its support classes"""

# TODO: stop doing this!
from Core import *
#from Core import Cat15Movie
from dimstimSkeletal import Movie

class BaseNeuron(object):
    """A Neuron object's spike data spans all the Experiments within a Recording.
    If different Recordings have Sorts with the same name, you can assume that the
    same spike template was used for all of those Recordings, and that therefore
    the neuron ids are the same"""

    from Sort import Sort

    def __init__(self, id=None, name=None, parent=Sort): # neuron names don't include the '.spk' ending, neuron filenames do
        self.level = 5 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.sort = parent() # init parent Sort object
        except TypeError: # parent is an instance, not a class
            self.sort = parent # save parent Sort object
        if id is not None:
            name = self.id2name(self.sort.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'Neuron id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.sort.path
    '''
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.sort.writetree(string)
    '''
    def id2name(self, path, id):
        #if len(str(id)) == 1: # if id is only 1 digit long
        #    id = '0'+str(id) # add a leading zero
        name = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(path)
                 if os.path.isfile(os.path.join(path, fname))
                 and
                 (    '_t%s.spk' % id in fname
                   or '_t0%s.spk' % id in fname
                   or '_t00%s.spk' % id in fname )
                 ] # have to deal with leading zero ids, go up to 3 digit ids, should really use a re to do this properly...
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Neuron id: %s' % id
        else:
            return name[0] # pull the string out of the list

    def name2id(self, name):
        try:
            id = name[name.rindex('_t')+2::] # everything from just after the last '_t' to the end of the neuron name, rindex() raises ValueError if it can't be found
        except ValueError:
            raise ValueError, 'Badly formatted Neuron name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id

    def load(self):
        #treestr = self.level*TAB + self.name + '/'
        #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        f = file(os.path.join(self.path, self.name) + '.spk', 'rb') # open the spike file for reading in binary mode
        self.spikes = np.fromfile(f, dtype=np.int64) # read in all spike times in us
        f.close()
        self.nspikes = len(self.spikes)
        if self.nspikes == 0: # this Neuron has no spikes, not much of a Neuron as far as this Recording is concerened
            self.spikes = np.append(self.nspikes, None) # add None to empty spike list. Why None? Why not just leave it empty?
        #self.results = {} # a dictionary to store results in
        self.trange = self.spikes[0], self.spikes[-1]
        # then, maybe add something that loads the template for this neuron, as well as its modelled (or just guesstimated) location in 3D (or just 2D) coordinates in um, as well as cell type potentially

    def cut(self, *args):
        """Returns a view of the Neuron's spike times where tstart <= spikes <= tend
        *args can be: nothing (returns all spikes), None, tstart, or (tstart, tend)
        None and 0 for tstart are shorthand for 'from first spike'
        None and -1 for tend are shorthand for 'to last spike'"""
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
        tstart, tend = self._parse_cutargs(*args)
        cutspikes = self.cut(tstart, tend)
        tstart = np.int64(round(tstart)) # let's keep all the returned spikes as integers, us is more than accurate enough
        return cutspikes - tstart

    def _parse_cutargs(self, *args):
        """Takes set of args as passed to cut or cutrel and returns (tstart, tend)"""
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
                raise ValueError, 'sequence is too long'
        elif len(args) == 2: # passed tstart and tend as separate args
            tstart = args[0]
            tend = args[1]
        else:
            raise ValueError, 'too many arguments'
        if tstart in [None, 0]: # shorthand for "from first spike" - would be problematic if a spike existed at t=0
            tstart = self.spikes[0]
        if tend in [None, -1]: # shorthand for "to last spike" - would be problematic if a spike existed at t=-1
            tend = self.spikes[-1]
        return (tstart, tend)

    def copy(self):
        """Returns a copy of the Neuron"""
        return copy(self)

    def append(self, others):
        """Maybe should split this up into an .append() and an .extend() method.

        Appends the spike times of self and other Neurons and returns a copy of the uberneuron.

        MAYBE I SHOULD MAKE THIS AN IN-PLACE operation, as self.append() would imply. In other words,
        don't return a copy? But that would modify the original Neuron. Maybe it should be up to you
        to make a copy of the Neuron before appending other Neurons to it. That sounds best.

        All Neurons must have the same id and be of the same Track. The user needs to ensure
        that the same template was used to sort all the Neurons

        If you want to recover the original spike times so as to compare with stimuli,
        you need to cut the spikes with the appropriate trange in tranges, and then subtract
        the appropriate offset in offsets.

        You can't simply append neuron spike times, need to add some kind of huge ass time offset between Recordings?
        Also need to have tranges instead of just a trange, to skip over/get around the huge ass time offset.
        Or maybe just have a 1 us offset tween consecutive Recordings? In the end, I decided to make the spike times from the nth
        Recording to be the sum of the last spike times of all the 0 to n-1 Recordings preceding it. This is like pretending
        that each consecutive Recording started immediately after the preceding one ended.

        HOW THE HELL DO YOU FIND THE INDEX OF A VALUE IN A NUMPY ARRAY? a==5 gives you a long boolean array,
        how do you find the indices where that boolean is True? use .nonzero() or np.where()

        AS OF OCT 11, 2006, DOESN'T LOOK LIKE THIS METHOD IS ACTUALLY BEING USED YET. IT'LL BE USEFUL ONCE I START
        DOING ANALYSES OVER MULTIPLE RECORDINGS TO INCREASE STATS SIGNIFICANCE"""

        others = makeiter(others)
        neurons = [self]
        neurons.extend(others) # put self and all others into a single list of Neurons
        for other in others:
            assert self.id == other.id, 'Neuron ids are different'
            #assert self.sort == other.sort, 'Sorts are different' # forget this, only the user can really know if they used the same template
            assert self.sort.r.t == other.sort.r.t, 'Tracks are different'
        uberneuron = self.copy() # create a copy of self
        uberneuron.spikes = array([], dtype=np.int64) # clear its spikes attrib
        uberneuron.name = '' # clear it
        uberneuron.path = '' # clear it
        uberneuron.sort = [] # clear, and init a list
        # sort the Neurons according to their Recording id
        rids = [ (n.sort.r.id, i) for i, n in enumerate(neurons) ]
        rids.sort() # sorts according to first entry in each tuple
        sortedis = [ i for r, i in rids ]
        neurons = [ neurons[i] for i in sortedis ]
        offset = 0
        uberneuron.offsets = [] # new field to store the offsets that have been added to the spike times for each Recording
        uberneuron.tranges = [] # new field to store the tranges of all the neurons appended to make the uberneuron
        for n in neurons:
            uberneuron.spikes = np.append(uberneuron.spikes, n.spikes+offset)
            # don't forget to add offset to tranges as well!
            try: # neuron is a ConstrainedNeuron with a tranges attrib
                uberneuron.tranges.extend([ (trange[0]+offset, trange[1]+offset) for trange in n.tranges ])
                uberneuron.offsets.extend(offset*len(n.tranges)) # add n.tranges copies of offset to offset list, one for each trange in this ConstrainedNeuron
            except AttributeError: # no tranges attrib, neuron is a normal Neuron
                uberneuron.tranges.append((n.trange[0]+offset, n.trange[1]+offset))
                uberneuron.offsets.append(offset)
            uberneuron.name += n.name + ', ' # keep it as a single string
            uberneuron.path += n.path + ', ' # keep it as a single string
            uberneuron.sort.append(n.sort) # concatenate all sorts into a list
            offset += n.spikes[-1] # for next Recordings spikes, inc offset by the last spike in this Recording
        uberneuron.spikes.sort() # make sure spiketimes remain sorted
        uberneuron.nspikes = len(uberneuron.spikes) # update it
        uberneuron.trange = uberneuron.spikes[0], uberneuron.spikes[-1] # update it
        uberneuron.name = uberneuron.name[0:-2] # get rid of trailing ', '
        uberneuron.path = uberneuron.path[0:-2]
        return uberneuron

    def isi(self, trange=None):
        """Returns the inter-spike intervals of the Neuron's spike train in trange"""
        return diff(self.cut(trange))

    def iisii(self, trange=None):
        """Returns the inter-ISI intervals (the differences between consecutive ISIs)
        of the Neuron's spike train in trange.
        See delta def'n in: 2002 Segev, et al - Long term behavior of lithographically..."""
        return diff(self.isi(trange))


class XCorr(object):
    def __init__(self, n1=None, n2=None, trange=(-100000, 100000)):
        """Cross-correlation object. n1 has to be a Neuron, n2 can be a Neuron or a Neuron id"""
        self.n1 = n1
        try: # assume n2 is a Neuron id from the same Sort as n1, get the associated Neuron object
            n2 = self.n1.sort.n[n2]
        except KeyError: # n2 is probably a Neuron object
            pass
        self.n2 = n2
        self.trange = trange

    def calc(self):
        dts = []
        for spike in self.n1.spikes:
            trangei = self.n2.spikes.searchsorted(spike+self.trange) # find where the trange around this spike would fit in other.spikes
            dt = self.n2.spikes[trangei[0]:trangei[1]] - spike # find dt between this spike and only those other.spikes that are in trange of this spike
            dts.extend(dt)
        self.dts = array(dts)
        return self.dts

        # could use some weave code to speed this up
        '''
        code = """
               #line 193 "Neuron.py" (This is only useful for debugging)
               double tmp, err, diff;
               err = 0.0;
               for (int spikei=0; spikei<=nspikes; i++) {
                   trangei = otherspikes.searchsorted(spike+trange(0), spike+trange(1))
                   out = trangei
               }
               return_val = out;
               """
        # compiler keyword only needed on windows with MSVC installed
        out = weave.inline(code,
                           [selfspikes, otherspikes, nspikes, trange],
                           type_converters=weave.converters.blitz,
                           compiler = 'msvc?')
        '''
    def plot(self, nbins=100, figsize=(6.5, 6.5), style='count'):
        """style can be 'count' or 'rate'"""
        f = figure(figsize=figsize)
        a = f.add_subplot(111)
        n, t = histogram(self.dts, bins=nbins)
        self.n = n
        self.t = t
        barwidth = (t.max()-t.min())/float(nbins)
        if style == 'rate':
            n = n / float(barwidth)
        bar(left=t, height=n, width=barwidth)
        gcfm().frame.SetTitle(lastcmd())
        #gcfm().frame.SetTitle('r%d.n[%d].xcorr(%d)' % (self.n1.sort.r.id, self.n1.id, self.n2.id))
        a.set_title('n%d spikes relative to n%d spikes' % (self.n2.id, self.n1.id))
        a.set_xlabel('time (msec)')
        if style == 'rate':
            a.set_ylabel('spike rate (Hz)')
        else:
            a.set_ylabel('bin count')
        xticks = a.get_xticks() / 1e3 # convert from us to ms
        xticklabels = []
        [ xticklabels.append('%d' % xtick) for xtick in xticks ] # truncate floats into ints
        a.set_xticklabels(xticklabels)


class NeuronXCorr(BaseNeuron):
    """Mix-in class that defines the xcorr Neuron method"""
    def xcorr(self, other=None, **kwargs):
        """Returns a cross-correlation object"""
        xco = XCorr(n1=self, n2=other, **kwargs)
        xco.calc()
        return xco
    xcorr.__doc__ += '\n\n**kwargs:'
    xcorr.__doc__ += '\n'+getargstr(XCorr.__init__)


class BaseCode(object):
    """Abstract spike code class. tranges is a list of time ranges to use
    when generating the code"""
    def __init__(self, neuron=None, tranges=None):
        self.neuron = neuron
        if tranges == None:
            self.tranges = [ self.neuron.trange ]
        else:
            self.tranges = tolist(tranges)

    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their c, t, and s attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd] for key in ['c', 't', 's'] if d.has_key(key) ]
        if type(self) == type(other) and selfd == otherd:
            return True
        else:
            return False

    def plot(self):
        # plot some kind of long grid of white and black elements?
        pass


class BinaryCode(BaseCode):
    """Quantizes the spike train, cut according to tranges, into a binary signal with a given time resolution.
    phase in us specifies where to start the codetrain in time, relative to the nearest multiple of tres before each trange.
    Phase is in degrees of a single bin period. -ve phase is leading (codetrain starts earlier in time),
    +ve is lagging (codetrain starts later in time)"""
    def __init__(self, neuron=None, tranges=None, tres=CODETRES, phase=CODEPHASE):
        super(BinaryCode, self).__init__(neuron=neuron, tranges=tranges)
        self.kind = 'binary'
        self.tres = tres
        self.phase = phase # in degrees of tres

    def calc(self):
        self.t = array([], dtype=np.int64) # set up empty arrays with correct dtypes (otherwise, when appending to them later, they'd default to float64s)
        self.s = array([], dtype=np.int64)
        self.c = array([], dtype=np.int8)
        for trange in self.tranges:
            # make the start of the timepoints be an even multiple of self.tres. Round down to the nearest multiple. This way, timepoints will line up for different code objects. Finally, add phase offset relative to this
            tstart = trange[0] - (trange[0] % self.tres) + self.phase/360.0*self.tres # left edge of first code bin
            t = np.int64(np.round(arange(tstart, trange[1]+self.tres, self.tres))) # t sequence demarcates left bin edges, add extra tres to end to make t end inclusive, keep 'em in us integers
            s = self.neuron.cut(trange) # spike times, cut over originally specified trange, not from start to end of newly generated code bin timepoints
            c = zeros(len(t), dtype=np.int8) # init binary code signal
            # searchsorted returns indices where s fits into t. Sometimes more than one spike will fit into the same time bin, which means searchsorted will return multiple occurences of the same index. You can set c at these indices to 1 a multiple number of times, or prolly more efficient, do an np.unique on it to only set each index to 1 once.
            c[np.unique(t.searchsorted(s)) - 1] = 1 # dec index by 1 so that you get indices that point to the most recent bin edge. For each bin that has at least 1 spike in it, set its value to 1
            self.t = np.append(self.t, t) # save 'em
            self.s = np.append(self.s, s)
            self.c = np.append(self.c, c)

    def plot(self):
        super(BinaryCode, self).plot()
        title('neuron %d - binary spike code' % self.neuron.id)


class NeuronCode(BaseNeuron):
    """Mix-in class that defines the spike code related Neuron methods"""
    def code(self, tranges=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns an existing Code object, or creates a new one if necessary"""
        try:
            self._codes
        except AttributeError: # self._codes doesn't exist yet
            self._codes = [] # create a list that'll hold Code objects for this Neuron
        if kind == 'binary':
            co = BinaryCode(neuron=self, tranges=tranges, tres=tres, phase=phase) # init a new BinaryCode object
        else:
            raise ValueError, 'Unknown kind: %r' % self.kind
        for code in self._codes:
            if co == code: # need to define special == method for class Code()
                return code # returns the first Code object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codes
        co.calc() # no matching Code was found, calculate it
        self._codes.append(co) # add it to the Code object list
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
        figure()
        plot(self.t, self.r)
        # diagnostic for comparing interpolated to raw nisi rate:
        #plot(self.rawt, self.rawr, 'r+')
        xlabel('t')
        ylabel('spike rate')


class BinRate(BaseRate):
    """Uses simple binning to calculate firing rate"""
    def __init__(self, neuron=None, trange=None, tres=100000):
        super(BinRate, self).__init__(neuron=neuron, trange=trange)
        self.kind = 'bin'
        self.tres = tres

    def calc(self):
        # make the start of the timepoints be an even multiple of self.tres. Round down to the nearest multiple. Do the same for the end of the timepoints. This way, timepoints will line up for different code objects
        tstart = self.trange[0] - (self.trange[0] % self.tres)
        tend   = self.trange[1] - (self.trange[1] % self.tres)
        t = arange( tstart, tend+self.tres, self.tres ) # t sequence demarcates left bin edges, add tres to trange[1] to make t end inclusive
        s = self.neuron.cut(self.trange) # spike times
        self.r, self.t = histogramSorted(self.neuron.spikes, bins=t) # assumes spikes are in chrono order
        self.r = self.r / float(self.tres) * 1000000 # spikes/sec

    def plot(self):
        super(BinRate, self).plot()
        title('neuron %d - binned spike rate' % self.neuron.id)


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
        self.t = arange(tstart, t[-1], self.tres) # new set of timepoints to interpolate over
        if self.interp == 'sah':
            self.r = sah(t, r, self.t, keep=False)
        elif self.interp == 'linear':
            f = sig.interpolate.interp1d(t, r, kind='linear') # returns an interpolation f'n
            self.r = f(self.t) # interpolate over the new timepoints
            # or maybe try sig.resample() instead
        else:
            raise ValueError, 'unknown interpolation method: %s' % self.interp

    def plot(self):
        super(nISIRate, self).plot()
        title('neuron %d - %d-inter-spike-interval spike rate, %s interpolation' % (self.neuron.id, self.nisi, self.interp))


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
        title('Gaussian sliding window spike rate')


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
        title('rectangular sliding window spike rate')


class RatePDF(object):
    """Firing rate probability distribution function"""
    def __init__(self, neuron=None, rrange=(0, 200), nbins=100, scale='log', normed='pmf', **kwargs):
        # rrange == rate range, ie limits of pdf x axis; nbins == number of rate bins
        self.neuron = neuron
        self.rrange = rrange
        self.nbins = nbins
        self.scale = scale
        self.normed = normed
        self.rate = neuron.rate(**kwargs) # returns a Rate object of some kind

    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their n and r and logrrange attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd] for key in ['n', 'r', 'logrrange'] if d.has_key(key) ]
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
            raise ValueError, 'Unknown scale: %r' % scale
        self.n, self.r = histogram(self.rate.r, bins=r, normed=self.normed)

    def plot(self):
        figure()
        if self.scale == 'log':
            barwidth = list(diff(self.r)) # each bar will have a different width, convert to list so you can append
            # need to add one more entry to barwidth to the end to get nbins of them:
            #barwidth.append(barwidth[-1]) # not exactly correct
            logbinwidth = (self.logrrange[1]-self.logrrange[0]) / float(self.nbins)
            barwidth.append(10**(self.logrrange[1]+logbinwidth) - self.r[-1]) # should be exactly correct
        elif self.scale == 'linear':
            barwidth = (self.rrange[1]-self.rrange[0]) / float(self.nbins)
        else:
            raise ValueError, 'Unknown scale: %r' % scale
        #hist(self.n, bins=self.r, normed=0, bottom=0, width=None, hold=False) # doesn't seem to work
        bar(left=self.r, height=self.n, width=barwidth)
        axes().set_xscale(self.scale, basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
        title('neuron %d - %s spike rate PDF' % (self.neuron.id, self.rate.kind))
        if self.normed:
            if self.normed == 'pmf': # it's a probability mass function
                ylabel('probability mass')
            else: # it's a probability density function
                ylabel('probability density')
        else:
            ylabel('count')
        xlabel('spike rate')


class NeuronRate(BaseNeuron):
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
            raise ValueError, 'Unknown kind: %r' % self.kind
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
            self.movie.load(asarray=True, flip=False) # Load as 3D array instead of as a list of 2D arrays, more convenient for analysis, although will cause memory problems for really big (>1GB movies). Don't flip the movie frames vertically for OpenGL's bottom left origin, since we aren't using OpenGL for analysis
        self.nt = nt # number of revcorr timepoints
        self.tis = range(0, nt, 1) # revcorr timepoint indices
        self.t = [ intround(ti * self.movie.dynamic.sweepSec * 1000) for ti in self.tis ] # revcorr timepoint values, stored in a list, not an array. Bad behaviour happens during __eq__ below if attribs are numpy arrays cuz comparing numpy arrays returns an array of booleans, not just a simple boolean
        self.ndinperframe = intround(self.movie.dynamic.sweepSec * 1000000 / self.experiment.REFRESHTIME)
        #self.movie.frames = asarray(self.movie.frames)
        self.width = self.movie.frames.shape[-1] # (nframes, height, width)
        self.height = self.movie.frames.shape[-2]
        self.done = False # hasn't yet successfully completed its calc() method

    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their rcdini and rf attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd] for key in ['rcdini', 'rf', 'done'] if d.has_key(key) ]
        if type(self) == type(other) and selfd == otherd:
            return True
        else:
            return False

    def calc(self):
        """General calc step that has to be performed for all kinds of reverse correlations"""
        spikes = self.neuron.cut(self.trange)
        self.rcdini = self.experiment.din[:, 0].searchsorted(spikes) - 1 # revcorr dini. Find where the spike times fall in the din, dec so you get indices that point to the most recent din value for each spike
        #self.din = self.experiment.din[rcdini, 1] # get the din (frame indices) at the rcdini

    def plot(self, interp='nearest', normed=True, title='ReceptiveFieldFrame', scale=2.0):
        """Plots the spatiotemporal RF as bitmaps in a wx.Frame"""
        rf = self.rf.copy() # create a copy to manipulate for display purposes, (nt, width, height)
        if normed: # normalize across the timepoints for this RevCorr
            norm = mpl.colors.normalize(vmin=rf.min(), vmax=rf.max(), clip=True) # create a single normalization object to map luminance to the range [0,1]
            rf = norm(rf) # normalize the rf the same way across all timepoints
        else: # don't normalize across timepoints, leave each one to autoscale
            for ti in range(self.nt):
                norm = mpl.colors.normalize(vmin=None, vmax=None, clip=True) # create a normalization object to map luminance to the range [0,1], autoscale
                rf[ti] = norm(rf[ti]) # normalize the rf separately at each timepoint
        cmap = mpl.cm.jet # get a colormap object
        rf = cmap(rf)[::, ::, ::, 0:3] # convert luminance to RGB via the colormap, throw away alpha channel (not used for now in ReceptiveFieldFrame)
        rf = rf * 255 # scale up to 8 bit values
        rf = rf.round().astype(np.uint8) # downcast from float to uint8 for feeding to ReceptiveFieldFrame
        frame = ReceptiveFieldFrame(title=title, rfs=[rf], neurons=[self.neuron], t=self.t, scale=scale)
        frame.Show()
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
        super(STA, self).calc() # run the base calc() steps first
        #sys.stdout.write('n%d' % self.neuron.id) # prevents trailing space and newline
        self.rf = zeros([self.nt, self.height, self.width], dtype=np.float64) # init a 3D matrix to store the STA at each timepoint. rf == 'receptive field'
        #frames = np.float64(self.movie.frames) # converting from uint8 to float64 seems to speed up mean() method a bit
        frames = self.movie.frames
        tstart = time.clock()
        pd = wx.ProgressDialog(title='n%d STA progress' % self.neuron.id, message='', maximum=self.tis[-1], style=1) # create a progress dialog
        for ti in self.tis:
            cont, skip = pd.Update(ti-1, newmsg='timepoint: %dms\nelapsed: %.1fs' % (self.t[ti], time.clock()-tstart))
            if not cont:
                #self.rf = zeros([self.nt, self.height, self.width], dtype=np.float64) # set back to zeros
                pd.Destroy()
                self.done = False
                return
            rcdini = self.rcdini - ti*self.ndinperframe # this can unintentionally introduce -ve valued indices at the left boundary
            rcdini = rcdini[rcdini >= 0] # remove any -ve valued indices. Is this the most efficient way to do this?
            frameis = self.experiment.din[rcdini, 1] # get the din values (frame indices) at the rcdini for this timepoint
            # in Cat 15, we erroneously duplicated the first frame of the mseq movies at the end, giving us one more frame (0 to 65535 for mseq32) than we should have had (0 to 65534 for mseq32). We're now using the correct movies, but the din for Cat 15 mseq experiments still have those erroneous frame indices (65535 and 16383 for mseq32 and mseq16 respectively), so we'll just ignore them for revcorr purposes.
            if MSEQ32 in self.movie.static.fname:
                frameis = frameis[frameis != 65535] # remove all occurences of 65535
            elif MSEQ16 in self.movie.static.fname:
                frameis = frameis[frameis != 16383] # remove all occurences of 16383
            # take the mean of all the frames at this timepoint:
            '''
            # slowest way:
            self.rf[ti] = frames[frameis].mean(axis=0)
            '''
            pickedframes = frames.take(frameis.astype(np.int32), axis=0) # collect the relevant frames for this timepoint, take is much faster than direct indexing, but have to typecast indices to int32, maybe cuz this machine is 32bit?
            '''
            # faster way:
            self.rf[ti] = frames.mean(axis=0)
            '''
            # much faster way:
            self.rf[ti] = mean_accum(pickedframes)
            #self.rf[ti] = mean_accum2(frames, frameis)
        pd.Destroy()
        self.done = True

    def plot(self, interp='nearest', normed=True, scale=2.0):
        super(STA, self).plot(interp=interp, normed=normed,
                              title=lastcmd(),
                              #title='r%d.n[%d].sta().plot(interp=%r, normed=%r, scale=%r)' %
                              #(self.experiment.r.id, self.neuron.id, interp, normed, scale),
                              scale=scale)
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


class NeuronRevCorr(BaseNeuron):
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
            try: # assume experiment is an Experiment id, get the associated object
                experiment = self.sort.r.e[experiment]
            except KeyError: # experiment is probably an Experiment object
                pass
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


class Neuron(NeuronRevCorr,
             NeuronRate,
             NeuronCode,
             NeuronXCorr,
             BaseNeuron):
    """Inherit all the Neuron classes into a single Neuron class"""
    pass


class ConstrainedNeuron(Neuron):
    """A Neuron with its spike times constrained to only those periods of
    time in its (grand)parent Recording during which an Experiment is in progress"""
    def load(self):
        super(ConstrainedNeuron, self).load()
        cspikes = array([], dtype=np.int64) # init a temporary array
        #self.tranges = np.empty((len(self.r.e), 2), dtype=np.int64) # give it shape(numexps, 2)
        self.tranges = []
        if self.sort.r.e != {}: # if it ain't empty, ie if there are Experiments in this Sort's Recording
            for exp in self.sort.r.e.values():
                cspikes = np.append(cspikes, self.cut(exp.trange), axis=0)
                self.tranges.append(exp.trange)
            cspikes.sort() # make sure they're sorted, in case Experiments weren't in temporal order for some reason
            self.tranges.sort() # ditto
            self.spikes = cspikes # save 'em
            self.nspikes = len(self.spikes) # update it
            if self.nspikes == 0: # this Neuron has no spikes, not much of a Neuron as far as this Recording is concerened
                self.spikes = np.append(self.spikes, None) # add None to empty spike list
            self.trange = (self.tranges[0][0], self.tranges[-1][1]) # overwrite inherited self.trange, use start of first trange and end of last trange
