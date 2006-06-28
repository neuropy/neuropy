"""Defines the Neuron class"""

print 'importing Neuron'

from Core import *

class Neuron(object):
    """A Neuron object\'s spike data spans all the Experiments within a Recording.
    If different Recordings have Rips with the same name, you can assume that the
    same spike template was used for all of those Recordings, and that therefore
    the neuron ids are the same"""
    from Rip import Rip
    def __init__(self, id=None, name=None, parent=Rip): # neuron names don't include the '.spk' ending, although neuron filenames do
        self.level = 5 # level in the hierarchy
        #self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.rip = parent() # init parent Rip object
        except TypeError: # parent is an instance, not a class
            self.rip = parent # save parent Rip object
        if id is not None:
            name = self.id2name(self.rip.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'neuron id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.rip.path
    '''
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.rip.writetree(string)
    '''
    def id2name(self, path, id):
        #if len(str(id)) == 1: # if id is only 1 digit long
        #    id = '0'+str(id) # add a leading zero
        name = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(path) if os.path.isfile(path+fname) and \
               ( fname.find('_t'+str(id)+'.spk')!=-1 or fname.find('_t0'+str(id)+'.spk')!=-1 or fname.find('_t00'+str(id)+'.spk')!=-1 ) ] # have to deal with leading zero ids, go up to 3 digit ids, should really use a re to do this properly...
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Neuron id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        try:
            id = name[name.rindex('_t')+2::] # everything from just after the last '_t' to the end of the neuron name, index() raises ValueError if it can't be found
        except ValueError:
            raise ValueError, 'Badly formatted Neuron name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self): # or loadspikes()?
        f = file(self.path + self.name + '.spk', 'rb') # open the spike file for reading in binary mode
        self.spikes = np.fromfile(f, dtype=np.int64) # read in all spike times in us
        f.close()
        self.nspikes = len(self.spikes)
        #self.results = {} # a dictionary to store results in
        #treestr = self.level*TAB + self.name + '/'
        #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
    def cut(self, *args):
        """Returns all of the Neuron's spike times where tstart <= spikes <= tend

        args can be: nothing (returns all spikes), tstart, or (tstart, tend)
        None and 0 are shorthand for 'from first spike'
        None and -1 are shorthand for 'to last spike'"""
        tstart = None
        tend = None
        if len(args) == 0: # passed nothing
            pass
        elif len(args) == 1: # passed None, or just tstart, or a (tstart, tend) tuple/list
            if not iterable(args[0]):
                if args[0] == None:
                    pass
                else: # just tstart was passed
                    tstart = args[0]
            elif len(args[0]) == 2: # it's a tuple/list
                tstart = args[0][0]
                tend = args[0][1]
            else:
                raise ValueError, 'tuple/list is too long'
        elif len(args) == 2: # passed tstart and tend as separate args
            tstart = args[0]
            tend = args[1]
        else:
            raise ValueError, 'too many arguments'
        if tstart in [None, 0]: # shorthand for "from first spike" - would be problematic if a spike existed at t=0
            tstart = self.spikes[0]
        if tend in [None, -1]: # shorthand for "to last spike" - would be problematic if a spike existed at t=-1
            tend = self.spikes[-1]
        '''
        # this is what we're trying to do:
        return self.spikes[ (self.spikes >= tstart) & (self.spikes <= tend) ]
        self.searchsorted(values) method does it faster. It returns an index where the value would fit in self. The index is such that self[index-1] < value <= self[index]. In this formula self[self.size]=inf and self[-1]= -inf
        Another possibility could be to use a masked array instead?
        '''
        lo, hi = self.spikes.searchsorted([tstart, tend]) # returns indices where tstart and tend would fit in spikes
        if tend  == self.spikes[min(hi, self.nspikes-1)]: # if tend matches a spike (protect from going out of index bounds when checking)
            hi += 1 # inc to include a spike if it happens to exactly equal tend. This gives us end inclusion
            hi = min(hi, self.nspikes) # limit hi to max slice index (==max value index + 1)
        return self.spikes[ lo : hi ] # slice it

    def cutrel(self, *args):
        """Cuts Neuron spike times and returns them relative to tstart"""
        if len(args) == 0: # passed nothing
            tstart = self.spikes[0]
            tend = self.spikes[-1]
        elif len(args) == 1: # passed None or a tuple
            if args[0] == None:
                tstart = self.spikes[0]
                tend = self.spikes[-1]
            else: # it's a tuple
                tstart = args[0][0]
                tend = args[0][1]
        elif len(args) == 2: # passed tstart and tend as separate args
            tstart = args[0]
            tend = args[1]
        else:
            raise RuntimeError, 'Too many arguments'
        #if tstart == 0: # shorthand for "from first spike"
        #   tstart = self.spikes[0]
        if tend == -1: # shorthand for "to last spike"
            tend = self.spikes[-1]
        tstart = np.int64(round(tstart)) # let's keep all the returned spikes as integers, us is more than accurate enough
        return self.cut(tstart, tend) - tstart

    def isi(self, trange=None):
        """Returns the inter-spike interval of the Neuron's spike train in trange"""
        #try:
        #    return self._isi # see if it's already been calculated
        #except AttributeError:
        #    self._isi = np.diff(self.spikes) # store it
        #    return self._isi
        return np.diff(self.cut(trange))

    def iisii(self, trange=None):
        """Returns the inter-ISI interval (the differences between consecutive ISIs)
        of the Neuron's spike train in trange"""
        #try:
        #    return self._iisii # see if it's already been calculated
        #except AttributeError:
        #    self._iisii = np.diff(self.isi()) # store it
        #    return self._iisii
        return np.diff(self.isi(trange))

    class Code(object):
        """Abstract spike code class"""
        def __init__(self, neuron=None, trange=None):
            self.neuron = neuron
            if trange == None:
                trange = self.neuron.spikes[0], self.neuron.spikes[-1]
            self.trange = trange
        def __eq__(self, other):
            selfd = self.__dict__.copy()
            otherd = other.__dict__.copy()
            # Delete their c and t attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
            [ d.__delitem__(key) for d in [selfd, otherd] for key in ['c', 't'] if d.has_key(key) ]
            if type(self) == type(other) and selfd == otherd:
                return True
            else:
                return False
        def plot(self):
            # plot some kind of long grid of white and black elements?
            pl.figure()

    class BinaryCode(Code):
        """Quantizes the spike train into a binary signal with a given time resolution"""
        def __init__(self, neuron=None, trange=None, tres=20000):
            super(Neuron.BinaryCode, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'binary'
            self.tres = tres
        def calc(self):
            # make the start of the timepoints be an even multiple of self.tres. Round down to the nearest multiple. Do the same for the end of the timepoints, but round up. This way, timepoints will line up for different code objects
            tstart = self.trange[0] - (self.trange[0] % self.tres)
            tend   = self.trange[1] - (self.trange[1] % self.tres)
            self.t = np.arange( tstart, tend+self.tres, self.tres ) # t sequence demarcates left bin edges, add tres to trange[1] to make t end inclusive
            s = self.neuron.cut(self.trange) # spike times
            self.c = np.zeros(len(self.t)) # binary code signal
            self.c[np.unique(self.t.searchsorted(s)) - 1] = 1 # dec index by 1 so that you get indices that point to the most recent bin edge. For each bin that has 1 or more spikes in it, set its bit to 1
        def plot(self):
            super(Neuron.BinaryCode, self).plot()
            pl.title('neuron %d - binary spike code' % self.neuron.id)

    def code(self, kind='binary', **kwargs):
        """Returns an existing Code object, or creates a new one if necessary"""
        try:
            self._codes
        except AttributeError: # self._codes doesn't exist yet
            self._codes = [] # create a list that'll hold Code objects
        if kind == 'binary':
            co = Neuron.BinaryCode(neuron=self, **kwargs) # init a new BinaryCode object
        else:
            raise ValueError, 'Unknown kind: %s' % repr(self.kind)
        for code in self._codes:
            if co == code: # need to define special == method for class Code()
                return code # returns the first Code object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codes
        co.calc() # no matching Rate was found, calculate it
        self._codes.append(co) # add it to the Code object list
        return co
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nbinary: '+getargstr(BinaryCode.__init__)

    class Rate(object):
        """Abstract firing rate class"""
        def __init__(self, neuron=None, trange=None):
            self.neuron = neuron
            if trange == None:
                trange = self.neuron.spikes[0], self.neuron.spikes[-1]
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

    class BinRate(Rate):
        """Uses simple binning to calculate firing rate"""
        def __init__(self, neuron=None, trange=None, tres=100000):
            super(Neuron.BinRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'bin'
            self.tres = tres
        def calc(self):
            # make the start of the timepoints be an even multiple of self.tres. Round down to the nearest multiple. Do the same for the end of the timepoints. This way, timepoints will line up for different code objects
            tstart = self.trange[0] - (self.trange[0] % self.tres)
            tend   = self.trange[1] - (self.trange[1] % self.tres)
            t = np.arange( tstart, tend+self.tres, self.tres ) # t sequence demarcates left bin edges, add tres to trange[1] to make t end inclusive
            s = self.neuron.cut(trange) # spike times
            self.r, self.t = histogramSorted(self.neuron.spikes, bins=t) # assumes spikes are in chrono order
            self.r = self.r / float(self.tres) * 1000000 # spikes/sec
        def plot(self):
            super(Neuron.BinRate, self).plot()
            pl.title('neuron %d - binned spike rate' % self.neuron.id)

    class nISIRate(Rate):
        """Uses nisi inter spike intervals to calculate firing rate, with a fixed number of spikes nisi+1 per bin. nisi == 1 is the ISI rate"""
        def __init__(self, neuron=None, trange=None, nisi=3, tres=50000, interp='sah'):
            super(Neuron.nISIRate, self).__init__(neuron=neuron, trange=trange)
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
                raise ValueError, 'unknown interpolation method: %s' % self.interp

        def plot(self):
            super(Neuron.nISIRate, self).plot()
            pl.title('neuron %d - %d-inter-spike-interval spike rate, %s interpolation' % (self.neuron.id, self.nisi, self.interp))

    class WnISIRate(Rate):
        """Uses a weighted sum of various n inter spike intervals to calculate rate"""
        pass
        #def plot(self)

    class IDPRate(Rate):
        """Instantaneous discharge probability. See Pauluis and Baker, 2000"""
        def __init__(self, neuron=None, IDP_a=4, tres=50000):
            super(Neuron.IDPRate, self).__init__(neuron=neuron, trange=trange)
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

            # Step 2: extend back half an ISI from sudden rate increase, and forward half an ISI from sudden rate decrease

            # Step 3: sample our new ISI rate with extended intervals with sufficient resolution
            isi = self.neuron.rate(kind='nisi', nisi=1, interp=None) # ISI rate object with with no interpolation
            isi.r
            isi.t
            isi = sah(isi)

            # use np.piecewise() ?
        pass
        #def plot(self)

    class GaussRate(Rate):
        """Uses a sliding Gaussian window to calculate firing rate"""
        def __init__(self, neuron=None, trange=None, width=200000):
            super(Neuron.GaussRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'gauss'
            self.width = width
        def calc(self):
            pass
        def plot(self):
            super(Neuron.GaussRate, self).plot()
            pl.title('Gaussian sliding window spike rate')

    class RectRate(Rate):
        """Uses a sliding rectangular window to calculate firing rate"""
        def __init__(self, neuron=None, trange=None, width=200000):
            super(Neuron.RectRate, self).__init__(neuron=neuron, trange=trange)
            self.kind = 'rect'
            self.width = width
        def calc(self):
            pass
        def plot(self):
            super(Neuron.RectRate, self).plot()
            pl.title('rectangular sliding window spike rate')

    def rate(self, kind='nisi', **kwargs):
        """Returns an existing Rate object, or creates a new one if necessary"""
        try:
            self._rates
        except AttributeError: # self._rates doesn't exist yet
            self._rates = [] # create a list that'll hold Rate objects
        if kind == 'bin':
            ro = Neuron.BinRate(neuron=self, **kwargs) # init a new BinRate object
        elif kind == 'nisi':
            ro = Neuron.nISIRate(neuron=self, **kwargs) # init a new nISIRate object
        elif kind == 'wnisi':
            ro = Neuron.WnISIRate(neuron=self, **kwargs) # init a new WnISIRate object
        elif kind == 'idp':
            ro = Neuron.IDPRate(neuron=self, **kwargs) # init a new IDPRate object
        elif kind == 'gauss':
            ro = Neuron.GaussRate(neuron=self, **kwargs) # init a new GaussRate object
        elif kind == 'rect':
            ro = Neuron.RectRate(neuron=self, **kwargs) # init a new RectRate object
        else:
            raise ValueError, 'Unknown kind: %s' % repr(self.kind)
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

    class RatePDF(object):
        """Firing rate probability distribution function object"""
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
                for (i,r) in enumerate(safe):
                    if r == 0:
                        safe[i] = 0.1 # set to 0.1 Hz
                self.logrrange = np.log10(tuple(safe)) # convert back to tuple, is now safe to take log
                # r sequence demarcates left rate bin edges
                r = np.logspace(start=self.logrrange[0], stop=self.logrrange[1], num=self.nbins, endpoint=True, base=10.0)
            elif self.scale == 'linear':
                r = np.linspace(start=self.rrange[0], stop=self.rrange[1], num=self.nbins, endpoint=True)
            else:
                raise ValueError, 'Unknown scale: %s' % repr(scale)
            self.n, self.r = np.histogram(self.rate.r, bins=r, normed=self.normed)
        def plot(self):
            pl.figure()
            if self.scale == 'log':
                barwidth = list(np.diff(self.r)) # each bar will have a different width, convert to list so you can append
                # need to add one more entry to barwidth to the end to get nbins of them:
                #barwidth.append(barwidth[-1]) # not exactly correct
                logbinwidth = (self.logrrange[1]-self.logrrange[0]) / float(self.nbins)
                barwidth.append(10**(self.logrrange[1]+logbinwidth)-self.r[-1]) # should be exactly correct
            elif self.scale == 'linear':
                barwidth = (self.rrange[1]-self.rrange[0]) / float(self.nbins)
            else:
                raise ValueError, 'Unknown scale: %s' % repr(scale)
            #pl.hist(self.n, bins=self.r, normed=0, bottom=0, width=None, hold=False) # doesn't seem to work
            pl.bar(left=self.r, height=self.n, width=barwidth)
            pl.axes().set_xscale(self.scale, basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
            pl.title('neuron %d - %s spike rate PDF' % (self.neuron.id, self.rate.kind))
            if self.normed:
                if self.normed == 'pmf': # it's a probability mass function
                    pl.ylabel('probability mass')
                else: # it's a probability density function
                    pl.ylabel('probability density')
            else:
                pl.ylabel('count')
            pl.xlabel('spike rate')

    def ratepdf(self, **kwargs):
        """Returns an existing RatePDF object, or creates a new one if necessary"""
        try:
            self._ratepdfs
        except AttributeError: # self._ratepdfs doesn't exist yet
            self._ratepdfs = [] # create a list that'll hold RatePDF objects
        rpdf = Neuron.RatePDF(neuron=self, **kwargs)
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

    def raster(self):
        # use pylab.vlines()
        pass

