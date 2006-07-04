"""Defines the Experiment class.  This could be split up into multiple base classes
(one for basic init and loading ExperimentBase), one for rates (ExperimentRate), one for codes (ExperimentCodes), etc...), and then combined
using multiple inheritance into a single child Class called Experiment"""

print 'importing Experiment'

from Core import *
from Dimstim.Core import buildSweepTable

class Experiment(object):
    """An Experiment corresponds to a single contiguous VisionEgg stimulus session.
    It contains information about the stimulus during that session, including
    the DIN values, the text header, and any Movies that were involved"""
    from Recording import Recording
    from Neuron import Neuron
    from Neuron import BinaryCode
    from Neuron import RatePDF

    def __init__(self, id=None, name=None, parent=Recording):
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.r = parent() # init parent Recording object
        except TypeError: # parent is an instance, not a class
            self.r = parent # save parent Recording object
        if name is None:
            raise ValueError, 'experiment name can\'t be None'
        self.id = id # not really used by the Experiment class, just there for user's info
        self.name = name
        self.path = self.r.path
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)
    # doesn't need a id2name or name2id method, neither can really be derived from the other in an easy way (although could use re), the id is just chronological (which is also alphabetical) order, at least for now
    def load(self):
        from Recording import Recording # not sure why this has to be imported here again (see above), but seems that it does
        from Movie import Movie, MSEQ32, MSEQ16
        f = file(self.path + self.name + '.din', 'rb') # open the din file for reading in binary mode
        self.din = np.fromfile(f, dtype=np.int64).reshape(-1,2) # reshape to nrows x 2 columns
        f.close()
        try:
            f = file(self.path + self.name + '.textheader', 'r') # open the textheader file for reading
            self.textheader = f.read() # read it all in
            f.close()
        except IOError:
            if type(self.r) is Recording: # parent is a Recording, which normally have textheaders in their Experiments
                warn('Error reading: <%s.textheader>, text header not loaded' % self.name)
            elif type(self.r) is Run: # parent is a Run, which don't have textheaders in their Experiments, so don't print a warning
                pass
            else:
                raise RuntimeError, 'parent is invalid type: %s %s' % (type(self.r), Run)
            self.textheader = '' # set to empty

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

        if self.textheader: # if it isn't empty
            names1 = locals().copy() # namespace before execing the textheader
            exec(self.textheader) # execute the textheader as Python code, maybe add some checks here to prevent changes to filesystem from accidental code, unload the os or sys modules or something?
            names2 = locals().copy() # namespace after
            newnames = [ n2 for n2 in names2 if not names1.has_key(n2) and n2 != 'names1' ] # names that were added to the namespace, excluding the 'names1' name itself
            for newname in newnames:
                self.__setattr__(newname, eval(newname)) # for each variable that was defined in the textheader, bind it as an attribute of this Experiment

            try:
                self.stims = unique(self.playlist) # self.stims is a non-repeating list of object oriented stim objects (Movie is the only possible kind right now) in this Experiment
            except AttributeError: # this was a simple non object-oriented stim, has no playlist
                self.stims = []
            #if len(self.stims) == 1:
            #   self.stims = self.stims[0] # get rid of the list
            '''If you inited a stim object(s) (like a movie) while execing the textheader, you didn't have a chance to pass this exp as the parent in the init. So just set the attribute manually:'''
            for s in self.stims:
                s.e = self
                try: # this'll probably only apply to Movies stim, cuz others won't have fnames
                    if s.name == None:
                        s.name = s.fname # fname should've been defined when loading in the textheader
                except AttributeError: # probably not a Movie stim
                    pass
                # Search self.moviepath string (from textheader) for 'Movies' word (preferably case insensitive). Everything after that is the relative path to your base movies folder. Eg, if self.moviepath = 'C:\\Desktop\\Movies\\reliability\\e\\single\\', then set self.relpath = 'reliability\\e\\single\\'
                try:
                    s.moviepath = s.moviepath.replace('\\','/') # replace annoying double backslashes with single forward slashes, which seem to work
                    s.relpath = s.moviepath[ s.moviepath.index('Movies/')+len('Movies/') :: ]
                    s.path = path + s.relpath
                except AttributeError: # this Movie was manually inited, not loaded from a textheader. s.moviepath doesn't exist, use s.path instead. Or it might not even be a Movie
                    pass
                '''
                # Also, if you inited a stim that needs to be loaded (like a movie), maybe you should also load it now (this wasn't done when execing the textheader)
                try:
                    s.load()
                except:
                    pass
                '''
            # Generate the sweeptable here, no need to load if from files anymore...
            #self.sweeptable = {[]} # dictionary of lists, ie sweeptable={'ori',[0,45,90],'sfreq',[1,1,1]}
            # so you index into it with self.sweeptable['var'][sweepi]
            # vars = self.sweeptable.keys()

            if self.stims: # this Experiment has object-oriented stim(s)
                for s in self.stims:
                    varvals={} # init a dictionary that will contain variable values
                    for var in s.varlist:
                        varvals[var]=eval('s.'+var) # generate a dictionary with var:val entries, to pass to buildSweepTable
                    s.sweepTable = buildSweepTable(s.varlist, varvals, s.nruns, s.shuffleRuns, s.blankSweep, s.shuffleBlankSweeps, makeSweepTableText=0)[0] # passing varlist by reference, dim indices end up being modified
            else: # this is a simple stim (not object oriented)
                varvals={} # init a dictionary that will contain variable values
                for var in self.varlist:
                    varvals[var]=eval('self.'+var) # generate a dictionary with var:val entries, to pass to buildSweepTable
                self.sweepTable = buildSweepTable(self.varlist, varvals, self.nruns, self.shuffleRuns, self.blankSweep, self.shuffleBlankSweeps, makeSweepTableText=0)[0] # passing varlist by reference, dim indices end up being modified

            '''
            # Old code for creating a sweeptable file (used by Matlab and NVS):
            sweeptabletext = sweeptabletext.replace('[','') # get rid of brackets and ' in first line, these demarcate dimensions, but aren't needed in matlab
            sweeptabletext = sweeptabletext.replace(']','')
            sweeptabletext = sweeptabletext.replace('\'','')
            sweeptabletext = sweeptabletext.replace(',','') # also, replace any commas or spaces (in the first line) with tabs for delimiting
            sweeptabletext = sweeptabletext.replace(' ','\t')

            fname = string.replace(sys.argv[0],sys.path[0]+'\\','') # name of file that launched Python
            fname = fname.replace('.textheader','') # remove .textheader part of filename
            #fname = fname.replace('.py','') # remove .py part of filename
            fname += '.sweeptable' # add .sweeptable extension

            fullpathfname = sys.path[0]+'\\'+fname

            print 'Writing to file:', fullpathfname

            f = file(fullpathfname,'w')
            f.write(sweeptabletext)
            f.close()
            '''

            for defaultm in [MSEQ32, MSEQ16]: # check if this Experiment uses specific default movies
                for s in self.stims: # for all stims inited by the textheader
                    if s.name == defaultm.name and isinstance(s, Movie):
                        if defaultm.data == None: # see if this default movie has yet to be loaded
                            defaultm.load() # load this default movie
                        s.data = defaultm.data # point this movie's data to default movie data
                        treestr = s.level*TAB + s.name
                        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

        try:
            self.REFRESHTIME = int(round(1/float(self.REFRESHRATE)*1000000)) # in us, keep 'em integers
        except AttributeError:
            self.REFRESHTIME = self.din[1,0] - self.din[0,0] # use the time difference between the first two din instead
        #self.buildsweepranges()

    def buildsweepranges(self):
        self.sweepranges = []

    def code(self, neuron, trange=None, **kwargs):
        """Returns a Neuron.Code object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        if trange == None:
            trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
        try:
            return neuron.code(trange=trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].code(trange=trange, **kwargs) # neuron is probably a Neuron id
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nNeuron.code: '+getargstr(Neuron.code)
    code.__doc__ += '\nbinary: '+getargstr(BinaryCode.__init__)

    DEFAULTCODEBITLENGTH = 10

    # All these spike code method thingies should really be more OOP. see 2006 Schneidmann figs 1e and 1f

    def codes(self, nis, **kwargs):
        """Returns a 2D array where each row is a neuron code constrained to the time range of this Experiment"""
        neurons = self.r.n
        codes = []
        for ni in nis:
            codes.append( [ self.code(neurons[ni], **kwargs).c ] ) # each is a nested list (ie, 2D)
        codes = tuple(codes) # required for concatenate
        return cat(codes)

    def intcodes(self, nis, **kwargs):
        """Returns an array of the integer representation of the neuronal population code for each time bin"""
        codes = self.codes(nis, **kwargs) # 2D array of binary words. rows are neuron codes, columns are words for each time bin
        return binaryarray2int(codes)

    def intcodesPDF(self, nis, **kwargs):
        """Returns the pdf across all possible population code words"""
        intcodes = self.intcodes(nis=nis, **kwargs)
        nbits = len(nis)
        p, bins = histogram(intcodes, bins=arange(2**nbits), normed='pmf')
        return p, bins

    def intcodesFPDF(self, nis, **kwargs):
        """Returns the probability of getting each population code word, assuming independence between neurons, taking into account each neuron's spike probability"""
        nbits = len(nis)
        intcodes = arange(2**nbits)
        codes = self.codes(nis=nis, **kwargs)
        spikeps = [] # list spike probabilities for all neurons
        for neuroncode in codes: # for each neuron, ie each row
            spikeps.append( neuroncode.sum() / float(neuroncode.size) ) # calc the average p of getting a spike for this neuron, within any time bin
        spikeps = array(spikeps, ndmin = 2)
        nospikeps = 1 - spikeps
        #print 'spikesps: ', spikeps
        #print 'nospikesps: ', nospikeps
        pon = getbinarytable(nbits)*spikeps.transpose()
        poff = (1 - getbinarytable(nbits))*nospikeps.transpose()
        #print 'pon', pon
        #print 'poff', poff
        x = pon + poff
        #print 'x', x
        intcodeps = x.prod(axis=0)
        return intcodeps, intcodes

    def plot_codeProbScatter(self, nbits=DEFAULTCODEBITLENGTH, randomneurons=False, **kwargs):
        """Scatterplots the expected probabilities of all possible population codes (y axis) vs their observed probabilities (x axis)"""
        # pick which and how many cells to include
        neurons = self.r.n
        nis = neurons.keys()
        if nbits == None: # use all cells
            nbits = len(nis)
        if randomneurons:
            nis = random.sample(nis, nbits) # randomly sample nbits of the nis
        else:
            nis.sort() # make sure they're in increasing order
            nis = nis[:nbits] # use just the first nbits neurons to make your words
        print 'neurons:', nis
        print 'try colouring each point in the scatter according to the number of 1s in the spike word'
        pobserved, observedwords = self.intcodesPDF(nis=nis, **kwargs)
        pexpected, expectedwords = self.intcodesFPDF(nis=nis, **kwargs) # expected, assuming independence
        figure()
        plot([10**-6, 1], [10**-6, 1], 'b-') # plot an x=y line
        hold(True)
        # pl.scatter(pobserved, pexpected), followed by setting the x and y axes to log scale freezes the figure and runs 100% cpu
        # gca().set_xscale('log')
        # gca().set_yscale('log')
        # use loglog() instead
        loglog(pobserved, pexpected, 'k.')
        xlabel('observed population code probability')
        ylabel('expected population code probability')
        title('population code probabilities - experiment %d - %s' % (self.id, self.name))

    def nspikingPDF(self, nbits=None, **kwargs):
        """Returns the PDF of observing n cells spiking in the same population code time bin"""
        print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
        pobserved, observedwords = self.intcodesPDF(nbits=nbits, **kwargs)
        nspiking = [] # collect observances of the number of cells spiking for each pop code time bin
        for observedword in observedwords: # slow hack
            nspiking.append( np.binary_repr(observedword).count('1') ) # convert words to binary, count the number of 1s in each
        pnspiking, bins = histogram(nspiking, bins=arange(nbits), normed='pmf') # histogram 'em
        return pnspiking, bins

    def nspikingFPDF(self, nbits=None, **kwargs):
        """Returns the PDF of observing n cells spiking in the same population code time bin, assuming independence by shuffling each cell's code train"""
        print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
        pass

    def plot_nspikingPDFandFPDF(self, nbits=None, **kwargs):
        print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
        nspikingPDF
        nspikingFPDF

















    class CodeWords(object):
        def __init__(self, trange=None):
            if trange == None:
                self.trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
            else:
                self.trange = trange
            self.neurons = self.r.n
        def calc(self):
            self.codes = []
        def codes(self, nis, **kwargs):
            """Returns a 2D array where each row is a neuron code constrained to the time range of this Experiment"""
            neurons = self.r.n
            codes = []
            for ni in nis:
                codes.append( [ self.code(neurons[ni], **kwargs).c ] ) # each is a nested list (ie, 2D)
            codes = tuple(codes) # required for concatenate
            return cat(codes)

        def intcodes(self, nis, **kwargs):
            """Returns an array of the integer representation of the neuronal population code for each time bin"""
            codes = self.codes(nis, **kwargs) # 2D array of binary words. rows are neuron codes, columns are words for each time bin
            return binaryarray2int(codes)

        def intcodesPDF(self, nis, **kwargs):
            """Returns the pdf across all possible population code words"""
            intcodes = self.intcodes(nis=nis, **kwargs)
            nbits = len(nis)
            p, bins = histogram(intcodes, bins=arange(2**nbits), normed='pmf')
            return p, bins

        def intcodesFPDF(self, nis, **kwargs):
            """Returns the probability of getting each population code word, assuming independence between neurons, taking into account each neuron's spike probability"""
            nbits = len(nis)
            intcodes = arange(2**nbits)
            codes = self.codes(nis=nis, **kwargs)
            spikeps = [] # list spike probabilities for all neurons
            for neuroncode in codes: # for each neuron, ie each row
                spikeps.append( neuroncode.sum() / float(neuroncode.size) ) # calc the average p of getting a spike for this neuron, within any time bin
            spikeps = array(spikeps, ndmin = 2)
            nospikeps = 1 - spikeps
            #print 'spikesps: ', spikeps
            #print 'nospikesps: ', nospikeps
            pon = getbinarytable(nbits)*spikeps.transpose()
            poff = (1 - getbinarytable(nbits))*nospikeps.transpose()
            #print 'pon', pon
            #print 'poff', poff
            x = pon + poff
            #print 'x', x
            intcodeps = x.prod(axis=0)
            return intcodeps, intcodes

        def plot_codeProbScatter(self, nbits=DEFAULTCODEBITLENGTH, randomneurons=False, **kwargs):
            """Scatterplots the expected probabilities of all possible population codes (y axis) vs their observed probabilities (x axis)"""
            # pick which and how many cells to include
            neurons = self.r.n
            nis = neurons.keys()
            if nbits == None: # use all cells
                nbits = len(nis)
            if randomneurons:
                nis = random.sample(nis, nbits) # randomly sample nbits of the nis
            else:
                nis.sort() # make sure they're in increasing order
                nis = nis[:nbits] # use just the first nbits neurons to make your words
            print 'neurons:', nis
            print 'try colouring each point in the scatter according to the number of 1s in the spike word'
            pobserved, observedwords = self.intcodesPDF(nis=nis, **kwargs)
            pexpected, expectedwords = self.intcodesFPDF(nis=nis, **kwargs) # expected, assuming independence
            figure()
            plot([10**-6, 1], [10**-6, 1], 'b-') # plot an x=y line
            hold(True)
            # pl.scatter(pobserved, pexpected), followed by setting the x and y axes to log scale freezes the figure and runs 100% cpu
            # gca().set_xscale('log')
            # gca().set_yscale('log')
            # use loglog() instead
            loglog(pobserved, pexpected, 'k.')
            xlabel('observed population code probability')
            ylabel('expected population code probability')
            title('population code probabilities - experiment %d - %s' % (self.id, self.name))

        def nspikingPDF(self, nbits=None, **kwargs):
            """Returns the PDF of observing n cells spiking in the same population code time bin"""
            print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
            pobserved, observedwords = self.intcodesPDF(nbits=nbits, **kwargs)
            nspiking = [] # collect observances of the number of cells spiking for each pop code time bin
            for observedword in observedwords: # slow hack
                nspiking.append( np.binary_repr(observedword).count('1') ) # convert words to binary, count the number of 1s in each
            pnspiking, bins = histogram(nspiking, bins=arange(nbits), normed='pmf') # histogram 'em
            return pnspiking, bins

        def nspikingFPDF(self, nbits=None, **kwargs):
            """Returns the PDF of observing n cells spiking in the same population code time bin, assuming independence by shuffling each cell's code train"""
            print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
            pass

        def plot_nspikingPDFandFPDF(self, nbits=None, **kwargs):
            print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
            nspikingPDF
            nspikingFPDF

        class plot(object):
            def scatter():
                pass
            def pdf():
                pass
















    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code objects
        Uses naive corr() f'n defined by me. SLOWWWWWWWWWWWW!!!!!!!!!!!!!!!!!!!!!!!!!"""
        code1 = self.code(neuron1, **kwargs)
        code2 = self.code(neuron2, **kwargs)
        return corr(code1.c, code2.c)
    codecorr.__doc__ += '\n\n**kwargs:'
    codecorr.__doc__ += '\nNeuron.code: '+getargstr(Neuron.code)
    codecorr.__doc__ += '\nbinary: '+getargstr(BinaryCode.__init__)

    class CodeCorrPDF(object):
        """A PDF of the correlations of the codes of all cell pairs in this Experiment"""
        def __init__(self, experiment=None, crange=None, nbins=100, normed='pmf', **kwargs):
            self.e = experiment
            self.r = self.e.r
            if crange: # range of corrs for pdf to span
                self.crange = crange
            self.nbins = nbins
            self.normed = normed
            self.kwargs = kwargs
        def __eq__(self, other):
            selfd = self.__dict__.copy()
            otherd = other.__dict__.copy()
            # Delete their n and c attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
            [ d.__delitem__(key) for d in [selfd, otherd] for key in ['n', 'c'] if d.has_key(key) ]
            if type(self) == type(other) and selfd == otherd:
                return True
            else:
                return False
        def calc(self):
            neurons = self.r.n.keys()
            nneurons = len(neurons)
            corrs = [ self.e.codecorr(neurons[ni1], neurons[ni2], **self.kwargs) for ni1 in range(0,nneurons) for ni2 in range(ni1,nneurons) if ni1 != ni2 ]
            try: # figure out the bin edges
                c = np.linspace(start=self.crange[0], stop=self.crange[1], num=self.nbins, endpoint=True)
            except AttributeError: # self.crange doesn't exist, let histogram() figure out the bin edges
                c = self.nbins
            self.n, self.c = histogram(corrs, bins=c, normed=self.normed)
        def plot(self):
            figure()
            try:
                barwidth = (self.crange[1] - self.crange[0]) / float(self.nbins)
            except AttributeError: # self.crange doesn't exist
                barwidth = self.c[1] - self.c[0]
            #hist(self.n, bins=self.r, normed=0, bottom=0, width=None, hold=False) # doesn't seem to work
            bar(left=self.c, height=self.n, width=barwidth, bottom=0, color='k', yerr=None, xerr=None, ecolor='k', capsize=3)
            try:
                xlim(self.crange)
            except AttributeError: # self.crange doesn't exist
                pass
            title('neuron pair code correlation pdf - experiment %d - %s' % (self.e.id, self.e.name))
            if self.normed:
                if self.normed == 'pmf':
                    ylabel('probability mass')
                else:
                    ylabel('probability density')
            else:
                ylabel('count')
            xlabel('correlation coefficient')

    def codecorrpdf(self, **kwargs):
        """Returns an existing CodeCorrPDF object, or creates a new one if necessary"""
        try:
            self._codecorrpdfs
        except AttributeError: # doesn't exist yet
            self._codecorrpdfs = [] # create a list that'll hold CodeCorrPDF objects
        co = self.CodeCorrPDF(experiment=self, **kwargs) # init a new one
        for ccpdf in self._codecorrpdfs:
            if co == ccpdf: # need to define special == method for class CodeCorrPDF()
                return ccpdf # returns the first object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codecorrpdfs
        co.calc() # no matching object was found, calculate it
        self._codecorrpdfs.append(co) # add it to the object list
        return co
    codecorrpdf.__doc__ += '\n\n**kwargs:'
    codecorrpdf.__doc__ += '\nCodeCorrPDF: '+getargstr(CodeCorrPDF.__init__)
    codecorrpdf.__doc__ += '\nNeuron.code: '+getargstr(Neuron.code)
    codecorrpdf.__doc__ += '\nbinary: '+getargstr(BinaryCode.__init__)

    def rate(self, neuron, **kwargs):
        """Returns a Neuron.Rate object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
        try:
            return neuron.rate(trange=trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].rate(trange=trange, **kwargs) # neuron is probably a Neuron id
    rate.__doc__ += '\n\n**kwargs:'
    rate.__doc__ += Neuron._rateargs

    def ratepdf(self, neuron, **kwargs):
        """Returns a Neuron.RatePDF object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
        try:
            return neuron.ratepdf(trange=trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].ratepdf(trange=trange, **kwargs) # neuron is probably a Neuron id
    ratepdf.__doc__ += '\n\n**kwargs:'
    ratepdf.__doc__ += '\nNeuron.RatePDF: '+getargstr(RatePDF.__init__)
    ratepdf.__doc__ += '\nNeuron.rate: '+getargstr(Neuron.rate)
    ratepdf.__doc__ += Neuron._rateargs
