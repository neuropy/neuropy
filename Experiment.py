"""Defines the Experiment class and all of its support classes."""

#print 'importing Experiment'

from Core import *
from Dimstim.Core import buildSweepTable
import Neuron
from Recording import PopulationRaster

class BaseExperiment(object):
    """An Experiment corresponds to a single contiguous VisionEgg stimulus session.
    It contains information about the stimulus during that session, including
    the DIN values, the text header, and any Movies that were involved"""

    from Recording import Recording

    def __init__(self, id=None, name=None, parent=Recording):
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.r = parent() # init parent Recording object
        except TypeError: # parent is an instance, not a class
            self.r = parent # save parent Recording object
        if name is None:
            raise ValueError, 'Experiment name can\'t be None'
        self.id = id # not really used by the Experiment class, just there for user's info
        self.name = name
        self.path = self.r.path
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)

    # doesn't need a id2name or name2id method, neither can really be derived from the other in an easy way (although could use re), the id is just chronological (which is also alphabetical) order, at least for now

    def load(self):

        from Recording import Recording
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
            exec(self.textheader) # execute the textheader as Python code. TODO: maybe add some checks here to prevent changes to filesystem from accidental code, unload the os or sys modules or something? But that wouldn't prevent 'accidental' code from 'accidentally' re-importing such modules
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
            # self.sweeptable = {[]} # dictionary of lists, ie sweeptable={'ori',[0,45,90],'sfreq',[1,1,1]}
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

        self.trange = (self.din[0,0], self.din[-1,0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off

    def buildsweepranges(self):
        print 'INCOMPLETE!!!!!!!!!!!!!!!!'
        self.sweepranges = {}


class CodeCorrPDF(object):
    """A PDF of the correlations of the codes of all cell pairs in this Experiment
    See 2006 Schneidman fig 1d"""
    def __init__(self, experiment=None, **kwargs):
        self.e = experiment
        self.r = self.e.r
        self.kwargs = kwargs
    def __eq__(self, other):
        selfd = self.__dict__.copy()
        otherd = other.__dict__.copy()
        # Delete their n and c attribs, if they exist, to prevent comparing them below, since those attribs may not have yet been calculated
        [ d.__delitem__(key) for d in [selfd, otherd] for key in ['corrs', 'n', 'c', 'crange', 'nbins', 'normed'] if d.has_key(key) ]
        if type(self) == type(other) and selfd == otherd:
            return True
        else:
            return False
    def calc(self):
        neurons = self.r.n.keys()
        nneurons = len(neurons)
        self.corrs = [ self.e.codecorr(neurons[ni1], neurons[ni2], **self.kwargs) for ni1 in range(0,nneurons) for ni2 in range(ni1,nneurons) if ni1 != ni2 ]
    def plot(self, figsize=(7.5, 6.5), crange=None, nbins=100, normed='pdf'):
        self.crange = crange
        self.nbins = nbins
        self.normed = normed
        f = figure(figsize=figsize)
        a = f.add_subplot(111)
        try: # figure out the bin edges
            c = np.linspace(start=self.crange[0], stop=self.crange[1], num=self.nbins, endpoint=True)
        except TypeError: # self.crange is None, let histogram() figure out the bin edges
            c = self.nbins
        self.n, self.c = histogram(self.corrs, bins=c, normed=self.normed)
        try:
            barwidth = (self.crange[1] - self.crange[0]) / float(self.nbins)
        except TypeError: # self.crange is None, take width of first bin in self.c
            barwidth = self.c[1] - self.c[0]
        a.bar(left=self.c, height=self.n, width=barwidth, bottom=0, color='k', yerr=None, xerr=None, ecolor='k', capsize=3)
        try:
            a.set_xlim(self.crange)
        except TypeError: # self.crange is None
            pass
        gcfm().frame.SetTitle('r%d.e[%d].codecorrpdf(nbins=%d)' % (self.r.id, self.e.id, self.nbins))
        a.set_title('neuron pair code correlation pdf')
        if self.normed:
            if self.normed == 'pmf':
                a.set_ylabel('probability mass')
            else:
                a.set_ylabel('probability density')
        else:
            a.set_ylabel('count')
        a.set_xlabel('correlation coefficient')


class Codes(object):
    """A 2D array where each row is a neuron code, and each column
    is a binary population word for that time bin"""
    def __init__(self, neurons=None, trange=None, **kwargs):
        self.trange = trange
        self.neurons = neurons
        self.kwargs = kwargs
    def calc(self):
        self.c = []
        for neuron in self.neurons.values():
            self.c.append( [ neuron.code(trange=self.trange, **self.kwargs).c ] ) # each is a nested list (ie, 2D)
        self.c = tuple(self.c) # required for concatenate
        self.c = cat(self.c)
        #self.c.reshape(len(self.neurons), -1) # reshape to 2D array
'''
class CodeWords(object):
    """What's this supposed to do?"""
    pass
'''

class Schneidman(object):
    """see 2006 Schneidman figs 1e and 1f"""
    def __init__(self, experiment=None):
        self.e = experiment
        self.neurons = self.e.r.n

    def intcodes(self, nis=None, **kwargs):
        """Returns an array of the integer representation of the neuronal population code for each time bin"""
        if nis == None:
            nis = self.neurons.keys()
        neurons = {}
        [ neurons.__setitem__(ni, self.neurons[ni]) for ni in nis ]
        codeso = self.e.codes(neurons=neurons, **kwargs) # 2D array of binary words. rows are neuron codes, columns are words for each time bin
        return binaryarray2int(codeso.c)

    def intcodesPDF(self, nis=None, **kwargs):
        """Returns the pdf across all possible population code words"""
        if nis == None:
            nis = self.neurons.keys()
        intcodes = self.intcodes(nis=nis, **kwargs)
        nbits = len(nis)
        p, bins = histogram(intcodes, bins=arange(2**nbits), normed='pmf')
        return p, bins

    def intcodesFPDF(self, nis=None, **kwargs):
        """Returns the probability of getting each population code word, assuming independence between neurons, taking into account each neuron's spike probability"""
        if nis == None:
            nis = self.neurons.keys()
        nbits = len(nis)
        intcodes = arange(2**nbits)
        neurons = {}
        [ neurons.__setitem__(ni, self.neurons[ni]) for ni in nis ]
        codeso = self.e.codes(neurons=neurons, **kwargs)
        spikeps = [] # list spike probabilities for all neurons
        for neuroncode in codeso.c: # for each neuron, ie each row
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

    def scatter(self, nbits=DEFAULTCODEWORDLENGTH, randomize=False, **kwargs):
        """Scatterplots the expected probabilities of all possible population codes (y axis) vs their observed probabilities (x axis)"""
        # pick which and how many cells to include
        nis = self.neurons.keys()
        if nbits == None: # use all cells
            nbits = len(nis)
        if randomize:
            nis = random.sample(nis, nbits) # randomly sample nbits of the nis
        else:
            nis.sort() # make sure they're in increasing order
            nis = nis[:nbits] # use just the first nbits neurons to make your words
        print 'neurons:', nis
        pobserved, observedwords = self.intcodesPDF(nis=nis, **kwargs)
        pexpected, expectedwords = self.intcodesFPDF(nis=nis, **kwargs) # expected, assuming independence
        figure()
        plot([10**-6, 1], [10**-6, 1], 'b-') # plot an x=y line
        hold(True)
        # pl.scatter(pobserved, pexpected), followed by setting the x and y axes to log scale freezes the figure and runs 100% cpu
        # gca().set_xscale('log')
        # gca().set_yscale('log')
        # use loglog() instead
        inds = []
        for nspikes in range(0,5):
            inds.append([])
            [ inds[nspikes].append(i) for i in range(0,2**nbits) if bin(i).count('1') == nspikes ]
        pobserved1 = pobserved[inds[1]]; pexpected1 = pexpected[inds[1]]
        pobserved2 = pobserved[inds[2]]; pexpected2 = pexpected[inds[2]]
        pobserved3 = pobserved[inds[3]]; pexpected3 = pexpected[inds[3]]
        pobserved4 = pobserved[inds[4]]; pexpected4 = pexpected[inds[4]]
        pobserved[inds[1]], pexpected[inds[1]] = None, None # remove all these
        pobserved[inds[2]], pexpected[inds[2]] = None, None
        pobserved[inds[3]], pexpected[inds[3]] = None, None
        pobserved[inds[4]], pexpected[inds[4]] = None, None
        loglog(pobserved, pexpected, 'k.') # plots what's left in black
        loglog(pobserved4, pexpected4, 'm.')
        loglog(pobserved3, pexpected3, 'c.')
        loglog(pobserved2, pexpected2, 'y.')
        loglog(pobserved1, pexpected1, 'r.')
        gcfm().frame.SetTitle('r%d.e[%d].schneidman.scatter(nbits=%s, randomize=%s)' % (self.e.r.id, self.e.id, nbits, randomize))
        title('neurons: %s' % repr(nis))
        xlabel('observed population code probability')
        ylabel('expected population code probability')

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

    def plot_pdf(self, nbits=None, **kwargs):
        """Plots nspikingPDF and nspikingFPDF together. See 2006 Schneidman fig 1e"""
        print 'INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!'
        nspikingPDF
        nspikingFPDF
    '''
    class plot(object):
        """Would allow you to do r92.e[0].schneidman().plot().scatter()"""
        def scatter():
            pass
        def pdf():
            pass
    '''

class ExperimentCode(BaseExperiment):
    """Mix-in class that defines the spike code related Experiment methods"""
    def code(self, neuron=None, **kwargs):
        """Returns a Neuron.Code object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        try:
            return neuron.code(trange=self.trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].code(trange=self.trange, **kwargs) # neuron is probably a Neuron id
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code)
    code.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__)

    def codes(self, neurons=None, **kwargs):
        """Returns a 2D array where each row is a neuron code constrained to the time range of this Experiment
        INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
        if neurons == None:
            neurons = self.r.n
        codeso = Codes(neurons=neurons, trange=self.trange, **kwargs)
        codeso.calc()
        return codeso
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nCodes: '+getargstr(Codes.__init__)
    code.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code)
    code.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__)

    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code objects
        Uses naive corrcoef() f'n defined by me. SLOWWWWWWWWWWWW!!!!!!!!!!!!!!!!!!!!!!!!!"""
        code1 = self.code(neuron1, **kwargs)
        code2 = self.code(neuron2, **kwargs)
        return corrcoef(code1.c, code2.c)
    codecorr.__doc__ += '\n\n**kwargs:'
    codecorr.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code)
    codecorr.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__)

    def codecorrpdf(self, **kwargs):
        """Returns an existing CodeCorrPDF object, or creates a new one if necessary"""
        try:
            self._codecorrpdfs
        except AttributeError: # doesn't exist yet
            self._codecorrpdfs = [] # create a list that'll hold CodeCorrPDF objects
        cco = CodeCorrPDF(experiment=self, **kwargs) # init a new one
        for ccpdf in self._codecorrpdfs:
            if cco == ccpdf: # need to define special == method for class CodeCorrPDF()
                return ccpdf # returns the first object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codecorrpdfs
        cco.calc() # no matching object was found, calculate it
        self._codecorrpdfs.append(cco) # add it to the object list
        return cco
    codecorrpdf.__doc__ += '\n\n**kwargs:'
    codecorrpdf.__doc__ += '\nCodeCorrPDF: '+getargstr(CodeCorrPDF.__init__)
    codecorrpdf.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code)
    codecorrpdf.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__)

    def schneidman(self):
        """Returns a Schneidman object"""
        so = Schneidman(experiment=self)
        return so

    '''
    def codewords(self, **kwargs):
        cw = CodeWords(trange=self.trange)
        cw.calc()
        return cw
    '''

class ExperimentRate(BaseExperiment):
    """Mix-in class that defines the spike rate related Experiment methods"""
    def rate(self, neuron, **kwargs):
        """Returns a Neuron.Rate object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        try:
            return neuron.rate(trange=self.trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].rate(trange=self.trange, **kwargs) # neuron is probably a Neuron id
    rate.__doc__ += '\n\n**kwargs:'
    rate.__doc__ += Neuron.Neuron._rateargs

    def ratepdf(self, neuron, **kwargs):
        """Returns a Neuron.RatePDF object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        try:
            return neuron.ratepdf(trange=self.trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].ratepdf(trange=self.trange, **kwargs) # neuron is probably a Neuron id
    ratepdf.__doc__ += '\n\n**kwargs:'
    ratepdf.__doc__ += '\nNeuron.RatePDF: '+getargstr(Neuron.RatePDF.__init__)
    ratepdf.__doc__ += '\nNeuron.rate: '+getargstr(Neuron.Neuron.rate)
    ratepdf.__doc__ += Neuron.Neuron._rateargs


class RevCorrs(object):
    """Base class for doing reverse correlation of multiple neurons to a simulus"""
    def __init__(self, neurons=None, experiment=None, trange=None, nt=10):
        self.neurons = neurons
        self.experiment = experiment
        if trange == None:
            self.trange = self.experiment.trange
        else:
            self.trange = trange
        self.movie = self.experiment.stims[0]
        self.nt = nt # number of revcorr timepoints
        self.tis = range(0, nt, 1) # revcorr timepoint indices
        self.t = [ int(round(ti * self.movie.sweeptimeMsec)) for ti in self.tis ] #list(array(self.tis) * self.movie.sweeptimeMsec) # revcorr timepoint values,
    def plot(self, interp='nearest', normed=True, title='ReceptiveFieldFrame', scale=2.0, **kwargs):
        """Plots the RFs as bitmaps in a wx.Frame. normed = 'global'|True|False"""
        rfs = [] # list of receptive fields to pass to ReceptiveFieldFrame object
        if normed == 'global': # normalize across all timepoints for all neurons
            vmin = min([ sta.rf.min() for sta in self.stas ]) # global min
            vmax = max([ sta.rf.max() for sta in self.stas ]) # global max
        for ni, sta in enumerate(self.stas):
            rf = sta.rf.copy() # create a copy to manipulate for display purposes, (nt, width, height)
            if normed: # either 'global' or True
                if normed == True: # normalize across the timepoints for this Neuron
                    vmin, vmax = rf.min(), rf.max()
                norm = mpl.colors.normalize(vmin=vmin, vmax=vmax, clip=True) # create a single normalization object to map luminance
                rf = norm(rf) # normalize the rf the same way across all timepoints
            else: # don't normalize across timepoints, leave each one to autoscale
                for ti in range(self.nt):
                    norm = mpl.colors.normalize(vmin=None, vmax=None, clip=True) # create a normalization object to map luminance to the range [0,1], autoscale
                    rf[ti] = norm(rf[ti]) # normalize the rf separately at each timepoint
            cmap = mpl.cm.jet # get a colormap object
            rf = cmap(rf)[::,::,::,0:3] # convert normalized luminance to RGB via the colormap, throw away alpha channel (not used for now in ReceptiveFieldFrame)
            rf = rf * 255 # scale up to 8 bit values
            rf = rf.round().astype(np.uint8) # downcast from float to uint8 for feeding to ReceptiveFieldFrame
            rfs.append(rf)
        frame = ReceptiveFieldFrame(title=title, rfs=rfs, neurons=self.neurons, t=self.t, scale=scale, **kwargs)
        frame.Show()

class STAs(RevCorrs):
    """Just a container class for multiple Neuron.STA objects. The plot() method is unique
    though: it plots all the Neuron.STA objects in a single figure"""
    def calc(self):
        self.stas = [] # store STAs in a list
        for neuron in self.neurons:
            stao = neuron.sta(experiment=self.experiment, trange=self.trange, nt=self.nt)
            self.stas.append(stao)
    def plot(self, interp='nearest', normed=True, scale=2.0, **kwargs):
        super(STAs, self).plot(interp=interp, normed=normed,
                               title='STA: r[%d], e[%d], interp=%s, normed=%s, scale=%s' %
                               (self.experiment.r.id, self.experiment.id, repr(interp), repr(normed), repr(scale)),
                               scale=scale,
                               **kwargs)
    plot.__doc__ = RevCorrs.plot.__doc__


class STCs(RevCorrs):
    """Just a container class for multiple Neuron.STC objects. The plot() method is unique
    though: it plots all the Neuron.STC objects in a single figure"""
    def calc(self):
        self.stcs = [] # store STCs in a list
        for neuron in self.neurons:
            stco = neuron.stc(experiment=self.experiment, trange=self.trange, nt=self.nt)
            self.stcs.append(stco)
    def plot(self, interp='nearest', normed=True, scale=2.0, **kwargs):
        super(STCs, self).plot(interp=interp, normed=normed,
                               title='STC: r[%d], e[%d], interp=%s, normed=%s, scale=%s' %
                               (self.experiment.r.id, self.experiment.id, repr(interp), repr(normed), repr(scale)),
                               scale=scale,
                               **kwargs)
    plot.__doc__ = RevCorrs.plot.__doc__


class ExperimentRevCorr(BaseExperiment):
    """Mix-in class that defines the reverse correlation related Experiment methods"""
    def sta(self, neurons=None, **kwargs):
        """Returns an STAs RevCorrs object and plots it"""
        if neurons == None: # no Neurons were passed, use all the Neurons from the default Rip for this experiment's Recording
            keyvals = self.r.n.items() # get key val pairs in a list of tuples
            keyvals.sort() # make sure they're sorted by key
            neurons = []
            for key, val in keyvals:
                neurons.append(val)
        else:
            try: # assume neurons is a list of Neuron ids, get the associated Neuron objects from the default Rip for this experiment's Recording
                neurons = [ self.r.n[ni] for ni in neurons ]
            except KeyError: # neurons is probably a list of Neuron objects
                pass
        staso = STAs(neurons=neurons, experiment=self, **kwargs) # init a new STAs object
        staso.calc()
        return staso
    sta.__doc__ += '\n\n**kwargs:\n'
    sta.__doc__ += getargstr(STAs.__init__)

    def stc(self, neurons=None, **kwargs):
        """Returns an STCs RevCorrs object and plots it"""
        if neurons == None: # no Neurons were passed, use all the Neurons from the default Rip for this experiment's Recording
            keyvals = self.r.n.items() # get key val pairs in a list of tuples
            keyvals.sort() # make sure they're sorted by key
            neurons = []
            for key, val in keyvals:
                neurons.append(val)
        else:
            try: # assume neurons is a list of Neuron ids, get the associated Neuron objects from the default Rip for this experiment's Recording
                neurons = [ self.r.n[ni] for ni in neurons ]
            except KeyError: # neurons is probably a list of Neuron objects
                pass
        stcso = STCs(neurons=neurons, experiment=self, **kwargs) # init a new STCs object
        stcso.calc()
        return stcso
    stc.__doc__ += '\n\n**kwargs:\n'
    stc.__doc__ += getargstr(STCs.__init__)


class ExperimentPopulationRaster(PopulationRaster):
    """A population raster limited to a single Experiment"""
    def __init__(self, experiment, sortby='id'):
        super(ExperimentPopulationRaster, self).__init__(recording=experiment.r, sortby=sortby)
        self.e = {experiment.id: experiment} # overwrite its e dictionary with just this experiment
        self.t0 = experiment.trange[0] # overwrite t0
        gcfm().frame.SetTitle('r%d.e[%d].raster(sortby=%s)' % (experiment.r.id, experiment.id, repr(self.sortby))) # set appropriate caption


class ExperimentRaster(BaseExperiment):
    """Mix-in class that defines the raster related Experiment methods"""
    def raster(self, **kwargs):
        """Creates a population spike raster plot"""
        sortby = kwargs.pop('sortby', 'id')
        pr = ExperimentPopulationRaster(experiment=self, sortby=sortby)
        pr.plot(**kwargs)
    raster.__doc__ += '\n\n'+ExperimentPopulationRaster.__doc__
    raster.__doc__ += '\n\n**kwargs:'
    raster.__doc__ += '\n__init__: '+getargstr(ExperimentPopulationRaster.__init__)
    raster.__doc__ += '\n    plot: '+getargstr(ExperimentPopulationRaster.plot)


class Experiment(ExperimentRaster,
                 ExperimentRevCorr,
                 ExperimentRate,
                 ExperimentCode,
                 BaseExperiment):
    """Inherits all the Experiment classes into a single Experiment class"""
    pass


