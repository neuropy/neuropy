"""Defines the Recording class and all of its support classes"""

#print 'importing Recording'

from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class BaseRecording(object):
    """A Recording corresponds to a single SURF file, ie everything recorded between when
    the user hits record and when the user hits stop and closes the SURF file, including any
    pauses in between Experiments within that Recording. A Recording can have multiple Experiments,
    and multiple spike extractions, called Rips"""
    def __init__(self, id=None, name=None, parent=None):

        from Track import Track

        self.level = 3 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        if parent == None:
            try:
                self.t = _data.c[DEFAULTCATID].t[DEFAULTTRACKID] # see if the default Track has already been init'd
            except KeyError:
                self.t = Track() # init the default Track...
                _data.c[DEFAULTCATID].t[self.t.id] = self.t  # ...and add it to the default Cat object's list of Tracks
        else:
            self.t = parent # save parent Track object
        if id is not None:
            name = self.id2name(self.t.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'Recording id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.t.path + self.name + SLASH
        self.t.r[self.id] = self # add/overwrite this Recording to its parent's dict of Recordings, in case this Recording wasn't loaded by its parent
        self.e = {} # store Experiments in a dictionary
        self.rip = {} # store Rips in a dictionary
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.t.writetree(string)
    def id2name(self, path, id):
        if len(str(id)) == 1: # if id is only 1 digit long
            id = '0'+str(id) # add a leading zero
        name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith(str(id)+' - ') ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Recording id: %s' % id
        else:
            name = name[0] # pull the string out of the list
        return name
    def name2id(self, name):
        try:
            id = name[0:name.index(' - ')] # everything before the first ' - ', index() raises ValueError if it can't be found
        except ValueError:
            raise ValueError, 'Badly formatted Recording name: %s' % name
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, leave it as a string
        return id
    def load(self):

        from Experiment import Experiment
        from Rip import Rip

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.din') ] # returns din filenames without their .din extension
        for (experimentid, experimentName) in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
            experiment = Experiment(id=experimentid, name=experimentName, parent=self) # pass both the id and the name
            experiment.load() # load the Experiment
            self.e[experiment.id] = experiment # save it
        #if len(self.e) == 1:
        #   self.e = self.e.values[0] # pull it out of the dictionary
        ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
        defaultRipNames = [ ripName for ripName in ripNames for ripkeyword in RIPKEYWORDS if ripName.count(ripkeyword) ]
        if len(defaultRipNames) < 1:
            warn('Couldn\'t find a default Rip for Recording(%s)' % self.id)
        if len(defaultRipNames) > 1: # This could just be a warning instead of an exception, but really, some folder renaming is in order
            raise RuntimeError, 'More than one Rip folder in Recording(%s) has a default keyword: %s' %(self.id, defaultRipNames)
        for (ripid, ripName) in enumerate(ripNames): # ripids will be according to alphabetical order of ripNames
            rip = Rip(id=ripid, name=ripName, parent=self) # pass both the id and the name
            rip.load() # load the Rip
            self.rip[rip.name] = rip # save it
            # make the Neurons from the default Rip (if it exists in the Recording path) available in the Recording, so you can access them via r.n[nid] instead of having to do r.rip[name].n[nid]. Make them just another pointer to the data in r.rip[ripName].n
            for ripkeyword in RIPKEYWORDS[::-1]: # reverse the keywords so first one gets processed last
                if rip.name.count(ripkeyword): # if the keyword is in the ripName
                    self.n = self.rip[rip.name].n # make it the default Rip
                    self.cn = self.rip[rip.name].cn # make it the default Rip for ConstrainedNeurons too
        #if len(self.rip) == 1:
        #   self.rip = self.rip.values[0] # pull it out of the dictionary
        try:
            firstexp = min(self.e.keys())
            lastexp = max(self.e.keys())
            self.trange = self.e[firstexp].trange[0], self.e[lastexp].trange[1] # start of the first experiment to end of the last one
        except ValueError: # self.e is empty, no Experiments in this Recording, use first and last spike across all Neurons
            tranges = [ n.trange for n in self.n.values() ]
            self.trange = min(tranges[:][0]), max(tranges[:][1])
        # then, maybe add other info about the Recording, stored in the same folder, like skull coordinates, angles, polytrode name and type...

    '''
    # What should be done with this????????????

    def plot_codeProbScatter(self, nbits=DEFAULTCODEWORDLENGTH, randomneurons=False, **kwargs):
        """Scatterplots the expected probabilities of all possible population codes (y axis) vs their observed probabilities (x axis)"""
        neurons = self.n
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
        # set trange from first din of first experiment to last din of last experiment
        mint = np.inf
        maxt = 0
        for exp in self.e.values():
            if exp.din[0,0] < mint:
                mint = exp.din[0,0]
            if exp.din[-1,0] > maxt:
                maxt = exp.din[-1,0]
        trange = (mint, maxt+self.e[0].REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off
        # call the intcodes methods bound to the first (or any) experiment, cuz you have to have 'em bound.
        pobserved, observedwords = self.e[0].intcodesPDF(nis=nis, trange=trange, **kwargs)
        pexpected, expectedwords = self.e[0].intcodesFPDF(nis=nis, trange=trange, **kwargs) # expected, assuming independence
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
        title('population code probabilities - recording %d - %s' % (self.id, self.name))
    '''

class PopulationRaster(object):
    """A population spike raster plot. 'sortby' is the neuron attribute name to sort the raster by.
    Useful attributes to sort by: 'id', 'nspikes', 'trange'"""
    def __init__(self, recording=None, experiments=None, sortby='id'):
        self.r = recording
        if experiments == None:
            self.e = recording.e # dictionary
        else:
            self.e = experiments # should also be a dict
        firstexp = min(self.e.keys())
        self.t0 = self.e[firstexp].trange[0]
        experimentmarkers = [] # a list of all experiment start and stop times, in sorted order
        for e in self.e.values():
            experimentmarkers.extend(e.trange)
        self.experimentmarkers = asarray(experimentmarkers) - self.t0 # make 'em relative to t0
        self.experimentmarkers.sort() # just in case exps weren't in sorted order for some reason
        self.sortby = sortby
        self.neurons = list(self.r.n.values()) # convert to a list to allow sorting
        self.sort()
        self.f = figure(figsize=(14, 6))
        self.a = self.f.add_subplot(111)
        self.a.xaxis.set_major_locator(neuropyAutoLocator()) # better behaved tick locator
        self.a.xaxis.set_major_formatter(neuropyScalarFormatter()) # better behaved tick label formatter
        gcfm().frame.SetTitle(lastcmd())
        #gcfm().frame.SetTitle('r%d.raster(sortby=%s)' % (self.r.id, repr(self.sortby)))
        self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100)) # create a long tooltip with newline to get around bug where newlines aren't recognized on subsequent self.tooltip.SetTip() calls
        self.tooltip.Enable(False) # leave disabled for now
        self.tooltip.SetDelay(0) # set popup delay in ms
        gcfm().canvas.SetToolTip(self.tooltip) # connect the tooltip to the canvas
        self.a.set_xlabel('time (msec)')
        self.a.set_yticks([]) # turn off y axis
        self.yrange = (0, len(self.neurons))
        self.a.set_ylim(self.yrange)
        self.a.set_position([0.02, 0.1, 0.96, 0.88])
        self.f.canvas.mpl_connect('motion_notify_event', self.onmotion)
        self.f.canvas.mpl_connect('key_press_event', self.onkeypress)
    def sort(self):
        """Sorts self.neurons according to the neuron attribute specified by self.sortby"""
        if self.sortby != None:
            self.neurons.sort(key=lambda n: n.__getattribute__(self.sortby))
            print 'sorted by %s: %s' % (self.sortby, repr([ n.__getattribute__(self.sortby) for n in self.neurons ]))
    def plot(self, left=0, width=200000):
        """Plots the raster, units are us wrt beginning of first experiment"""
        self.left = left
        self.width = width
        # plot experiment start and endpoints
        for e in self.e.values():
            estart = e.trange[0]-self.t0
            eend = e.trange[1]-self.t0
            if left <=  estart and estart <= left+width: # experiment start point is within view
                startlines = self.a.vlines(x=estart/1000.0, ymin=self.yrange[0], ymax=self.yrange[1], fmt='k-') # marks exp start, convert to ms
                startlines[0].set_color((0, 1, 0)) # set to bright green
            if left <= eend and eend <= left+width: # experiment end point is within view
                endlines = self.a.vlines(x=eend/1000.0, ymin=self.yrange[0], ymax=self.yrange[1], fmt='k-') # marks exp end, convert to ms
                endlines[0].set_color((1, 0, 0)) # set to bright red
        # plot the rasters
        for nii, neuron in enumerate(self.neurons):
            x = (neuron.cut((self.t0+left, self.t0+left+width)) - self.t0) / 1000.0 # make spike times always relative to t0, convert to ms
            self.a.vlines(x=x, ymin=nii, ymax=nii+1, fmt='k-')
        self.a.set_xlim(left/1000.0, (left+width)/1000.0) # convert from us to ms
    def panx(self, npages=None, left=None):
        """Pans the raster along the x axis by npages, or to position left"""
        self.a.lines=[] # first, clear all the vlines, this is easy but a bit innefficient, since we'll be redrawing most of the ones we just cleared
        if left != None: # use left
            self.plot(left=left, width=self.width)
        else: # use npages instead
            self.plot(left=self.left+self.width*npages, width=self.width)
        self.f.canvas.draw() # redraw the figure
    def zoomx(self, factor):
        """Zooms the raster along the x axis by factor"""
        self.a.lines=[] # first, clear all the vlines, this is easy but a bit innefficient, since we'll be redrawing most of the ones we just cleared
        centre = (self.left + self.left+self.width) / 2.0
        width = self.width / factor
        left = centre - width / 2.0
        self.plot(left=left, width=width)
        self.f.canvas.draw() # redraw the figure
    def onmotion(self, event):
        """Called during mouse motion over figure. Pops up neuron and
        experiment info in a tooltip when hovering over a neuron row."""
        if event.xdata != None and event.ydata != None: # if mouse is inside the axes
            nii = int(math.floor(event.ydata)) # use ydata to get index into sorted list of neurons
            currentexp = None
            for e in self.e.values(): # for all experiments
                estart = (e.trange[0]-self.t0)/1000.0
                eend = (e.trange[1]-self.t0)/1000.0
                if estart < event.xdata  < eend:
                    currentexp = e
                    break # don't need to check any of the other experiments
            tip = 'n%d: %d spikes' % (self.neurons[nii].id, self.neurons[nii].nspikes)
            if currentexp == None:
                tip += '\nno experiment'
            else:
                tip += '\nexperiment %s: %s' % (currentexp.id, repr(currentexp.name))
            self.tooltip.SetTip(tip) # update the tooltip
            self.tooltip.Enable(True) # make sure it's enabled
        else: # mouse is outside the axes
            self.tooltip.Enable(False) # disable the tooltip
    def onkeypress(self, event):
        """Called during a figure keypress"""
        key = event.guiEvent.GetKeyCode() # wx dependent
        # you can also just use the backend-neutral event.key, but that doesn't recognize as many keypresses, like pgup, pgdn, etc.
        if not event.guiEvent.ControlDown(): # wx dependent
            if key == wx.WXK_RIGHT:
                self.panx(+0.1)
            elif key == wx.WXK_LEFT:
                self.panx(-0.1)
            elif key == wx.WXK_UP:
                self.zoomx(1.2)
            elif key == wx.WXK_DOWN:
                self.zoomx(1/1.2)
            elif key == wx.WXK_NEXT: # PGDN (page right)
                self.panx(+1)
            elif key == wx.WXK_PRIOR: # PGUP (page left)
                self.panx(-1)
            elif key == wx.WXK_HOME: # go to start of first Experiment
                self.panx(left=0)
            elif key == wx.WXK_END: # go to end of last Experiment
                lastexp = max(self.e.keys())
                self.panx(left=self.e[lastexp].trange[1]-self.t0-self.width)
        else: # Ctrl key is down, skip backwards or forwards to next experiment marker
            if key == wx.WXK_LEFT:
                i = self.experimentmarkers.searchsorted(self.left, side='left') # current position of left edge of the window in experimentmarkers list
                i = max(0, i-1) # decrement by 1, do bounds checking
                self.panx(left=self.experimentmarkers[i])
            elif key == wx.WXK_RIGHT:
                i = self.experimentmarkers.searchsorted(self.left, side='right') # current position of left edge of the window in experimentmarkers list
                i = min(i, len(self.experimentmarkers)-1) # bounds checking
                self.panx(left=self.experimentmarkers[i])
            elif key == wx.WXK_UP: # zoom in faster
                self.zoomx(3.0)
            elif key == wx.WXK_DOWN: # zoom out faster
                self.zoomx(1/3.0)


class RecordingRaster(BaseRecording):
    """Mix-in class that defines the raster related Recording methods"""
    def raster(self, **kwargs):
        """Creates a population spike raster plot"""
        sortby = kwargs.pop('sortby', 'id')
        pr = PopulationRaster(recording=self, sortby=sortby)
        pr.plot(**kwargs)
    raster.__doc__ += '\n\n'+PopulationRaster.__doc__
    raster.__doc__ += '\n\n**kwargs:'
    raster.__doc__ += '\n__init__: '+getargstr(PopulationRaster.__init__)
    raster.__doc__ += '\n    plot: '+getargstr(PopulationRaster.plot)


class Codes(object):
    """A 2D array where each row is a neuron code, and each column
    is a binary population word for that time bin"""
    def __init__(self, neurons=None, tranges=None, **kwargs):
        self.tranges = tranges
        self.neurons = neurons
        self.kwargs = kwargs
    def calc(self):
        self.c = []
        for neuron in self.neurons.values():
            self.c.append( [ neuron.code(tranges=self.tranges, **self.kwargs).c ] ) # each is a nested list (ie, 2D)
        self.c = tuple(self.c) # required for concatenate
        self.c = cat(self.c)
        #self.c.reshape(len(self.neurons), -1) # reshape to 2D array

'''
class CodeWords(object):
    """What's this supposed to do?"""
    pass
'''

class CodeCorrPDF(object):
    """A PDF of the correlations of the codes of all cell pairs in this Recording
    See 2006 Schneidman fig 1d"""
    def __init__(self, recording=None, experiments=None, **kwargs):
        self.r = recording
        if experiments != None:
            assert experiments.__class__ == dict
        self.e = experiments # save it, should be a dict if not None
        if self.e != None: # specific experiments were specified
            self.tranges = [ e.trange for e in self.e.values() ]
        else:
            self.tranges = [ self.r.trange ] # use the Recording's trange
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
        cnis = self.r.cn.keys() # ConstrainedNeuron indices
        ncneurons = len(cnis)
        self.corrs = [ self.r.codecorr(cnis[cnii1], cnis[cnii2], tranges=self.tranges, **self.kwargs)
                       for cnii1 in range(0,ncneurons) for cnii2 in range(cnii1,ncneurons)
                       if cnii1 != cnii2 ]
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
        gcfm().frame.SetTitle(lastcmd())
        #gcfm().frame.SetTitle('r%d.codecorrpdf(nbins=%d)' % (self.r.id, self.nbins))
        titlestring = 'neuron pair code correlation pdf'
        if self.e != None:
            print self.e
            titlestring += '\nexperiments: %s' % repr(self.e.keys())
        a.set_title(titlestring)
        if self.normed:
            if self.normed == 'pmf':
                a.set_ylabel('probability mass')
            else:
                a.set_ylabel('probability density')
        else:
            a.set_ylabel('count')
        a.set_xlabel('correlation coefficient')


class RecordingCode(BaseRecording):
    """Mix-in class that defines the spike code related Recording methods"""
    def code(self, cneuron=None, **kwargs):
        """Returns a ConstrainedNeuron.Code object, constrained to the time
        ranges of the Experiments in this Recording, as well as by tranges. Takes either a
        ConstrainedNeuron object or just a ConstrainedNeuron id"""
        try:
            return cneuron.code(**kwargs) # see if cneuron is a ConstrainedNeuron
        except AttributeError:
            return self.cn[cneuron].code(**kwargs) # cneuron is probably a ConstrainedNeuron id
    code.__doc__ += '\n\n**kwargs:'
    #code.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code) # causes import problems
    #code.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__) # causes import problems

    def codes(self, neurons=None, **kwargs):
        """Returns a 2D array where each row is a neuron code constrained to the time range of this Experiment
        INCOMPLETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
        if neurons == None:
            neurons = self.n
        codeso = Codes(neurons=neurons, tranges=self.tranges, **kwargs)
        codeso.calc()
        return codeso
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nCodes: '+getargstr(Codes.__init__)
    #code.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code) # causes import problems
    #code.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__) # causes import problems

    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code (or ConstrainedNeuron.Code)
        objects. Uses naive corrcoef() f'n defined by me. SLOWWWWWWWWWWWW!!!!!!!!!!!!!!!!!!!!!!!!!"""
        code1 = self.code(neuron1, **kwargs)
        code2 = self.code(neuron2, **kwargs)
        return corrcoef(code1.c, code2.c)
    codecorr.__doc__ += '\n\n**kwargs:'
    #codecorr.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code) # causes import problems
    #codecorr.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__) # causes import problems

    def codecorrpdf(self, experiments=None, **kwargs):
        """Returns an existing CodeCorrPDF object, or creates a new one if necessary"""
        try:
            self._codecorrpdfs
        except AttributeError: # doesn't exist yet
            self._codecorrpdfs = [] # create a list that'll hold CodeCorrPDF objects
        cco = CodeCorrPDF(recording=self, experiments=experiments, **kwargs) # init a new one
        for ccpdf in self._codecorrpdfs:
            if cco == ccpdf: # need to define special == method for class CodeCorrPDF()
                return ccpdf # returns the first object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codecorrpdfs
        cco.calc() # no matching object was found, calculate it
        self._codecorrpdfs.append(cco) # add it to the object list
        return cco
    codecorrpdf.__doc__ += '\n\n**kwargs:'
    codecorrpdf.__doc__ += '\nCodeCorrPDF: '+getargstr(CodeCorrPDF.__init__)
    #codecorrpdf.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code) # causes import problems
    #codecorrpdf.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__) # causes import problems

    '''
    def codewords(self, **kwargs):
        cw = CodeWords(tranges=self.tranges)
        cw.calc()
        return cw
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

    def scatter(self, nbits=DEFAULTCODEWORDLENGTH, randomneurons=False, shufflecodes=False, **kwargs):
        """Scatterplots the expected probabilities of all possible population codes (y axis) vs their observed probabilities (x axis)"""
        # pick which and how many cells to include
        nis = self.neurons.keys()
        if nbits == None: # use all cells
            nbits = len(nis)
        if randomneurons:
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
        gcfm().frame.SetTitle(lastcmd())
        #gcfm().frame.SetTitle('r%d.e[%d].schneidman.scatter(nbits=%s, randomneurons=%s, shufflecodes=%s)' % (self.e.r.id, self.e.id, nbits, randomneurons, shufflecodes))
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

class RecordingSchneidman(BaseRecording):
    """Mix-in class that defines the spike code related Schneidman methods"""
    def schneidman(self, experiments=None):
        """Returns a Schneidman object"""
        if experiments == None:
            experiments = self.e # all experiments in this Recording
        so = Schneidman(experiments=experiments)
        return so


class Recording(RecordingRaster,
                RecordingCode,
                RecordingSchneidman,
                BaseRecording):
    """Inherits all the Recording classes into a single Recording class"""
    pass
