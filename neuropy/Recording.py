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
        self.treebuf = cStringIO.StringIO() # create a string buffer to print tree hierarchy to
        if parent == None:
            try:
                self.t = _data.a[ANIMALNAME].t[TRACKID] # see if the default Track has already been init'd
            except KeyError:
                self.t = Track() # init the default Track...
                _data.a[ANIMALNAME].t[self.t.id] = self.t  # ...and add it to the default Animal's list of Tracks
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
        self.path = os.path.join(self.t.path, self.name)
        self.t.r[self.id] = self # add/overwrite this Recording to its parent's dict of Recordings, in case this Recording wasn't loaded by its parent
        self.e = dictattr() # store Experiments in a dictionary with attrib access
        self.rip = dictattr() # store Rips in a dictionary with attrib access
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.t.writetree(string)
    def id2name(self, path, id):
        name = [ dirname for dirname in os.listdir(path)
                 if os.path.isdir(os.path.join(path, dirname))
                 and dirname.startswith(str(id)) ]
        if len(name) != 1:
            raise NameError, 'Ambiguous or non-existent Recording id: %s' % id
        else:
            return name[0] # pull the string out of the list
    def name2id(self, name):
        id = name.split()[0] # return the first word in the name, using whitespace as separators
        try:
            int(id[0]) # first character of Recording id better be an integer
        except ValueError:
            raise ValueError, 'First character of Recording name %r should be a number' % name
        try:
            id = int(id) # convert entire id to int if possible
        except ValueError:
            pass # it's alphanumeric (but starts with a number), leave it as a string
        return id
    def load(self):

        from Experiment import Experiment
        from Rip import Rip

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path)
                            if os.path.isfile(os.path.join(self.path, fname))
                            and fname.endswith('.din') ] # returns din filenames without their .din extension
        for experimentid, experimentName in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
            experiment = Experiment(id=experimentid, name=experimentName, parent=self) # pass both the id and the name
            experiment.load() # load the Experiment
            self.e[experiment.id] = experiment # save it
        ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path)
                     if os.path.isdir(os.path.join(self.path, dirname))
                     and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
        defaultRipNames = [ ripName for ripName in ripNames for ripkeyword in RIPKEYWORDS
                            if ripkeyword in ripName ]
        if len(defaultRipNames) < 1:
            warn('Couldn\'t find a default Rip for Recording(%s)' % self.id)
        if len(defaultRipNames) > 1: # This could just be a warning instead of an exception, but really, some folder renaming is in order
            raise RuntimeError, 'More than one Rip folder in Recording(%s) has a default keyword: %s' % (self.id, defaultRipNames)
        for ripid, ripName in enumerate(ripNames): # ripids will be according to alphabetical order of ripNames
            rip = Rip(id=ripid, name=ripName, parent=self) # pass both the id and the name
            rip.load() # load the Rip
            self.rip[rip.name] = rip # save it
            # make the Neurons from the default Rip (if it exists in the Recording path) available in the Recording, so you can access them via r.n[nid] instead of having to do r.rip[name].n[nid]. Make them just another pointer to the data in r.rip[ripName].n
            for ripkeyword in RIPKEYWORDS[::-1]: # reverse the keywords so last one gets processed first
                if ripkeyword in rip.name:
                    self.n = self.rip[rip.name].n # make it the default Rip
                    self.cn = self.rip[rip.name].cn # make it the default Rip for ConstrainedNeurons too
        if len(self.rip) == 1: # there's only one rip, make it the default Rip, even if it doesn't have a ripkeyword in it
            self.n = self.rip.values()[0].n # make it the default Rip
            self.cn = self.rip.values()[0].cn # make it the default Rip for ConstrainedNeurons too
        try:
            firstexp = min(self.e.keys())
            lastexp = max(self.e.keys())
            self.trange = self.e[firstexp].trange[0], self.e[lastexp].trange[1] # start of the first experiment to end of the last one
        except ValueError: # self.e is empty, no Experiments in this Recording, use first and last spike across all Neurons
            tranges = [ n.trange for n in self.n.values() ]
            self.trange = min(tranges[:][0]), max(tranges[:][1])
        # then, maybe add other info about the Recording, stored in the same folder, like skull coordinates, angles, polytrode name and type...


class PopulationRaster(object):
    """A population spike raster plot. nis are indices of neurons to
    plot in the raster, in order from bottom to top.
    jumpts are a sequence of timepoints (in us) that can then be quickly cycled
    through in the plot using keyboard controls.
    Defaults to absolute time origin (when acquisition began)"""
    def __init__(self, recording=None, experiments=None, nis=None,
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
        self.experimentmarkers = asarray(experimentmarkers) - self.t0 # make 'em relative to t0
        self.experimentmarkers.sort() # just in case exps weren't in sorted order for some reason

        self.jumpts = asarray(jumpts)
        self.jumpts.sort() # just in case jumpts weren't in sorted order for some reason

        units2tconv = {'usec': 1e0, 'msec': 1e3, 'sec': 1e6}
        self.units = units
        self.tconv = units2tconv[units]

        self.publication = publication

        self.neurons = self.r.n # still a dictattr
        if nis != None:
            self.nis = nis
        else:
            self.nis = self.r.n.keys()
            self.nis.sort() # keep it tidy

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
            figheight = 1.25+0.2*len(self.nis)
            self.f = figure(figsize=(14, figheight))
            self.a = self.f.add_subplot(111)
            self.formatter = neuropyScalarFormatter() # better behaved tick label formatter
            self.formatter.thousandsSep = ',' # use a thousands separator
            if not self.publication:
                self.a.xaxis.set_major_locator(neuropyAutoLocator()) # better behaved tick locator
                self.a.xaxis.set_major_formatter(self.formatter)
                self.a.set_yticks([]) # turn off y axis
            gcfm().frame.SetTitle(lastcmd())
            self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100)) # create a long tooltip with newline to get around bug where newlines aren't recognized on subsequent self.tooltip.SetTip() calls
            self.tooltip.Enable(False) # leave disabled for now
            self.tooltip.SetDelay(0) # set popup delay in ms
            gcfm().canvas.SetToolTip(self.tooltip) # connect the tooltip to the canvas
            self.a.set_xlabel('time (%s)' % self.units)
            if not self.publication:
                self.yrange = (0, len(self.nis))
            else:
                self.a.set_ylabel('cell index') # not really cell id, it's the ii (the index into the id)
                self.yrange = (-0.5, len(self.nis)-0.5)
            self.a.set_ylim(self.yrange)
            #aheight = min(0.025*len(self.nis), 1.0)
            bottominches = 0.75
            heightinches = 0.15+0.2*len(self.nis)
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
                startlines = self.a.vlines(x=estart/self.tconv, ymin=self.yrange[0], ymax=self.yrange[1], fmt='k-') # marks exp start, convert to ms
                startlines[0].set_color((0, 1, 0)) # set to bright green
            if left <= eend and eend <= left+width: # experiment end point is within view
                endlines = self.a.vlines(x=eend/self.tconv, ymin=self.yrange[0], ymax=self.yrange[1], fmt='k-') # marks exp end, convert t
                endlines[0].set_color((1, 0, 0)) # set to bright red
        # plot the bin edges. Not taking into account self.t0 for now, assuming it's 0
        if self.plotbinedges:
            leftbinedge = (left // self.binwidth + 1)*self.binwidth
            binedges = arange(leftbinedge, left+width, self.binwidth)
            binlines = self.a.vlines(x=binedges/self.tconv, ymin=self.yrange[0], ymax=self.yrange[1], fmt='b:') # convert t
        # plot the rasters
        for nii, ni in enumerate(self.nis):
            neuron = self.neurons[ni]
            x = (neuron.cut((self.t0+left, self.t0+left+width)) - self.t0) / self.tconv # make spike times always relative to t0, convert t
            if not self.publication:
                self.a.vlines(x=x, ymin=nii, ymax=nii+1, fmt='k-')
            else:
                self.a.vlines(x=x, ymin=nii-0.5, ymax=nii+0.5, fmt='k-')
        self.a.set_xlim(left/self.tconv, (left+width)/self.tconv) # convert t
    def _panx(self, npages=None, left=None):
        """Pans the raster along the x axis by npages, or to position left"""
        self.a.lines=[] # first, clear all the vlines, this is easy but a bit innefficient, since we'll probably be redrawing most of the ones we just cleared
        if left != None: # use left
            self.plot(left=left, width=self.width)
        else: # use npages instead
            self.plot(left=self.left+self.width*npages, width=self.width)
        self.f.canvas.draw() # redraw the figure
    def _zoomx(self, factor):
        """Zooms the raster along the x axis by factor"""
        self.a.lines=[] # first, clear all the vlines, this is easy but a bit innefficient, since we'll probably be redrawing most of the ones we just cleared
        centre = (self.left + self.left+self.width) / 2.0
        width = self.width / float(factor)
        left = centre - width / 2.0
        self.plot(left=left, width=width)
        self.f.canvas.draw() # redraw the figure
    def _goto(self):
        """Bring up a dialog box to jump to timepoint, mark it with a dotted line"""
        ted = wx.TextEntryDialog(parent=None, message='Go to timepoint (ms):', caption='Goto',
                                 defaultValue=str(int(round(self.left / self.tconv))), #wx.EmptyString,
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
        if event.xdata != None and event.ydata != None: # if mouse is inside the axes
            nii = int(math.floor(event.ydata)) # use ydata to get index into sorted list of neurons
            ni = self.nis[nii]
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

class RecordingRaster(BaseRecording):
    """Mix-in class that defines the raster related Recording methods"""
    def raster(self, **kwargs):
        """Creates a population spike raster plot"""
        #sortby = kwargs.pop('sortby', 'id')
        #pr = PopulationRaster(recording=self, sortby=sortby)
        pr = PopulationRaster(recording=self, **kwargs)
        return pr
    raster.__doc__ += '\n\n'+PopulationRaster.__doc__
    raster.__doc__ += '\n\n**kwargs:'
    raster.__doc__ += '\n__init__: '+getargstr(PopulationRaster.__init__)


class Codes(object):
    """A 2D array where each row is a neuron code, and each column
    is a binary population word for that time bin, sorted LSB to MSB from top to bottom.
    neurons is a list of Neurons, also from LSB to MSB. Order in neurons is preserved."""
    def __init__(self, neurons=None, kind=CODEKIND, tranges=None, tres=CODETRES, phase=0, shufflecodes=False):
        self.neurons = neurons
        self.kind = kind
        self.tranges = tolist(tranges)
        self.tres = tres
        self.phase = phase
        self.shufflecodes = shufflecodes
        self.nis = [ neuron.id for neuron in self.neurons ]
        self.nneurons = len(self.neurons)
        self.nis2niisdict = dict(zip(self.nis, range(self.nneurons))) # make a dict from keys:self.nis, vals:range(self.nneurons). This converts from nis to niis (from neuron indices to indices into the binary code array self.c)
    def nis2niis(self, nis=None):
        """Converts from nis to niis (from neuron indices to indices into the binary code array self.co.c).
        nis can be a sequence"""
        try:
            return [ self.nis2niisdict[ni] for ni in nis ]
        except TypeError: # iteration over non-sequence, nis is a scalar
            return self.nis2niisdict[nis]
    def calc(self):
        self.s = [] # stores the corresponding spike times for each neuron, just for reference
        self.c = [] # stores the 2D code array
        # append neurons in their order in self.neurons, store them LSB to MSB from top to bottom
        for neuron in self.neurons:
            codeo = neuron.code(kind=self.kind, tranges=self.tranges, tres=self.tres, phase=self.phase)
            self.s.append(codeo.s) # each is a nested list (ie, 2D), each row will have different length
            if self.shufflecodes:
                c = codeo.c.copy() # make a copy (wanna leave the codeo's codetrain untouched)
                np.random.shuffle(c) # shuffle each neuron's codetrain separately, in-place operation
            else:
                c = codeo.c # just a pointer
            self.c.append(c) # flat list
        self.t = codeo.t # stores the bin edges, just for reference. all bin times should be the same for all neurons, cuz they're all given the same trange. use the bin times of the last neuron
        nneurons = len(self.neurons)
        nbins = len(self.c[0]) # all entries in the list should be the same length
        self.c = cat(self.c).reshape(nneurons, nbins)
    def syncis(self):
        """Returns synch indices, ie the indices of the bins for which all the
        neurons in this Codes object have a 1 in them"""
        return self.c.prod(axis=0).nonzero()[0] # take product down all rows, only synchronous events across all cells will survive
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
        self.tranges = [ trange for codeso in codesos for trange in codeso.tranges ] # this tranges potentially holds multiple tranges from each codes objects, times the number of codes objects
        self.calc() # recalculate this code with its new set of tranges
    '''

class CodeCorrPDF(object):
    """A PDF of the correlations of the codes of all cell pairs in this Recording
    See 2006 Schneidman fig 1d"""
    def __init__(self, recording=None, experiments=None, **kwargs):
        self.r = recording
        if experiments != None:
            assert experiments.__class__ == dictattr
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
        [ d.__delitem__(key) for d in [selfd, otherd]
          for key in ['corrs', 'n', 'c', 'crange', 'nbins', 'normed']
          if d.has_key(key) ]
        if self.__class__ == other.__class__ and selfd == otherd:
            return True
        else:
            return False
    def calc(self):
        """Works on ConstrainedNeurons, but is constrained even further if experiments
        were passed and their tranges were used to generate self.tranges (see __init__)"""
        cnis = self.r.cn.keys() # ConstrainedNeuron indices
        ncneurons = len(cnis)
        # it's more efficient to precalculate the means and stds of each cell's codetrain,
        # and then reuse them in calculating the correlation coefficients:
        means = dict( ( cni, self.r.code(cni, tranges=self.tranges, **self.kwargs).c.mean() ) for cni in cnis ) # store each code mean in a dict
        stds  = dict( ( cni, self.r.code(cni, tranges=self.tranges, **self.kwargs).c.std() ) for cni in cnis ) # store each code std in a dict
        self.corrs = []
        for cnii1 in range(ncneurons):
            for cnii2 in range(cnii1+1, ncneurons):
                cni1 = cnis[cnii1]
                cni2 = cnis[cnii2]
                code1 = self.r.code(cni1, tranges=self.tranges, **self.kwargs).c
                code2 = self.r.code(cni2, tranges=self.tranges, **self.kwargs).c
                cc = ((code1 * code2).mean() - means[cni1] * means[cni2]) / (stds[cni1] * stds[cni2]) # (mean of product - product of means) / by product of stds
                self.corrs.append(cc)
        '''
        # simpler, but slower way:
        self.corrs = [ self.r.codecorr(cnis[cnii1], cnis[cnii2], tranges=self.tranges, **self.kwargs)
                       for cnii1 in range(0,ncneurons) for cnii2 in range(cnii1+1,ncneurons) ]
        '''
    def plot(self, figsize=(7.5, 6.5), crange=[-0.05, 0.3], nbins=40, normed='pdf'):
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
        titlestring += '\n%s' % lastcmd()
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
    def code(self, cneuron=None, kind=CODEKIND, tranges=None, tres=CODETRES, phase=0):
        """Returns a ConstrainedNeuron.Code object, constrained to the time
        ranges of the Experiments in this Recording, as well as by tranges. Takes either a
        ConstrainedNeuron object or just a ConstrainedNeuron id"""
        try:
            return cneuron.code(kind=kind, tranges=tranges, tres=tres, phase=phase) # see if cneuron is a ConstrainedNeuron
        except AttributeError:
            return self.cn[cneuron].code(kind=kind, tranges=tranges, tres=tres, phase=phase) # cneuron is probably a ConstrainedNeuron id

    def codes(self, neurons=None, experiments=None, kind=CODEKIND, tres=CODETRES, phase=0, shufflecodes=False):
        """Returns a Codes object, a 2D array where each row is a neuron code constrained to the time range of this Recording,
        or if specified, to the time ranges of Experiments in this Recording"""
        if neurons != None:
            if neurons.__class__ == list:
                try: # assume is a list of Neuron ids?
                    neurons = [ self.n[ni] for ni in neurons ] # build up list of Neurons, ordered according to the id list in neurons
                except: # assume is a list of Neurons
                    pass
            else: # assume neurons is a dict of neurons
                neurons = list(neurons.values()) # convert to list of Neurons
        else:
            neurons = list(self.n.values()) # list of Neurons
        if experiments != None:
            # need to preserve order of expids as specified
            if experiments.__class__ == list:
                try: # assume is a list of Experiment ids?
                    tranges = [ self.e[ei].trange for ei in experiments ]
                except: # assume is a list of Experiments
                    tranges = [ e.trange for e in experiments ]
            else: # assume experiments is a dict of Experiments
                tranges = [ e.trange for e in experiments.values() ]
        else: # no experiments specified, use whole Recording trange
            tranges = [self.trange]
        codeso = Codes(neurons=neurons, kind=kind, tranges=tranges, tres=tres, phase=phase, shufflecodes=shufflecodes)
        codeso.calc()
        return codeso
    codes.__doc__ += '\n\nCodes object:\n' + Codes.__doc__

    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code (or ConstrainedNeuron.Code)"""
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

class BaseNetstate(object):
    """Base class of Network state analyses.
    Implements a lot of the analyses on network states found in the 2006 Schneidman paper"""

    def __init__(self, recording, experiments=None, kind=CODEKIND, tres=CODETRES, phase=0):
        self.r = recording
        if experiments == None:
            self.tranges = [self.r.trange] # or should we check to see if this Recording has a tranges field due to appending Neurons?
        else:
            experiments = tolist(experiments)
            try:
                self.tranges = [ e.trange for e in experiments ] # is experiments a list of Experiments?
            except AttributeError:
                self.tranges = [ self.r.e[ei].trange for ei in experiments ] # assume experiments is a list of experiment ids
                experiments = [ self.r.e[ei] for ei in experiments ] # convert to a list of Experiments
        self.e = experiments # save list of Experiments (could potentially be None)
        self.kind = kind
        self.tres = tres
        self.phase = phase
        self.neurons = self.r.n
        self.nneurons = len(self.neurons)
        nis = self.neurons.keys() # get all neuron indices in this Recording
        nis.sort() # make sure they're sorted
        self.co = self.codes(nis=nis) # generate and save the Codes object for all the nis

    def codes(self, nis=None, shufflecodes=False):
        """Returns the appropriate Codes object, depending on the recording
        and experiments defined for this Netstate object"""
        cneurons = [ self.r.cn[ni] for ni in nis ] # build up list of ConstrainedNeurons, according to nis
        # get codes for this Recording constrained to when stimuli were on screen
        return self.r.codes(neurons=cneurons,
                            experiments=self.e,
                            kind=self.kind,
                            tres=self.tres,
                            phase=self.phase,
                            shufflecodes=shufflecodes)

    def wordts(self, nis=None, mis=None):
        """Returns word times, ie the times of the left bin edges for which all the
        neurons in the mis in this Netstate object have a 1 in them, and all
        the rest have a 0 in them. nis lists the total population of neuron ids"""
        if nis == None:
            nis = self.co.nis
        mis = toiter(mis)
        for mi in mis:
            assert mi in nis # make sure mis is a subset of nis
        co = self.codes(nis=nis) # make a new code object using the nis population
        nis2niis = co.nis2niis
        notmis = [ ni for ni in nis if ni not in mis ] # nis not in mis
        mis_high = co.c[nis2niis(mis)].prod(axis=0) == 1 # take product down all rows, only synchronous events across all mis cells will survive, boolean array
        notmis_low = co.c[nis2niis(notmis)].sum(axis=0) == 0 # boolean array
        i = (mis_high * notmis_low).nonzero()[0] # indices where mis are 1 and all the others are 0
        return co.t[i] # return the times at those indices

    def wordtsms(self, nis=None, mis=None):
        """Returns word times to the nearest msec, with the on bits specified in mis.
        nis lists the total population of neuron ids"""
        return np.int32(np.round(self.wordts(nis=nis, mis=mis) / 1e3))

    def intcodes(self, nis=None):
        """Given neuron indices (ordered LSB to MSB top to bottom), returns an array of the integer representation
        of the neuronal population binary code for each time bin"""
        assert self.kind == binary
        if nis == None:
            nis = random.sample(self.co.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
        return binarray2int(self.codes(nis=nis).c)

    def intcodesPDF(self, nis=None):
        """Returns the observed pdf across all possible population binary code words,
        labelled according to their integer representation"""
        if nis == None:
            nis = random.sample(self.co.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
        intcodes = self.intcodes(nis=nis)
        nbits = len(nis)
        p, bins = histogram(intcodes, bins=arange(2**nbits), normed='pmf')
        return p, bins

    def intcodesFPDF(self, nis=None):
        """the F stands for factorial. Returns the probability of getting each population binary code word, assuming
        independence between neurons, taking into account each neuron's spike (and no spike) probability"""
        if nis == None:
            nis = random.sample(self.co.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
        nbits = len(nis)
        intcodes = arange(2**nbits)
        #neurons = dict( (ni, self.neurons[ni]) for ni in nis ) # this is like dict comprehension, pretty awesome!
        codeso = self.codes(nis=nis)
        spikeps = [] # list spike probabilities for all neurons
        for neuroncode in codeso.c: # for each neuron, ie each row
            spikeps.append( neuroncode.mean() ) # calc the average p of getting a spike for this neuron, within any time bin.
        spikeps = array(spikeps, ndmin=2) # convert to an nbits*1 array, make sure it's explicitly treated as a 2D array that can be transposed, or something
        nospikeps = 1 - spikeps
        #print 'spikesps: ', spikeps.__repr__()
        #print 'nospikesps: ', nospikeps.__repr__()
        binarytable = getbinarytable(nbits)
        pon = binarytable * spikeps.transpose() # 2D array of probs of having a 1 in the right place for all possible population code words
        poff = (1 - binarytable) * nospikeps.transpose() # 2D array of probs of having a 0 in the right place for all possible population code words
        #print 'pon', pon.__repr__()
        #print 'poff', poff.__repr__()
        x = pon + poff # add the 2D arrays, each has zero prob values where the other has non-zero prob values
        #print 'x', x.__repr__()
        intcodeps = x.prod(axis=0) # take the product along the 0th axis (the columns) to get the prob of each population code word
        return intcodeps, intcodes

    def ising(self, nis=None, algorithm='CG'):
        """Returns an Ising maximum entropy model that takes into account pairwise correlations within neuron codes
        algorithm can be 'CG', 'BFGS', 'LBFGSB', 'Powell', or 'Nelder-Mead'"""
        if nis == None:
            nis = self.co.nis[0:CODEWORDLEN]
        #print 'nis:', nis.__repr__()
        codeso = self.codes(nis=nis)
        #c = codeso.c
        # convert values in codes object from [0, 1] to [-1, 1] by mutliplying by 2 and subtracting 1
        c = codeso.c.copy() # don't modify the original
        c = c*2 - 1 # this should be safe to do cuz c is a 2D array of signed int8 values
        #print 'c:', c.__repr__()
        means = [ row.mean() for row in c ] # iterate over rows of codes in c
        nrows = c.shape[0]
        pairmeans = [ (c[i]*c[j]).mean() for i in range(0, nrows) for j in range(i+1, nrows) ] # take a pair of rows, find the mean of their elementwise product
        isingo = Ising(means=means, pairmeans=pairmeans, algorithm=algorithm)
        return isingo

class NetstateIsingHist(BaseNetstate):
    """Netstate Ising parameter histograms. See Schneidman 2006 Fig 3b"""
    def calc(self, nbits=10, ngroups=5, algorithm='CG'):
        """Collects hi and Jij parameter values computed from ising models
        of ngroups subgroups of cells of size nbits"""
        self.nbits = nbits
        self.ngroups = ngroups
        self.algorithm = algorithm

        self.ims = [] # holds Ising Model objects
        self.his = []
        self.Jijs = []

        pd = wx.ProgressDialog(title='NetstateIsingHist progress', message='', maximum=ngroups, # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for groupi in range(self.ngroups): # for each group of nbits cells
            nis = random.sample(self.co.nis, nbits) # randomly sample nbits of the Netstate's nis attrib
            im = self.ising(nis=nis, algorithm=algorithm) # returns a maxent Ising model
            self.ims.append(im)
            self.his.append(im.hi)
            self.Jijs.append(im.Jij)
            cancel = not pd.Update(groupi, newmsg='groupi = %d' % groupi)
            if cancel:
                pd.Destroy()
                return
        pd.Destroy()

    def plot(self, nbins=50, hirange=(-2.5, 2.5), Jijrange=(-1.1, 1.1)):
        """Plots hi and Jij histograms in separate figures"""

        try: self.his, self.Jij
        except AttributeError: self.calc()

        # histogram them in linear space
        hibins  = np.linspace(start=hirange[0], stop=hirange[1], num=nbins, endpoint=True)
        Jijbins = np.linspace(start=Jijrange[0], stop=Jijrange[1], num=nbins, endpoint=True)
        nhi  = histogram(self.his, bins=hibins, normed='pdf')[0]
        nJij = histogram(self.Jijs, bins=Jijbins, normed='pdf')[0]

        # plot the hi histogram
        f1 = figure()
        a1 = f1.add_subplot(111)
        a1.hold(True)
        a1.bar(left=hibins, height=nhi, width=hibins[1]-hibins[0], color='g', edgecolor='g')
        gcfm().frame.SetTitle(lastcmd())
        a1.set_title('hi histogram\n%s' % lastcmd())
        a1.set_ylabel('probability density')
        a1.set_xlabel('hi')
        a1.set_xlim(hirange)

        # plot the Jij histogram
        f2 = figure()
        a2 = f2.add_subplot(111)
        a2.hold(True)
        a2.bar(left=Jijbins, height=nJij, width=Jijbins[1]-Jijbins[0], color='m', edgecolor='m')
        gcfm().frame.SetTitle(lastcmd())
        a2.set_title('Jij histogram\n%s' % lastcmd())
        a2.set_ylabel('probability density')
        a2.set_xlabel('Jij')
        a2.set_xlim(Jijrange)


class NetstateOtherStuff(BaseNetstate):
    """TEMPORARY"""

    def nspikingPMF(self, nis=None, shufflecodes=False, **kwargs):
        """Returns the PMF of observing n cells spiking in the same time bin, for either
        an unshuffled or shuffled Codes object"""
        observedwords = self.intcodes(nis=nis, shufflecodes=shufflecodes, **kwargs)
        # collect observances of the number of cells spiking for each pop code time bin
        nspiking = [ np.binary_repr(observedword).count('1') for observedword in observedwords ] # for all time bins, convert words to binary, count the number of 1s in each. np.binary_repr() is a bit faster than using neuropy.Core.bin()
        pnspiking, bins = histogram(nspiking, bins=arange(len(self.neurons)+1), normed='pmf') # histogram 'em, want all probs to add to 1, not their area, so use pmf
        return pnspiking, bins

    def plotnspikingPMFs(self, nis=None, xrange=[-0.5, 15], **kwargs):
        """Plots nspikingPMF, for both observed and shuffled (forcing independence) codes.
        See 2006 Schneidman fig 1e"""
        if nis == None:
            niswasNone = True
            nis = random.sample(self.co.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
            nis.sort()
            print 'nis = %r' % nis
        else:
            niswasNone = False
        observedpnspiking, observedbins = self.nspikingPMF(nis=nis, shufflecodes=False, **kwargs)
        indeppnspiking, indepbins = self.nspikingPMF(nis=nis, shufflecodes=True, **kwargs)
        assert (observedbins == indepbins).all() # paranoid schizo, just checking
        assert approx(observedpnspiking.sum(), 1.0), 'total observed probs: %f' % observedpnspiking.sum()
        assert approx(indeppnspiking.sum(), 1.0), 'total indep probs: %f' % indeppnspiking.sum()
        f = figure()
        a = f.add_subplot(111)
        a.hold(True)
        a.plot(observedbins, observedpnspiking, 'r.-')
        a.plot(indepbins, indeppnspiking, 'b.-')
        titlestr = 'PMF of observing n cells spiking in the same time bin'
        titlestr += '\n%s' % lastcmd()
        if niswasNone:
            titlestr += ', nis: %r' % nis
        a.set_title(titlestr)
        a.legend(('observed', 'indep (shuffled)'))
        a.set_yscale('log')
        a.set_xlim(xrange)
        gcfm().frame.SetTitle(lastcmd())
        a.set_xlabel('number of spiking cells in a bin')
        a.set_ylabel('probability')

    # if python ever gets class decorators, an inner class could be specified as:
    #@innerclass
    class Scatter(object):
        def __init__(self, nis=None, nbits=None, model='indep', randomneurons=True, shufflecodes=False,
                     algorithm='CG', **kwargs):
            """Netstate scatter analysis object. See Schneidman Figures 1f and 2a.
            Calculates the expected probabilities, assuming a model in ['indep', 'ising'],
            of all possible population codes vs their observed probabilities.
            nis are in LSB to MSB order."""
            if nis == None:
                nis = self.__outer__.co.nis # grab nis attrib from Codes object in outer Netstate object
            else:
                randomneurons = False # specific nis have been passed, don't use random neurons
                if nbits == None:
                    nbits = len(nis) # if nis is specified and nbits isn't, each ni gets its own bit
            if nbits == None: # nis and nbits were both passed as None
                nbits = CODEWORDLEN
            nbits = min(len(nis), nbits) # constrain nbits to be no more than the number of nis
            if randomneurons:
                nis = random.sample(nis, nbits) # randomly sample nbits of the nis
                nis.sort() # to keep things organized
            else:
                nis = nis[:nbits] # use just the first nbits neurons to make your words

            self.nis = nis # save it so that plot() method can access it
            self.nbits = nbits
            self.model = model

            #if randomneurons:
            #    print 'neurons:', nis # print 'em out if they were randomly selected
            self.scatterintcodes = self.intcodes(nis=nis, shufflecodes=shufflecodes, **kwargs)
            self.pobserved, self.observedwords = histogram(self.scatterintcodes, bins=arange(2**nbits), normed='pmf')
            #self.pobserved, self.observedwords = self.intcodesPDF(nis=nis, shufflecodes=shufflecodes, **kwargs) # potentially shuffle the observed codes
            if model == 'indep':
                self.pexpected, self.expectedwords = self.intcodesFPDF(nis=nis, **kwargs) # expected, assuming independence. don't potentially shuffle expected codes
            elif model == 'ising':
                ising = self.ising(nis=nis, algorithm=algorithm) # returns a maxent Ising model
                self.pexpected = ising.p # expected, assuming maxent Ising model
                self.expectedwords = ising.intsamplespace
            elif model == 'both':
                ising = self.ising(nis=nis, algorithm=algorithm) # returns a maxent Ising model
                self.pexpected = ising.p # expected, assuming maxent Ising model
                self.expectedwords = ising.intsamplespace
                self.pindepexpected, blah = self.intcodesFPDF(nis=nis, **kwargs) # expected, assuming independence. don't potentially shuffle expected codes
            else:
                raise ValueError, 'Unknown model %r' % model
            assert (self.observedwords == self.expectedwords).all() # make sure we're comparing apples to apples

        def plot(self, color=False):
            """Scatterplots the expected probabilities of all possible population codes (y axis)
            vs their observed probabilities (x axis).
            nis are in LSB to MSB order. See Schneidman Figures 1f and 2a"""
            f = figure()
            a = f.add_subplot(111)
            a.plot([10**-6, 1], [10**-6, 1], 'b-') # plot a y=x line
            a.hold(True)

            self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100)) # create a long tooltip with newline to get around bug where newlines aren't recognized on subsequent self.tooltip.SetTip() calls
            self.tooltip.Enable(False) # leave disabled for now
            self.tooltip.SetDelay(0) # set popup delay in ms
            gcfm().canvas.SetToolTip(self.tooltip) # connect the tooltip to the canvas
            f.canvas.mpl_connect('motion_notify_event', self._onmotion) # connect the mpl event to the action

            # pylab.scatter(pobserved, pexpected), followed by setting the x and y axes to log scale freezes the figure and runs 100% cpu
            # gca().set_xscale('log')
            # gca().set_yscale('log')
            # use loglog() instead

            # colour each scatter point according to how many 1s are in the population code word it represents.
            # This is done very nastily, could use a cleanup:
            if color:
                inds = []
                for nspikes in range(0, 5):
                    inds.append([])
                    [ inds[nspikes].append(i) for i in range(0, 2**self.nbits) if bin(i).count('1') == nspikes ]
                pobserved = self.pobserved.copy() # make local copies that are safe to modify for colour plotting and shit
                pexpected = self.pexpected.copy()
                pobserved1 = pobserved[inds[1]]; pexpected1 = pexpected[inds[1]]
                pobserved2 = pobserved[inds[2]]; pexpected2 = pexpected[inds[2]]
                pobserved3 = pobserved[inds[3]]; pexpected3 = pexpected[inds[3]]
                pobserved4 = pobserved[inds[4]]; pexpected4 = pexpected[inds[4]]
                pobserved[inds[1]], pexpected[inds[1]] = None, None # remove all these
                pobserved[inds[2]], pexpected[inds[2]] = None, None
                pobserved[inds[3]], pexpected[inds[3]] = None, None
                pobserved[inds[4]], pexpected[inds[4]] = None, None

            if self.model == 'both': # plot the indep model too, and plot it first
                a.loglog(self.pobserved, self.pindepexpected, color='lightsteelblue', marker='.', linestyle='None')
            # plot whichever model was specified
            if color:
                a.loglog(pobserved, pexpected, 'k.') # plots what's left in black
                a.loglog(pobserved4, pexpected4, 'm.')
                a.loglog(pobserved3, pexpected3, 'c.')
                a.loglog(pobserved2, pexpected2, 'y.')
                a.loglog(pobserved1, pexpected1, 'r.')
            else:
                a.loglog(self.pobserved, self.pexpected, 'r.')
            '''
            a.plot(pobserved, pexpected, 'k.')
            '''
            gcfm().frame.SetTitle(lastcmd())
            missingcodeis = (self.pobserved == 0).nonzero()[0]
            nmissing = len(missingcodeis)
            percentmissing = nmissing / float(2**self.nbits) * 100
            missingcodetext = ''
            if nmissing != 0:
                missingcodes = self.observedwords[missingcodeis]
                pexpectedmissing = self.pexpected[missingcodeis]
                maxpi = pexpectedmissing.argmax()
                maxp = pexpectedmissing[maxpi]
                maxpcode = self.expectedwords[missingcodeis[maxpi]]
                missingcodetext += '\n nmissingcodes: %d, maxpmissingcode: (%r, pexpected=%.3g)' % (nmissing, bin(maxpcode, minbits=self.nbits), maxp)
            a.set_title('%s\nneurons: %s' % (lastcmd(), self.nis))# + missingcodetext)
            a.set_xlabel('observed population code probability')
            a.set_ylabel('expected population code probability')
            a.set_ylim(ymin=10**-11, ymax=10**0) # this makes all plots consistent, some might be missing a scatter point or two
            if self.model =='both':
                DJSstring = '(%.4f, %.4f)' % (DJS(self.pobserved, self.pindepexpected), DJS(self.pobserved, self.pexpected))
            else:
                DJSstring = '%.4f' % DJS(self.pobserved, self.pexpected)
            a.text(0.99, 0.01, '%.1f%% missing\nDJS=%s' % (percentmissing, DJSstring), # add DJS to bottom right of plot
                transform = a.transAxes,
                horizontalalignment = 'right',
                verticalalignment = 'bottom')
            print 'nis = %r' % self.nis # print out the nis so they can easily be copied and pasted elsewhere

        def _onmotion(self, event):
            """Called during mouse motion over scatterplot figure. Pops up the corresponding
            population code word and its int representation when hovering over a neuron scatter point"""
            if event.xdata != None and event.ydata != None: # if mouse is inside the axes
                i  = approx(event.xdata, self.pobserved, rtol=5e-2, atol=0).nonzero()[0] # find for what indices (if any) xdata == pobserved
                ii = approx(event.ydata, self.pexpected[i], rtol=1e-1, atol=0).nonzero()[0] # for those above, find for what index (if any) ydata == pexpected
                codeis = i[ii]
                if codeis.size > 0:
                    #tip += 'i: %r' % i
                    #tip += '\nii: %r' % ii
                    #tip += '\ncodeis: %r' % codeis
                    intcodes = self.observedwords[codeis] # get the int rep for those indices from either self.observedwords[i] or self.expectedwords[i], doesn't matter which since they should be identical
                    codes = [ bin(intcode, minbits=self.nbits) for intcode in intcodes ]
                    tip =  'codes: %s' % repr(codes).replace('\'', '')
                    tip += '\nintcodes: %r' % list(intcodes)
                    activenis = [ list(asarray(self.nis)[::-1][charfind(code, '1')]) for code in codes ]
                    tip += '\nactivenis: %r' % activenis
                    tip += '\npattern counts: %r' % [ (self.scatterintcodes == intcode).sum() for intcode in intcodes ]
                    tip += '\npattern rates (Hz): %s' % repr([ '%.3g' % (p / self.__outer__.tres * 1e6) for p in self.pobserved[codeis] ]).replace('\'', '')
                    tip += '\n(pobserved, pexpected): %s' % repr(zip([ '%.3g' % val for val in self.pobserved[codeis] ], [ '%.3g' % val for val in self.pexpected[codeis] ])).replace('\'', '')
                    tip += '\npobserved / pexpected: %s' % repr([ '%.3g' % (o/float(e)) for o, e in zip(self.pobserved[codeis], self.pexpected[codeis]) ]).replace('\'', '')
                    tip += '\npexpected / pobserved: %s' % repr([ '%.3g' % (e/float(o)) for o, e in zip(self.pobserved[codeis], self.pexpected[codeis]) ]).replace('\'', '')
                    self.tooltip.SetTip(tip) # update the tooltip
                    self.tooltip.Enable(True) # make sure it's enabled
                else:
                    self.tooltip.Enable(False) # disable the tooltip
            else: # mouse is outside the axes
                self.tooltip.Enable(False) # disable the tooltip

    # overwrite Scatter class def'n as an inner class of the current outer class Netstate, can refer to __outer__ attrib
    Scatter = innerclass(Scatter)

    def I2vsIN(self, N=10, ngroups=15, tres=CODETRES):
        """Does Schneidman fig 2c. I2/IN vs IN, for ngroups of cells.
        This shows you what fraction of network correlation is accounted
        for by the maxent pairwise model.
        Plotting Icond-indep is different from S1 (see below),
        and sounds annoying and not worth it (see methods)"""
        niss = nCrsamples(objects=self.neurons.keys(),
                          r=N, # pick N neurons at random
                          nsamples=ngroups) # do it ngroups times
        I2s = []
        INs = []
        pd = wx.ProgressDialog(title='I2vsIN() progress', message='', maximum=ngroups, # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for groupi, nis in enumerate(niss):
            p1 = asarray(self.intcodesFPDF(nis=nis, tres=tres)[0]) # indep model
            p2 = self.ising(nis=nis).p # expected, assuming maxent Ising model
            pN = asarray(self.intcodesPDF(nis=nis, tres=tres)[0]) # observed word probs
            S1 = entropy_no_sing(p1) # ignore any singularities
            S2 = entropy_no_sing(p2)
            SN = entropy_no_sing(pN)
            IN = S1 - SN
            I2 = S1 - S2
            I2s.append(I2 / tres * 1e6) # convert to bits/sec
            INs.append(IN / tres * 1e6)
            cancel = not pd.Update(groupi, newmsg='groupi = %d' % (groupi+1))
            if cancel:
                pd.Destroy()
                return
        pd.Destroy()
        I2s = asarray(I2s)
        INs = asarray(INs)
        I2divIN = I2s / INs
        f = figure()
        gcfm().frame.SetTitle(lastcmd())
        a = f.add_subplot(111)
        a.plot(INs, I2divIN, 'r.')
        a.set_xlim(xmin=0.0)
        a.set_ylim(ymin=0.0, ymax=1.0)
        a.set_xlabel('IN (bits / sec)')
        a.set_ylabel('I2 / IN')
        a.set_title('%s' % lastcmd())
        a.text(0.99, 0.01, 'mean=%.3f, std=%.3f' % (I2divIN.mean(), I2divIN.std()), # add mean and std to bottom right
            transform = a.transAxes,
            horizontalalignment = 'right',
            verticalalignment = 'bottom')

        return dictattr(I2divIN=I2divIN, INs=INs, niss=niss, a=a)

    def DJShist(self, nbits=CODEWORDLEN, ngroups=5, models=['ising', 'indep'],
                logrange=(-3.667, -0.333), nbins=50, shufflecodes=False, algorithm='CG', **kwargs):
        """Plots Jensen-Shannon divergence histograms and a histogram of their ratios
        for ngroups random groups of cells, each of length nbits.
        See Schneidman figure 2b"""
        niss = [] # list of lists, each sublist is a group of neuron indices
        DJSs = {} # hold the Jensen-Shannon divergences for different models and different groups of neurons
        for model in models:
            DJSs[model] = [] # init a dict with the model names as keys, and empty lists as values
        pd = wx.ProgressDialog(title='DJShist progress', message='', maximum=ngroups*len(models), # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for groupi in range(ngroups): # for each group of nbits cells
            nis = random.sample(self.co.nis, nbits) # randomly sample nbits of the Netstate's nis attrib
            niss.append(nis)
            for modeli, model in enumerate(models): # for each model, use the same nis
                so = self.Scatter(nis=nis, model=model, randomneurons=False,
                                  shufflecodes=shufflecodes, algorithm=algorithm, **kwargs)
                DJSs[model].append(DJS(so.pobserved, so.pexpected))
                cancel = not pd.Update(groupi*len(models)+modeli, newmsg='groupi = %d\nmodel = %s' % (groupi, model))
                if cancel:
                    pd.Destroy()
                    return
        pd.Destroy()

        # now find the DJSratios between the two models, for each group of neurons
        if len(models) == 2: # do it only if there's 2 models, otherwise it's indeterminate which two to take ratio of
            DJSratios = asarray(DJSs.values()[1]) / asarray(DJSs.values()[0]) # 2nd model as a ratio of the 1st

        # histogram DJSs and DJSratios in logspace
        x = np.logspace(start=logrange[0], stop=logrange[1], num=nbins, endpoint=True, base=10.0)
        n = {} # stores a list of the bin heights in a separate key for each model
        for model in models:
            n[model] = histogram(DJSs[model], bins=x, normed=False)[0]
        nratios = histogram(DJSratios, bins=x, normed=False)[0] # bin heights for the DJSratios
        color = {'indep': 'b', 'ising': 'r'} # dict that maps from model name to color
        barwidths = list(diff(x)) # each bar will have a different width, convert to list so you can append
        # need to add one more entry to barwidth to the end to get nbins of them:
        barwidths.append(0) # don't display the last one
        logbinwidth = (logrange[1]-logrange[0]) / float(nbins)
        #barwidths.append(10**(logrange[1]+logbinwidth) - x[-1]) # should be exactly correct

        # plot DJSs of all models on the same axes
        f = figure()
        a = f.add_subplot(111)
        #a.hold(True)
        bars = {}
        heights = {}
        for model in models:
            heights[model] = n[model] / float(ngroups * logbinwidth)
            bars[model] = a.bar(left=x, height=heights[model], width=barwidths, color=color[model],
                                 edgecolor=color[model])
        a.set_xscale('log', basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
        a.set_xlim(xmin=10**logrange[0], xmax=10**logrange[1])
        gcfm().frame.SetTitle(lastcmd())
        #a.set_title('Jensen-Shannon divergence histogram\n%s' % lastcmd())
        #a.set_ylabel('number of groups of %d cells' % nbits)
        #a.set_xlabel('DJS (bits)')
        #a.set_ylabel('probability density (1 / log10(DJS))', size=20)

        a.set_xticklabels(['', '0.001', '0.01', '0.1', '']) # hack!

        for label in a.get_xticklabels():
            label.set_size(30)
        for label in a.get_yticklabels():
            label.set_size(30)

        a.legend([ bars[model][0] for model in models ], ['pairwise', 'independent'], prop=mpl.font_manager.FontProperties(size=20) ) # grab the first bar for each model, label it with the model name

        # plot DJSratios
        if len(models) == 2:
            f2 = figure()
            a2 = f2.add_subplot(111)
            a1 = a
            a = [a1, a2]
            a2.bar(left=x, height=nratios, width=barwidths, color='g', edgecolor='g')
            a2.set_xscale('log', basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
            gcfm().frame.SetTitle(lastcmd())
            a2.set_title('Jensen-Shannon divergence ratios histogram\n%s' % lastcmd())
            a2.set_ylabel('number of groups of %d cells' % nbits)
            a2.set_xlabel('DJS ratio (%s / %s)' % (models[1], models[0]))

        return dictattr(niss=niss, DJSs=DJSs, DJSratios=DJSratios, a=a, x=x, heights=heights)

    def S1INvsN(self, minN=4, maxN=15, maxnsamples=10, tres=CODETRES):
        """Plots the average independent cell entropy S1 and average network multi-information IN (IN = S1 - SN)
        vs network size N. IN is how much less entropy there is in the system due to correlated network activity.
        For each network size up to maxN, averages S1 and IN over maxnsamples (or less if that many aren't possible)
        number of groups at each value of N"""
        S1ss= [] # as f'n of N
        INss = []
        N = range(minN, maxN+1) # network sizes from minN up to maxN
        #tstart = time.clock()
        nsamples = [ min(nCr(self.nneurons, r), maxnsamples) for r in N ] # nsamples as a f'n of N. For each value of N, take up to maxnsamples of all the other neurons, if that many are even possible

        pd = wx.ProgressDialog(title='S1INvsN progress', message='', maximum=sum(nsamples), # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for ni, n in enumerate(N): # for all network sizes
            # get a list of lists of neuron indices
            niss = nCrsamples(objects=self.neurons.keys(),
                              r=n, # pick n neurons
                              nsamples=nsamples[ni] ) # do it at most maxnsamples times
            S1s = []
            INs = []
            for nisi, nis in enumerate(niss):
                cancel = not pd.Update(ni*maxnsamples+nisi, newmsg='N = %d\nsamplei = %d' % (n, nisi))
                if cancel:
                    pd.Destroy()
                    return
                #t2 = time.clock()
                p1 = asarray(self.intcodesFPDF(nis=nis, tres=tres)[0]) # indep model
                pN = asarray(self.intcodesPDF(nis=nis, tres=tres)[0]) # observed word probs
                #print 'calcing ps took: %f sec' % (time.clock()-t2)
                S1 = entropy_no_sing(p1) # ignore any singularities
                SN = entropy_no_sing(pN)
                assert S1 > SN or approx(S1, SN), 'S1 is %.20f, SN is %.20f' % (S1, SN) # better be, indep model assumes the least structure
                IN = S1 - SN
                #print S1, SN, IN
                S1s.append(S1 / tres * 1e6) # convert to bits/sec
                INs.append(IN / tres * 1e6)
            S1ss.append(S1s)
            INss.append(INs)
        pd.Destroy()
        S1mean = [ asarray(S1s).mean() for S1s in S1ss ]
        S1std = [ asarray(S1s).std() for S1s in S1ss ]
        S1sem = asarray(S1std) / sqrt(asarray(nsamples))
        INmean = [ asarray(INs).mean() for INs in INss ]
        INstd = [ asarray(INs).std() for INs in INss ]
        INsem = asarray(INstd) / sqrt(asarray(nsamples))
        f = figure()
        gcfm().frame.SetTitle(lastcmd())
        a = f.add_subplot(111)
        a.hold(True)
        for n, S1s in zip(N, S1ss): # plot all the samples before plotting the means with errorbars
            a.plot([n]*len(S1s), S1s, '_', markersize=4, color='lightblue')
        for n, INs in zip(N, INss):
            a.plot([n]*len(INs), INs, '_', markersize=4, color='pink')
        S1line = a.errorbar(N, S1mean, yerr=S1sem, fmt='b.')[0]
        INline = a.errorbar(N, INmean, yerr=INsem, fmt='r.')[0]
        # do some linear regression in log10 space
        mS1, bS1 = sp.polyfit(log10(N), log10(S1mean), 1) # returns slope and y intercept
        mIN, bIN = sp.polyfit(log10(N), log10(INmean), 1)
        xintersect = (bIN - bS1) / (mS1 - mIN)
        x = array([-1, 3]) # define x in log10 space, this is really [0.1, 1000]
        yS1 = mS1*x + bS1 # y = mx + b
        yIN = mIN*x + bIN
        a.plot(10.0**x, 10.0**yS1, 'b-') # raise them to the power to make up for the fact that both the x and y scales will be log
        a.plot(10.0**x, 10.0**yIN, 'r-')
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim(1e0, 1e3)
        a.set_ylim(1e-2, 1e4)
        a.set_xlabel('Number of cells')
        a.set_ylabel('bits / sec')
        a.set_title('S1 & IN vs N\n%s' % lastcmd())
        a.legend((S1line, INline), ('S1, slope=%.3f' % mS1, 'IN, slope=%.3f' % mIN), loc='lower right')
        a.text(0.99, 0.98, 'Nc=%d' % np.round(10**xintersect), # add text box to upper right corner of axes
            transform = a.transAxes,
            horizontalalignment = 'right',
            verticalalignment = 'top')
        return S1mean, INmean, N

    def NMmutualinfo(self, nis=None, mis=None, Nbinarray=None, Mbinarray=None, verbose=False):
        """Calculates information that N cells provide about M cells (ie,
        their mutual information), as a fraction of the M cells' marginal entropy.
        nis is neuron indices of length N, and is ordered LSB to MSB.
        mis is neuron indices of length M"""

        if nis != None:
            nis = toiter(nis)
            N = len(nis)
            Nintcodes = self.intcodes(nis=nis)
            #print 'first 100 Nintcodes\n', Nintcodes[:100].__repr__()
        elif Nbinarray != None:
            Nbinarray = to2d(Nbinarray) # make it 2D if it's 1D
            N = len(Nbinarray) # gets the number of rows
            Nintcodes = binarray2int(Nbinarray)
        else:
            raise ValueError, 'nis and Nbinarray args can''t both be None'

        if mis != None:
            mis = toiter(mis)
            M = len(mis)
            Mintcodes = self.intcodes(nis=mis)
            #print 'first 100 Mintcodes\n', Mintcodes[:100].__repr__()
        elif Mbinarray != None:
            Mbinarray = to2d(Mbinarray) # make it 2D if it's 1D
            M = len(Mbinarray) # gets the number of rows
            Mintcodes = binarray2int(Mbinarray)
        else:
            raise ValueError, 'mis and Mbinarray args can''t both be None'

        # build up joint pdf of all the possible N words, and the two possible N+1th values (0 and 1)
        xedges = arange(2**N+1) # values 0 to 2**N - 1, plus 2**N which is needed as the rightmost bin edge for histogram2d (annoying)
        yedges = arange(2**M+1)
        bins = [xedges, yedges]
        jpdf, xedgesout, yedgesout = histogram2d(Nintcodes, Mintcodes, bins, normed='pmf') # generate joint pdf
        #print 'jpdf\n', jpdf.__repr__()
        #print 'jpdf.sum()', jpdf.sum()
        assert (np.float64(xedges) == xedgesout).all()
        assert (np.float64(yedges) == yedgesout).all() # make sure we know what we're doing

        # pdf of N cells
        #Npdf, Nedges = histogram(Nintcodes, bins=range(2**N), normed='pmf')
        #print 'first 100 Npdf\n', Npdf[:100].__repr__()

        # pdf of M cells
        #Mpdf, Medges = histogram(Mintcodes, bins=arange(2**M), normed='pmf')
        #print 'first 100 Mpdf\n', Mpdf[:100].__repr__()

        marginalMpdf = jpdf.sum(axis=0)
        #assert approx(Mpdf, marginalMpdf).all() # make sure what you get from the joint is what you get when just building up the pdf straight up on its own

        I = mutualinfo(jpdf)

        IdivS = I / entropy(marginalMpdf) # return mutual info as fraction of entropy in M group of cells

        if verbose:
            print 'nis', nis
            print 'mis', mis
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
        return IdivS

    def NNplus1(self, Nplus1s=None, maxN=15, maxnsamples=10):
        """Does Schneidman Figure 5b. Averages over as many as maxnsamples different
        groups of N cells for each N+1th cell in Nplus1s,
        all done for different values of N up to maxN"""
        if Nplus1s == None: # list of all indices of neurons that will be treated as the N+1th neuron
            Nplus1s = self.co.nis
        else:
            Nplus1s = toiter(Nplus1s)
        nNplus1s = len(Nplus1s)
        dims = (maxN, nNplus1s, maxnsamples)
        mask = np.zeros(dims) # this will be converted to an array of Falses
        IdivS = np.ma.array(mask, mask=mask, fill_value=666) # masked array that holds the mutual info between N and N+1th cells, as a ratio of the N+1th cell's entropy. Index like: IdivS[ni, Nplus1i, samplei], ie group size, N+1th cell you're comparing to, and number of samples of size N taken from the possible combs
        N = range(1, maxN+1) # cell group size, excluding the N+1th neuron. This will be the x axis in the plot
        #N=[15]#.reverse() # for fun and pleasure
        nsamples = [ min(maxnsamples, nCr(nNplus1s-1, r)) for r in N ] # take up to maxnsamples of all the other neurons, if that many even exist (for the lower N values, the constraint may end up being the total number of possible combinations of cells), for each N+1th cell. Taken from nNplus1s-1 cuz you always have to exclude an N+1th neurons
        for ni, n in enumerate(N): # for all group sizes
            IdivS.mask[ni, :, nsamples[ni]::] = True # mask out the sampleis that are out of range for this value of N, if any
        maximum = nNplus1s*sum(nsamples)
        pd = wx.ProgressDialog(title='NNplus1 progress', message='', maximum=maximum, # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)

        # get the binary array for the whole population, then index into it appropriately in the sample loop, find the corresponding integer codes, and feed it to NMmutualinfo, so you don't have to unnecessarily re-generate it on every iteration
        nis = self.co.nis
        nis2niis = self.co.nis2niis
        counter = 0 # counts inner loop for the progress dialog
        for ni, n in enumerate(N):
            for Nplus1i, Nplus1 in enumerate(Nplus1s): # for each N+1th neuron to compare to
                mii = nis2niis(Nplus1)
                niscopy = copy(nis) # make a copy of neuron indices
                niscopy.remove(Nplus1) # keep just the indices of all the other neurons
                samples = nCrsamples(niscopy, n, nsamples[ni]) # returns nsamples random unique choices of n items from niscopy
                for samplei, sample in enumerate(samples): # collect nsamples different combinations of the N other cells
                    niis = np.array([ nis2niis(s) for s in toiter(sample) ]) # most of the time (for n>1), sample will be a sequence of nis. Build an array of niis out of it to use as indices into the binary code array. Sometimes (for n=1) sample will be a scalar, hence the need to push it through toiter()
                    IdivS[ni, Nplus1i, samplei] = self.NMmutualinfo(Nbinarray=self.co.c[niis], Mbinarray=self.co.c[mii]) # do it
                    cancel = not pd.Update(counter, newmsg='N: %d; N+1th neuron: %d; samplei: %d' % (n, Nplus1, samplei))
                    if cancel:
                        pd.Destroy()
                        return
                    counter += 1
        pd.Destroy()
        IdivS = IdivS.reshape(maxN, nNplus1s*maxnsamples) # reshape such that you collapse all Nplus1s and samples into a single dimension (columns). The N are still in the rows
        #logIdivS = log10(IdivS)
        IdivSmeans = IdivS.mean(axis=1) # average over all Nplus1s and all samples. Values that are masked are ignored
        IdivSstds = IdivS.std(axis=1) # find stdev for the same
        assert IdivSmeans.shape == (maxN,)
        assert IdivSstds.shape == (maxN,)
        IdivSsems = IdivSstds / sqrt(asarray(nsamples)*nNplus1s)

        # plot the figure with error bars
        f = figure()
        gcfm().frame.SetTitle(lastcmd())
        a = f.add_subplot(111)
        a.hold(True)
        for n, row in zip(N, IdivS): # underplot the samples for each value of N
            a.plot([n]*len(row), row, '_', markersize=4, color='deepskyblue')
        a.errorbar(N, IdivSmeans, yerr=IdivSsems, fmt='b.') # now plot the means and sems
        # do some linear regression in log10 space
        m, b = sp.polyfit(log10(N), log10(IdivSmeans), 1) # returns slope and y intercept
        x = array([log10(0.9), 3]) # define x in log10 space, this is really [0.9, 1000]
        y = m*x + b
        xintersect = (0-b) / m # intersection point of regression line with y=1=10**0 line
        plot(10.0**x, 10.0**y, 'b-') # raise them to the power to make up for the fact that both the x and y scales will be log
        plot(10.0**x, [1e0]*2, 'r--') # plot horizontal line at y=1
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim(10**log10(0.9), 1e3)
        a.set_ylim(1e-3, 1e1)
        a.set_xlabel('Number of cells')
        a.set_ylabel('mutualinfo(N, N+1th) / entropy(N+1th)')
        a.set_title('fraction of info that N cells provide about the N+1th cell\n%s' % lastcmd())
        a.text(0.99, 0.98, 'Nc=%d' % np.round(10**xintersect), # add text box to upper right corner of axes
            transform = a.transAxes,
            horizontalalignment = 'right',
            verticalalignment = 'top')
        a.text(0.99, 0.01, 'slope=%.3f' % m, # add slope of fit line to bottom right
            transform = a.transAxes,
            horizontalalignment = 'right',
            verticalalignment = 'bottom')
        '''
        # plot the distributions of IdivS
        for ni, n in enumerate(N):
            f = figure()
            gcfm().frame.SetTitle('%s IdivS distrib for N=%d' % (lastcmd(), n))

            notmaskedis = IdivS[ni].mask==False # indexes the non-masked entries in IdivS, for this ni

            a1 = f.add_subplot(211) # axes with linear bins
            heights, bins = histogram(IdivS[ni, notmaskedis], bins=arange(0, 1, 0.02))
            barwidth = bins[1]-bins[0]
            a1.bar(left=bins, height=heights, width=barwidth, bottom=0, color='k')
            #a1.set_xlabel('mutualinfo(N, N+1th) / entropy(N+1th)')
            a1.set_ylabel('count')
            a1.set_title('IdivS distrib for N=%d' % n)

            a2 = f.add_subplot(212) # axes with log bins
            start = log10(0.001)
            stop = log10(1)
            bins = np.logspace(start=start, stop=stop, num=50, endpoint=True, base=10.0)
            heights, bins = histogram(IdivS[ni, notmaskedis], bins=bins)
            barwidth = list(diff(bins)) # each bar will have a different width, convert to list so you can append
            # need to add one more entry to barwidth to the end to get nbins of them:
            #barwidth.append(barwidth[-1]) # not exactly correct
            logbinwidth = (log10(stop)-log10(stop)) / float(len(bins))
            barwidth.append(10**(log10(stop)+logbinwidth) - stop) # should be exactly correct
            a2.bar(left=bins, height=heights, width=barwidth, bottom=0, color='k')
            a2.set_xscale('log')
            a2.set_xlabel('mutualinfo(N, N+1th) / entropy(N+1th)')
            a2.set_ylabel('count')
        return
        '''
    def _checkcell(self, ni=None, othernis=None, shufflecodes=False):
        """Returns the joint pdf of cell ni activity and the number of cells in
        othernis being active at the same time. ni should not be in othernis"""
        assert ni not in othernis
        nis2niis = self.co.nis2niis
        nii = nis2niis(ni)
        otherniis = nis2niis(othernis)
        nothers = len(othernis)

        nicode = self.co.c[nii] # 0s and 1s, this picks out the row in the binary code array that corresponds to ni
        othercodes = self.co.c[otherniis]
        if shufflecodes:
            nicode = asarray(shuffle(nicode))
            othercodes = asarray(shuffle(othercodes))
        nothersactive = othercodes.sum(axis=0) # anywhere from 0 up to and including nothers

        # build up joint pdf of the nicode and nothersactive
        xedges = np.array([0, 1, 2]) # values 0 and 1, plus 2 which is needed as the rightmost bin edge for histogram2d (annoying)
        yedges = arange(nothers+2) # anywhere from 0 up to and including nothers, plus nothers+1 as the rightmost bin edge
        bins = [xedges, yedges]

        jpdf, xedgesout, yedgesout = histogram2d(nicode, nothersactive, bins, normed=False) # generate joint pdf, nicode are in the rows, nothersactive are in the columns, leave it unnormalized, just counts

        # now, normalize each column separately, so that say, for nothersactive==5, p(checkcell==0)+p(checkcell==1) == 1.0
        jpdf = np.float64(jpdf) # convert to floats, updated entries are trunc'd to ints
        for coli in range(jpdf.shape[-1]):
            jpdf[:, coli] = normalize(jpdf[:, coli]) # save the normalized column back to the jpdf
        return jpdf

    def checkcells(self, nis=None, othernis=None, nothers=None, nsamples=10, shufflecodes=False):
        """Plots the probability of each cell (in nis) being active vs. the number of
        other active cells (in the Recording) at that time. For each ni, an average over
        nsamples, each being a different sample of nothers from othernis.
        See Schneidman figure 5c"""
        if nis == None:
            nis = self.co.nis
        else:
            nis = toiter(nis)
        if othernis == None:
            othernis = self.co.nis
        else:
            othernis = toiter(nis)
        if nothers == None:
            less = 0
            while nCr(len(othernis), len(othernis)-less) < nsamples:
                less += 1
            nothers = len(othernis) - 1 - less # -1 to remove ni, -less again to allow for at least nsamples combos of othernis
        N = arange(nothers+1)
        saved_othernis = copy(othernis) # save a copy so we can mess with the original

        for ni in nis:
            if saved_othernis == None:
                othernis = copy(self.co.nis)
            else:
                othernis = copy(saved_othernis)
            try:
                othernis.remove(ni) # all the other possible nis, excluding the current ni
            except ValueError: # ni isn't in othernis, nothing to remove
                pass
            otherniss = nCrsamples(objects=othernis, r=nothers, nsamples=nsamples) # get nsamples unique random samples of length nothers from othernis

            jpdfs = []
            for othernis in otherniss: # collect jpdfs across all random samples
                jpdf = self._checkcell(ni=ni, othernis=othernis, shufflecodes=shufflecodes)
                jpdfs.append(jpdf)

            jpdfs = asarray(jpdfs) # this is an nsamples x 2 x (1+nothers) matrix
            jpdfmean = jpdfs.mean(axis=0) # find the mean jpdf across all nsamples jpdfs
            jpdfstd = jpdfs.std(axis=0) # find the stdev across all nsamples jpdfs
            jpdfsem = jpdfstd / sqrt(nsamples)

            # plot it
            f = figure()
            gcfm().frame.SetTitle('%s for ni=%d' % (lastcmd(), ni))
            a = f.add_subplot(111)
            a.hold(True)
            # plot all the samples first
            for jpdf in jpdfs: # iter over the hyperrows
                a.plot(N, jpdf[1], '_', markersize=4, color=0.6) # marginal pdf of getting a 1 for the check cell
            # plot the stdevs, means, and sems of marginal pdf of getting a 1 for the check cell
            #a.errorbar(N, jpdfmean[1], yerr=jpdfstd[1], fmt=None, capsize=0, ecolor='grey')
            a.errorbar(N, jpdfmean[1], yerr=jpdfsem[1], fmt='k.-')
            a.set_ylim(ymin=0, ymax=1)

            titlestr = '%s\nni=%d' % (lastcmd(), ni)
            if nsamples == 1:
                titlestr += ', othernis=%r' % othernis
            else:
                titlestr += ', nsamples=%d' % nsamples
            a.set_title(titlestr)
            a.set_xlabel('Number of other active cells')
            a.set_ylabel('Probability of cell ni being active')


class RecordingNetstate(BaseRecording):
    """Mix-in class that defines Netstate related Recording methods"""
    def ns_(self, experiments=None, kind=CODEKIND, tres=CODETRES, phase=0):
        """Returns a BaseNetstate object"""
        return BaseNetstate(recording=self, experiments=experiments, kind=kind, tres=tres, phase=phase)
    def ns_isinghist(self, experiments=None, kind=CODEKIND, tres=CODETRES, phase=0):
        """Returns a NetstateIsingHist object"""
        return NetstateIsingHist(recording=self, experiments=experiments, kind=kind, tres=tres, phase=phase)
    def ns_other(self, experiments=None, kind=CODEKIND, tres=CODETRES, phase=0):
        """Returns a NetstateOtherStuff object"""
        return NetstateOtherStuff(recording=self, experiments=experiments, kind=kind, tres=tres, phase=phase)



class Recording(RecordingRaster,
                RecordingCode,
                RecordingNetstate,
                BaseRecording):
    """Inherits all the Recording classes into a single Recording class"""
    pass
