"""Defines the Recording class and all of its support classes"""

#print 'importing Recording'

import Core
from Core import *
from Core import _data # ensure it's imported, in spite of leading _

# Good global setting for presentation plots:
pl.rcParams['axes.labelsize'] = 30
pl.rcParams['xtick.labelsize'] = 25
pl.rcParams['ytick.labelsize'] = 25
pl.rcParams['xtick.major.size'] = 7
pl.rcParams['ytick.major.size'] = 7
pl.rcParams['lines.markersize'] = 10
# use gca().set_position([0.15, 0.15, 0.8, 0.8]) or just the 'configure subplots' widget to make all the labels fit within the figure


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
    def raster(self, experiments=None, nis=None,
                     jumpts=None, binwidth=None, relativet0=False, units='msec',
                     publication=False):
        """Returns a population spike raster plot"""
        return PopulationRaster(recording=self, experiments=experiments, nis=nis,
                                                jumpts=jumpts, binwidth=binwidth, relativet0=relativet0, units=units,
                                                publication=publication)
    raster.__doc__ += '\n\n'+PopulationRaster.__doc__


class Codes(object):
    """A 2D array where each row is a neuron code, and each column
    is a binary population word for that time bin, sorted LSB to MSB from top to bottom.
    neurons is a list of Neurons, also from LSB to MSB. Order in neurons is preserved."""
    def __init__(self, neurons=None, tranges=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE, shufflecodes=False):
        self.neurons = neurons
        self.tranges = tolist(tranges)
        self.kind = kind
        self.tres = tres
        self.phase = phase
        self.shufflecodes = shufflecodes
        self.nis = [ neuron.id for neuron in self.neurons ]
        self.nneurons = len(self.neurons)
        self.nis2niisdict = dict(zip(self.nis, range(self.nneurons))) # make a dict from keys:self.nis, vals:range(self.nneurons). This converts from nis to niis (from neuron indices to indices into the binary code array self.c)
    def nis2niis(self, nis=None):
        """Converts from nis to niis (from neuron indices to indices into the binary code array self.c).
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
            codeo = neuron.code(tranges=self.tranges, kind=self.kind, tres=self.tres, phase=self.phase)
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
    """A PDF of the correlations of the codes of all cell pairs (or of all cell pairs within
    some radius in um) in this Recording. See Schneidman2006 fig 1d"""
    def __init__(self, recording=None, experiments=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
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
        self.kind = kind
        self.tres = tres
        self.phase = phase
    '''
    # this was used to save on calc time by seeing if a CCPDF object with the same attribs had already been calc'd, seems dumb and unsafe, commented out
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
    '''
    def calc(self, radius=None, shuffleids=False):
        """Works on ConstrainedNeurons, but is constrained even further if experiments
        were passed and their tranges were used to generate self.tranges (see __init__)"""
        self.radius = radius
        self.shuffleids = shuffleids
        cnis = self.r.cn.keys() # ConstrainedNeuron indices
        ncneurons = len(cnis)
        # it's more efficient to precalculate the means and stds of each cell's codetrain,
        # and then reuse them in calculating the correlation coefficients:
        means = dict( ( cni, self.r.code(cni, tranges=self.tranges,
                                              kind=self.kind,
                                              tres=self.tres,
                                              phase=self.phase).c.mean() ) for cni in cnis ) # store each code mean in a dict
        stds  = dict( ( cni, self.r.code(cni, tranges=self.tranges,
                                              kind=self.kind,
                                              tres=self.tres,
                                              phase=self.phase).c.std() ) for cni in cnis ) # store each code std in a dict
        if self.shuffleids:
            scnis = shuffle(cnis) # shuffled neuron ids, this is a control to see if it's the radius of neurons included in the analysis, or the number of neurons included that's important. Turns out, neither is. The distribs look pretty much the same regardless of radius or shuffling
        else:
            scnis = cnis

        self.corrs = []
        for cnii1 in range(ncneurons):
            for cnii2 in range(cnii1+1, ncneurons):
                cni1 = cnis[cnii1]; scni1 = scnis[cnii1]
                cni2 = cnis[cnii2]; scni2 = scnis[cnii2]
                if radius == None or dist(self.r.cn[scni1].pos, self.r.cn[scni2].pos) < self.radius:
                    code1 = self.r.code(cni1, tranges=self.tranges, kind=self.kind, tres=self.tres, phase=self.phase).c
                    code2 = self.r.code(cni2, tranges=self.tranges, kind=self.kind, tres=self.tres, phase=self.phase).c
                    cc = ((code1 * code2).mean() - means[cni1] * means[cni2]) / (stds[cni1] * stds[cni2]) # (mean of product - product of means) / by product of stds
                    self.corrs.append(cc)
        self.corrs = array(self.corrs)
        '''
        # simpler, but slower way:
        self.corrs = [ self.r.codecorr(cnis[cnii1], cnis[cnii2], tranges=self.tranges, kind=self.kind, self.tres, self.phase)
                       for cnii1 in range(0,ncneurons) for cnii2 in range(cnii1+1,ncneurons) ]
        '''
    def plot(self, figsize=(7.5, 6.5), crange=[-0.1, 0.5], limitstats=True, nbins=30, normed='pdf'):
        """Plots the corrs. If limitstats, the stats displayed exclude any corr values that fall outside of crange"""
        self.crange = crange
        self.nbins = nbins
        self.normed = normed
        f = figure(figsize=figsize)
        a = f.add_subplot(111)
        try: # figure out the bin edges
            bins = np.linspace(start=self.crange[0], stop=self.crange[1], num=self.nbins, endpoint=True)
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
        self.mean = mean(corrs)
        self.median = median(corrs)
        argmode = n.argmax()
        self.mode = mean([c[argmode], c[argmode + 1]]) # find middle of tallest bin

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
        #titlestr = 'neuron pair code correlation pdf'
        #titlestr += '\n%s' % lastcmd()
        titlestr = '%s' % lastcmd()
        a.set_title(titlestr)
        '''
        if self.normed:
            if self.normed == 'pmf':
                a.set_ylabel('probability mass')
            else:
                a.set_ylabel('probability density')
        else:
            a.set_ylabel('count')
        a.set_xlabel('correlation coefficient')
        '''
        a.text(0.99, 0.99, 'mean = %.3f\nmedian = %.3f\nmode = %.3f'
                            % (self.mean, self.median, self.mode), # add stuff to top right of plot
                            transform = a.transAxes,
                            horizontalalignment='right',
                            verticalalignment='top')



class RecordingCode(BaseRecording):
    """Mix-in class that defines the spike code related Recording methods"""
    def code(self, cneuron=None, tranges=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a ConstrainedNeuron.Code object, constrained to the time
        ranges of the Experiments in this Recording, as well as by tranges. Takes either a
        ConstrainedNeuron object or just a ConstrainedNeuron id"""
        try:
            return cneuron.code(tranges=tranges, kind=kind, tres=tres, phase=phase) # see if cneuron is a ConstrainedNeuron
        except AttributeError:
            return self.cn[cneuron].code(tranges=tranges, kind=kind, tres=tres, phase=phase) # cneuron is probably a ConstrainedNeuron id

    def codes(self, neurons=None, experiments=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE, shufflecodes=False):
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
        codeso = Codes(neurons=neurons, tranges=tranges, kind=kind, tres=tres, phase=phase, shufflecodes=shufflecodes)
        codeso.calc()
        return codeso
    codes.__doc__ += '\n\nCodes object:\n' + Codes.__doc__

    def codecorr(self, neuron1, neuron2, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Calculates the correlation coefficient of the codes of two neurons"""
        code1 = self.code(neuron1, kind=kind, tres=tres, phase=phase)
        code2 = self.code(neuron2, kind=kind, tres=tres, phase=phase)
        return corrcoef(code1.c, code2.c)

    def codecorrpdf(self, experiments=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE, radius=None, shuffleids=False):
        """Returns an existing CodeCorrPDF object, or creates a new one if necessary"""
        # 2007-10-26 Got rid of all this saved calc() stuff, calculating these really doesn't take too long, and it's annoying having to redo it manually, and possibly unsafe
        #try:
        #    self._codecorrpdfs
        #except AttributeError: # doesn't exist yet
        #    self._codecorrpdfs = [] # create a list that'll hold CodeCorrPDF objects
        cco = CodeCorrPDF(recording=self, experiments=experiments, kind=kind, tres=tres, phase=phase) # init a new one
        #for ccpdf in self._codecorrpdfs:
        #    if cco == ccpdf: # need to define special == method for class CodeCorrPDF()
        #        return ccpdf # returns the first object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codecorrpdfs
        cco.calc(radius, shuffleids) # no matching object was found, calculate it
        #self._codecorrpdfs.append(cco) # add it to the object list
        return cco


class BaseNetstate(object):
    """Base class of Network state analyses.
    Implements a lot of the analyses on network states found in the 2006 Schneidman paper

    WARNING!!!!!! not sure if self.tranges, which derives from self.experiments, is being used at all yet!!!!!!!!!!!!!

    """
    def __init__(self, recording, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
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
        if nis == None:
            self.niswasNone = True
            nis = self.neurons.keys() # get all neuron indices in this Recording
            nis.sort() # make sure they're sorted
        else:
            self.niswasNone = False
        self.cs = self.codes(nis=nis) # generate and save the Codes object for all the nis
        # if you need to retrieve the nis, you can get them from self.cs.nis. Leave self.nis open for subclasses to use for their own purposes

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

    def get_wordts(self, nis=None, mis=None):
        """Returns word times, ie the times of the left bin edges for which all the
        neurons in the mis in this Netstate object have a 1 in them, and all
        the rest have a 0 in them. nis lists the total population of neuron ids"""
        if nis == None:
            nis = self.cs.nis
        mis = toiter(mis)
        for mi in mis:
            assert mi in nis # make sure mis is a subset of nis
        cs = self.codes(nis=nis) # make a new codes object using the nis population
        nis2niis = cs.nis2niis
        notmis = [ ni for ni in nis if ni not in mis ] # nis not in mis
        mis_high = cs.c[nis2niis(mis)].prod(axis=0) == 1 # take product down all rows, only synchronous events across all mis cells will survive, boolean array
        notmis_low = cs.c[nis2niis(notmis)].sum(axis=0) == 0 # boolean array
        i = (mis_high * notmis_low).nonzero()[0] # indices where mis are 1 and all the others are 0
        return cs.t[i] # return the times at those indices

    def get_wordtsms(self, nis=None, mis=None):
        """Returns word times to the nearest msec, with the on bits specified in mis.
        nis lists the total population of neuron ids"""
        return np.int32(np.round(self.get_wordts(nis=nis, mis=mis) / 1e3))

    def get_intcodes(self, nis=None, shufflecodes=False):
        """Given neuron indices (ordered LSB to MSB top to bottom), returns an array of the integer representation
        of the neuronal population binary code for each time bin"""
        assert self.kind == 'binary'
        if nis == None:
            nis = random.sample(self.cs.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
        return binarray2int(self.codes(nis=nis, shufflecodes=shufflecodes).c)

    def intcodesPDF(self, nis=None):
        """Returns the observed pdf across all possible population binary code words,
        labelled according to their integer representation"""
        if nis == None:
            nis = random.sample(self.cs.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
        intcodes = self.get_intcodes(nis=nis)
        nbits = len(nis)
        p, bins = histogram(intcodes, bins=arange(2**nbits), normed='pmf')
        return p, bins

    def intcodesFPDF(self, nis=None):
        """the F stands for factorial. Returns the probability of getting each population binary code word, assuming
        independence between neurons, taking into account each neuron's spike (and no spike) probability"""
        if nis == None:
            nis = random.sample(self.cs.nis, CODEWORDLEN) # randomly sample CODEWORDLEN bits of the nis
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

    def ising(self, nis=None, radius=None, algorithm='CG'):
        """Returns an Ising maximum entropy model that takes into account pairwise correlations within neuron codes
        algorithm can be 'CG', 'BFGS', 'LBFGSB', 'Powell', or 'Nelder-Mead'"""
        if nis == None:
            nis = self.cs.nis[0:CODEWORDLEN]
        #print 'nis:', nis.__repr__()
        codeso = self.codes(nis=nis)
        #c = codeso.c
        # convert values in codes object from [0, 1] to [-1, 1] by mutliplying by 2 and subtracting 1
        c = codeso.c.copy() # don't modify the original
        c = c*2 - 1 # this should be safe to do cuz c is a 2D array of signed int8 values
        #print 'c:', c.__repr__()
        means = [ row.mean() for row in c ] # iterate over rows of codes in c
        nrows = c.shape[0]
        pairmeans = []
        for i in range(0, nrows):
            for j in range(i+1, nrows):
                if radius == None or dist(self.r.n[nis[i]].pos, self.r.n[nis[j]].pos) < radius:
                    pairmeans.append((c[i]*c[j]).mean()) # take a pair of rows, find the mean of their elementwise product
                else:
                    pairmeans.append(None) # pair are outside the radius, ignore their pairmeans
        isingo = Core.Ising(means=means, pairmeans=pairmeans, algorithm=algorithm)
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

        pd = wx.ProgressDialog(title='NetstateIsingHist progress', message='', maximum=self.ngroups, # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for groupi in range(self.ngroups): # for each group of nbits cells
            nis = random.sample(self.cs.nis, nbits) # randomly sample nbits of the Netstate's nis attrib
            im = self.ising(nis=nis, algorithm=algorithm) # returns a maxent Ising model
            self.ims.append(im)
            self.his.append(im.hi)
            self.Jijs.append(im.Jij)
            cont, skip = pd.Update(groupi, newmsg='groupi = %d' % groupi)
            if not cont:
                pd.Destroy()
                return
        pd.Destroy()

        return self

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
        a1.set_title('hi histogram\n%s, nbits=%d, ngroups=%d, algorithm=%s' % (lastcmd(), self.nbits,
                                                                               self.ngroups, self.algorithm))
        a1.set_ylabel('probability density')
        a1.set_xlabel('hi')
        a1.set_xlim(hirange)

        # plot the Jij histogram
        f2 = figure()
        a2 = f2.add_subplot(111)
        a2.hold(True)
        a2.bar(left=Jijbins, height=nJij, width=Jijbins[1]-Jijbins[0], color='m', edgecolor='m')
        gcfm().frame.SetTitle(lastcmd())
        a2.set_title('Jij histogram\n%s, nbits=%d, ngroups=%d, algorithm=%s' % (lastcmd(), self.nbits,
                                                                                self.ngroups, self.algorithm))
        a2.set_ylabel('probability density')
        a2.set_xlabel('Jij')
        a2.set_xlim(Jijrange)

        self.f = {1:f1, 2:f2}
        self.a = {1:a1, 2:a2}
        return self


class NetstateNspikingPMF(BaseNetstate):
    """Netstate PMF of number of cells spiking in the same bin. See 2006 Schneidman fig 1e"""
    def calc(self, nbits=CODEWORDLEN):
        """Calcs the PMF of observing n cells spiking in the same time bin,
        as well as the PMF for indep cells (shuffled codes)"""
        if self.niswasNone:
            self.nis = random.sample(self.cs.nis, nbits) # randomly sample nbits of the nis
            self.nis.sort()
            self.nbits = nbits
        else:
            self.nis = self.cs.nis
            self.nbits = len(self.nis)
        self.words = {}
        self.nspiking = {}
        self.pnspiking = {}
        self.bins = {}

        for shufflecodes in (False, True):
            self.words[shufflecodes] = self.get_intcodes(nis=self.nis, shufflecodes=shufflecodes)
            # collect observances of the number of cells spiking for each pop code time bin
            self.nspiking[shufflecodes] = [ np.binary_repr(word).count('1') for word in self.words[shufflecodes] ] # convert the word at each time bin to binary, count the number of 1s in it. np.binary_repr() is a bit faster than using neuropy.Core.bin()
            self.pnspiking[shufflecodes], self.bins[shufflecodes] = histogram(self.nspiking[shufflecodes],
                                                                              bins=arange(self.nneurons+1),
                                                                              normed='pmf') # want all probs to add to 1, not their area, so use pmf

        assert (self.bins[False] == self.bins[True]).all() # paranoid, just checking
        self.bins = self.bins[False] # since they're identical, get rid of the dict and just keep one
        assert approx(self.pnspiking[False].sum(), 1.0), 'total observed probs: %f' % self.pnspiking[False].sum()
        assert approx(self.pnspiking[True].sum(), 1.0), 'total indep probs: %f' % self.pnspiking[True].sum()

        return self

    def plot(self, nbits=CODEWORDLEN, xlim=(-0.5, 15.5), ylim=(10**-6, 10**0)):
        """Plots nspikingPMF, for both observed and shuffled (forcing independence) codes"""

        try: self.pnspiking, self.bins
        except AttributeError: self.calc(nbits=nbits)

        f = figure()
        a = f.add_subplot(111)
        a.hold(True)
        a.plot(self.bins, self.pnspiking[False], 'r.-')
        a.plot(self.bins, self.pnspiking[True], 'b.-')
        titlestr = ''#'PMF of observing n cells spiking in the same time bin'
        titlestr += '\n%s' % lastcmd()
        if self.niswasNone:
            titlestr += '\nnis: %r' % self.nis
        a.set_title(titlestr)
        a.legend(('observed', 'indep (shuffled)'))
        a.set_yscale('log')
        if xlim:
            a.set_xlim(xlim)
        if ylim:
            a.set_ylim(ylim)
        gcfm().frame.SetTitle(lastcmd())
        a.set_xlabel('number of spiking cells in a bin')
        a.set_ylabel('probability')

        self.f = f
        self.a = a
        return self


class NetstateScatter(BaseNetstate):
    """Netstate scatter analysis object. See Schneidman Figures 1f and 2a"""
    def calc(self, nbits=None, model='both', radius=None, shufflecodes=False, algorithm='CG'):
        """Calculates the expected probabilities, assuming a model in ['indep', 'ising', 'both'],
        of all possible population codes vs their observed probabilities.
        self's nis are treated in LSB to MSB order"""
        self.nbits = nbits
        if self.nbits == None:
            self.nbits = CODEWORDLEN
        self.model = model
        self.radius = radius
        self.shufflecodes = shufflecodes
        self.algorithm = algorithm

        if self.niswasNone: # nis weren't specified in __init__
            self.nis = random.sample(self.cs.nis, self.nbits) # randomly sample nbits of the nis
            self.nis.sort()
        else: # nis were specified in __init__
            self.nis = self.cs.nis
            self.nbits = min(len(self.nis), self.nbits) # make sure nbits isn't > len(nis)

        self.intcodes = self.get_intcodes(nis=self.nis, shufflecodes=self.shufflecodes)
        self.pobserved, self.observedwords = histogram(self.intcodes, bins=arange(2**self.nbits), normed='pmf')
        if self.model == 'indep':
            self.pexpected, self.expectedwords = self.intcodesFPDF(nis=self.nis) # expected, assuming independence
        elif self.model == 'ising':
            ising = self.ising(nis=self.nis, radius=self.radius, algorithm=self.algorithm) # returns a maxent Ising model
            self.pexpected = ising.p # expected, assuming maxent Ising model
            self.expectedwords = ising.intsamplespace
        elif self.model == 'both':
            ising = self.ising(nis=self.nis, radius=self.radius, algorithm=self.algorithm) # returns a maxent Ising model
            self.pexpected = ising.p # expected, assuming maxent Ising model
            self.expectedwords = ising.intsamplespace
            self.pindepexpected = self.intcodesFPDF(nis=self.nis)[0] # expected, assuming independence
        else:
            raise ValueError, 'Unknown model %r' % self.model
        assert (self.observedwords == self.expectedwords).all() # make sure we're comparing apples to apples
        return self

    def plot(self, model='both', scale='freq',
             xlim=(10**-4, 10**2), ylim=(10**-11, 10**2), color=False):
        """Scatterplots the expected probabilities of all possible population codes (y axis)
        vs their observed probabilities (x axis). nis are in LSB to MSB order"""

        try: self.pobserved, self.pexpected
        except AttributeError: self.calc(model=model)

        f = figure()
        a = f.add_subplot(111)
        lo = min(xlim[0], ylim[0])
        hi = max(xlim[1], ylim[1])
        a.plot((lo, hi), (lo, hi), 'b-') # plot a y=x line
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
        if scale == 'freq':
            norm = self.tres / 1e6 # convert scale to patter freq in Hz
        elif scale == 'prob':
            norm = 1 # leave scale as pattern probabilities
        else:
            raise ValueError, 'Unknown scale %r' % scale
        self.norm = norm

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

        colorguide = ''
        if self.model == 'both': # plot the indep model too, and plot it first
            a.loglog(self.pobserved/norm, self.pindepexpected/norm, color='blue', marker='.', linestyle='None')
            colorguide = ' red: ising\n' + \
                         'blue: indep\n'
        # plot whichever model was specified
        if color:
            a.loglog(pobserved/norm, pexpected/norm, '.', color='black') # plots what's left in black
            a.loglog(pobserved4/norm, pexpected4/norm, '.', color='magenta')
            a.loglog(pobserved3/norm, pexpected3/norm, '.', color='blue')
            a.loglog(pobserved2/norm, pexpected2/norm, '.', color=(0, 1, 0))
            a.loglog(pobserved1/norm, pexpected1/norm, '.', color='red')
            colorguide = '    red: 1 spike patterns\n' + \
                         '  green: 2 spike patterns\n' + \
                         '   blue: 3 spike patterns\n' + \
                         'magenta: 4 spike patterns\n' + \
                         '  black: other patterns  \n'
        else:
            a.loglog(self.pobserved/norm, self.pexpected/norm, 'r.')
        '''
        a.plot(pobserved/norm, pexpected/norm, 'k.')
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
        titletext = lastcmd()
        if self.niswasNone:
            titletext += '\nneurons: %s' % self.nis # + missingcodetext)
        a.set_title(titletext)
        if scale == 'freq':
            labelend = 'state frequency (Hz)'
        elif scale == 'prob':
            labelend = 'state probability'
        a.set_xlabel('observed ' + labelend)
        a.set_ylabel('predicted ' + labelend)
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        if self.model =='both':
            DJSstring = '(%.4f, %.4f)' % (DJS(self.pobserved, self.pindepexpected), DJS(self.pobserved, self.pexpected))
        else:
            DJSstring = '%.4f' % DJS(self.pobserved, self.pexpected)
        a.text(0.99, 0.01, ('%s' +
                            '%.1f%% missing\n' +
                            'DJS=%s') % (colorguide, percentmissing, DJSstring), # add stuff to bottom right of plot
                            transform = a.transAxes,
                            horizontalalignment='right',
                            verticalalignment='bottom')

        self.f = f
        self.a = a
        return self

    def _onmotion(self, event):
        """Called during mouse motion over scatterplot figure. Pops up the corresponding
        population code word and its int representation when hovering over a neuron scatter point"""
        if event.xdata != None and event.ydata != None: # if mouse is inside the axes
            i  = approx(event.xdata, self.pobserved/self.norm, rtol=1e-1, atol=0).nonzero()[0] # find for what indices (if any) xdata == pobserved
            ii = approx(event.ydata, self.pexpected[i]/self.norm, rtol=1e-1, atol=0).nonzero()[0] # for those above, find for what index (if any) ydata == pexpected
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
                tip += '\npattern counts: %r' % [ (self.intcodes == intcode).sum() for intcode in intcodes ]
                tip += '\npattern freqs (Hz): %s' % repr([ '%.3g' % (p / self.tres * 1e6) for p in self.pobserved[codeis] ]).replace('\'', '')
                tip += '\n(pobserved, pexpected): %s' % repr(zip([ '%.3g' % val for val in self.pobserved[codeis] ], [ '%.3g' % val for val in self.pexpected[codeis] ])).replace('\'', '')
                tip += '\npobserved / pexpected: %s' % repr([ '%.3g' % (o/float(e)) for o, e in zip(self.pobserved[codeis], self.pexpected[codeis]) ]).replace('\'', '')
                tip += '\npexpected / pobserved: %s' % repr([ '%.3g' % (e/float(o)) for o, e in zip(self.pobserved[codeis], self.pexpected[codeis]) ]).replace('\'', '')
                self.tooltip.SetTip(tip) # update the tooltip
                self.tooltip.Enable(True) # make sure it's enabled
            else:
                self.tooltip.Enable(False) # disable the tooltip
        else: # mouse is outside the axes
            self.tooltip.Enable(False) # disable the tooltip


class NetstateI2vsIN(BaseNetstate):
    """Netstate I2/IN vs IN (fraction of pairwise correlated entropy vs all correlated entropy) analysis.
    See Schneidman fig 2c"""
    def calc(self, N=10, ngroups=15):
        """Computes I2/IN vs IN, for ngroups of cells.
        This shows you what fraction of network correlation is accounted for by the maxent pairwise model.
        Plotting Icond-indep is different from S1 (see below),
        and sounds annoying and not worth it (see methods in Schneidman 2006)"""
        self.N = N
        self.ngroups = ngroups

        self.niss = nCrsamples(objects=self.neurons.keys(),
                          r=self.N, # pick N neurons at random
                          nsamples=self.ngroups) # do it ngroups times
        I2s = []
        INs = []
        pd = wx.ProgressDialog(title='I2vsIN() progress', message='', maximum=self.ngroups, # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for groupi, nis in enumerate(self.niss):
            p1 = asarray(self.intcodesFPDF(nis=nis)[0]) # indep model
            p2 = self.ising(nis=nis).p # expected, assuming maxent Ising model
            pN = asarray(self.intcodesPDF(nis=nis)[0]) # observed word probs
            S1 = entropy_no_sing(p1) # ignore any singularities
            S2 = entropy_no_sing(p2)
            SN = entropy_no_sing(pN)
            IN = S1 - SN
            I2 = S1 - S2
            I2s.append(I2 / self.tres * 1e6) # convert to bits/sec
            INs.append(IN / self.tres * 1e6)
            cont, skip = pd.Update(groupi, newmsg='groupi = %d' % (groupi+1))
            if not cont:
                pd.Destroy()
                return
        pd.Destroy()
        self.I2s = asarray(I2s)
        self.INs = asarray(INs)
        self.I2divIN = self.I2s / self.INs

        return self

    def plot(self, xlim=(0.0, None), ylim=(0.0, 1.0)):
        """Plots I2/IN vs IN"""

        try: self.I2s, self.INs, self.I2divIN
        except AttributeError: self.calc()

        f = figure()
        gcfm().frame.SetTitle(lastcmd())
        a = f.add_subplot(111)
        a.plot(self.INs, self.I2divIN, 'r.')
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xlabel('IN (bits / sec)')
        a.set_ylabel('I2 / IN')
        a.set_title('%s' % lastcmd())
        a.text(0.99, 0.01, 'mean=%.3f, std=%.3f' % (self.I2divIN.mean(), self.I2divIN.std()), # add mean and std to bottom right
               transform=a.transAxes,
               horizontalalignment='right',
               verticalalignment='bottom')

        self.f = f
        self.a = a
        return self


class NetstateDJSHist(BaseNetstate):
    """Jensen-Shannon histogram analysis. See Schneidman 2006 figure 2b"""
    def calc(self, nbits=CODEWORDLEN, ngroups=5, models=['ising', 'indep'],
                   shufflecodes=False, algorithm='CG'):
        """Calculates Jensen-Shannon divergences and their ratios
        for ngroups random groups of cells, each of length nbits."""
        self.nbits = nbits
        self.ngroups = ngroups
        self.models = models
        self.shufflecodes = shufflecodes
        self.algorithm = algorithm

        self.niss = [] # list of lists, each sublist is a group of neuron indices
        self.DJSs = {} # hold the Jensen-Shannon divergences for different models and different groups of neurons
        for model in  self.models:
            self.DJSs[model] = [] # init a dict with the model names as keys, and empty lists as values
        pd = wx.ProgressDialog(title='DJShist progress', message='', maximum=self.ngroups*len(self.models),
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME) # create a progress dialog
        for groupi in range(self.ngroups): # for each group of nbits cells
            nis = random.sample(self.cs.nis, self.nbits) # randomly sample nbits of the Netstate Codes' nis attrib
            self.niss.append(nis)
            for modeli, model in enumerate(self.models): # for each model, use the same nis
                so = self.r.ns_scatter(experiments=self.e, nis=nis, kind=self.kind, tres=self.tres, phase=self.phase) # netstate scatter object
                so.calc(nbits=self.nbits, model=model, shufflecodes=self.shufflecodes, algorithm=self.algorithm)

                self.DJSs[model].append(DJS(so.pobserved, so.pexpected))
                cont, skip = pd.Update(groupi*len(self.models)+modeli, newmsg='groupi = %d\nmodel = %s' % (groupi, model))
                if not cont:
                    pd.Destroy()
                    return
        pd.Destroy()

        # now find the DJSratios between the two models, for each group of neurons
        if len(self.models) == 2: # do it only if there's 2 models, otherwise it's indeterminate which two to take ratio of
            self.DJSratios = asarray(self.DJSs.values()[1]) / asarray(self.DJSs.values()[0]) # 2nd model as a ratio of the 1st

        return self

    def plot(self, ngroups=5, logrange=(-3.667, -0.333), nbins=50, publication=False):
        """Plots histogram DJSs and DJSratios in logspace"""

        try: self.niss, self.DJSs
        except AttributeError: self.calc(ngroups=ngroups)

        x = np.logspace(start=logrange[0], stop=logrange[1], num=nbins, endpoint=True, base=10.0)
        n = {} # stores a list of the bin heights in a separate key for each model
        for model in self.models:
            n[model] = histogram(self.DJSs[model], bins=x, normed=False)[0]
        color = {'indep': 'b', 'ising': 'r'} # dict that maps from model name to color
        barwidths = list(diff(x)) # each bar will have a different width, convert to list so you can append
        # need to add one more entry to barwidth to the end to get nbins of them:
        barwidths.append(0) # don't display the last one
        logbinwidth = (logrange[1]-logrange[0]) / float(nbins)
        #barwidths.append(10**(logrange[1]+logbinwidth) - x[-1]) # should be exactly correct

        # plot DJSs of all models on the same axes
        f1 = figure()
        a1 = f1.add_subplot(111)
        #a1.hold(True)
        bars = {}
        heights = {}
        for model in self.models:
            heights[model] = n[model] / float(self.ngroups * logbinwidth)
            bars[model] = a1.bar(left=x, height=heights[model], width=barwidths, color=color[model],
                                 edgecolor=color[model])
        a1.set_xscale('log', basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
        a1.set_xlim(xmin=10**logrange[0], xmax=10**logrange[1])
        gcfm().frame.SetTitle(lastcmd())
        a1.set_title('%s' % lastcmd())
        if publication:
            a1.set_xticklabels(['', '0.001', '0.01', '0.1', '']) # hack!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for label in a1.get_xticklabels():
                label.set_size(30)
            for label in a1.get_yticklabels():
                label.set_size(30)
            a1.legend([ bars[model][0] for model in self.models ], ['pairwise', 'independent'], prop=mpl.font_manager.FontProperties(size=20) ) # grab the first bar for each model, label it with the model name
        else:
            a1.set_ylabel('number of groups of %d cells' % self.nbits)
            a1.set_xlabel('DJS (bits)')
            a1.set_ylabel('probability density (1 / log10(DJS))')
            a1.legend([ bars[model][0] for model in self.models ], ['pairwise', 'independent'])
        # plot DJSratios
        if len(self.models) == 2:
            f2 = figure()
            a2 = f2.add_subplot(111)
            nratios = histogram(self.DJSratios, bins=x, normed=False)[0] # bin heights for the DJSratios
            a2.bar(left=x, height=nratios, width=barwidths, color='g', edgecolor='g')
            a2.set_xscale('log', basex=10) # need to set scale of x axis AFTER bars have been plotted, otherwise autoscale_view() call in bar() raises a ValueError for log scale
            gcfm().frame.SetTitle(lastcmd())
            a2.set_title('Jensen-Shannon divergence ratios histogram\n%s' % lastcmd())
            a2.set_ylabel('number of groups of %d cells' % self.nbits)
            a2.set_xlabel('DJS ratio (%s / %s)' % (self.models[1], self.models[0]))

        self.f = {1:f1, 2:f2}
        self.a = {1:a1, 2:a2}
        return self


class NetstateS1INvsN(BaseNetstate):
    """Analysis of uncorrelated entropy and reduction by correlated entropy for increasing network size N"""
    def calc(self, minN=4, maxN=15, maxnsamples=10):
        """Calculates the average independent (uncorrelated) cell entropy S1
        and average network multi-information IN (IN = S1 - SN) vs network size N.
        IN is how much the correlated entropy reduces the total entropy of the system.
        For each network size up to maxN, averages S1 and IN over maxnsamples (or less if that many aren't possible)
        number of groups at each value of N"""
        self.minN = minN
        self.maxN = maxN
        self.maxnsamples = maxnsamples

        self.S1ss = [] # as f'n of N
        self.INss = []
        self.N = range(self.minN, self.maxN+1) # network sizes from minN up to maxN
        #tstart = time.clock()
        self.nsamples = [ min(nCr(self.nneurons, r), self.maxnsamples) for r in self.N ] # nsamples as a f'n of N. For each value of N, take up to maxnsamples of all the other neurons, if that many are even possible

        pd = wx.ProgressDialog(title='S1INvsN progress', message='', maximum=sum(self.nsamples), # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)
        for ni, n in enumerate(self.N): # for all network sizes
            # get a list of lists of neuron indices
            niss = nCrsamples(objects=self.neurons.keys(),
                              r=n, # pick n neurons
                              nsamples=self.nsamples[ni] ) # do it at most maxnsamples times
            S1s = []
            INs = []
            for nisi, nis in enumerate(niss):
                cont, skip = pd.Update(ni*self.maxnsamples+nisi, newmsg='N = %d\nsamplei = %d' % (n, nisi))
                if not cont:
                    pd.Destroy()
                    return
                #t2 = time.clock()
                p1 = asarray(self.intcodesFPDF(nis=nis)[0]) # indep model
                pN = asarray(self.intcodesPDF(nis=nis)[0]) # observed word probs
                #print 'calcing ps took: %f sec' % (time.clock()-t2)
                S1 = entropy_no_sing(p1) # ignore any singularities
                SN = entropy_no_sing(pN)
                assert S1 > SN or approx(S1, SN), 'S1 is %.20f, SN is %.20f' % (S1, SN) # better be, indep model assumes the least structure
                IN = S1 - SN
                #print S1, SN, IN
                S1s.append(S1 / self.tres * 1e6) # convert to bits/sec
                INs.append(IN / self.tres * 1e6)
            self.S1ss.append(S1s)
            self.INss.append(INs)
        pd.Destroy()
        self.S1mean = [ asarray(S1s).mean() for S1s in self.S1ss ]
        self.S1std = [ asarray(S1s).std() for S1s in self.S1ss ]
        self.S1sem = asarray(self.S1std) / sqrt(asarray(self.nsamples))
        self.INmean = [ asarray(INs).mean() for INs in self.INss ]
        self.INstd = [ asarray(INs).std() for INs in self.INss ]
        self.INsem = asarray(self.INstd) / sqrt(asarray(self.nsamples))

        return self

    def plot(self, xlim=(1e0, 1e3), ylim=(1e-2, 1e4)):
        """Plots the average independent (uncorrelated) cell entropy S1
        and average network multi-information IN (IN = S1 - SN) vs network size N."""

        try: self.S1ss
        except AttributeError: self.calc()

        f = figure()
        gcfm().frame.SetTitle(lastcmd())
        a = f.add_subplot(111)
        a.hold(True)
        for n, S1s in zip(self.N, self.S1ss): # plot all the samples before plotting the means with errorbars
            a.plot([n]*len(S1s), S1s, '_', markersize=4, color='lightblue')
        for n, INs in zip(self.N, self.INss):
            a.plot([n]*len(INs), INs, '_', markersize=4, color='pink')
        S1line = a.errorbar(self.N, self.S1mean, yerr=self.S1sem, fmt='b.')[0]
        INline = a.errorbar(self.N, self.INmean, yerr=self.INsem, fmt='r.')[0]
        # do least squares polynomial fit in log10 space
        mS1, bS1 = sp.polyfit(log10(self.N), log10(self.S1mean), 1) # returns slope and y intercept
        mIN, bIN = sp.polyfit(log10(self.N), log10(self.INmean), 1)
        xintersect = (bIN - bS1) / (mS1 - mIN)
        x = array([-1, 3]) # define x in log10 space, this is really [0.1, 1000]
        self.yS1 = mS1*x + bS1 # y = mx + b
        self.yIN = mIN*x + bIN
        a.plot(10.0**x, 10.0**self.yS1, 'b-') # take their power to make up for both the x and y scales being log
        a.plot(10.0**x, 10.0**self.yIN, 'r-')
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xlabel('Number of cells')
        a.set_ylabel('bits / sec')
        a.set_title('S1 & IN vs N\n%s' % lastcmd())
        a.legend((S1line, INline), ('S1, slope=%.3f' % mS1, 'IN, slope=%.3f' % mIN), loc='lower right')
        a.text(0.99, 0.98, 'Nc=%d' % np.round(10**xintersect), # add text box to upper right corner of axes
                           transform = a.transAxes,
                           horizontalalignment = 'right',
                           verticalalignment = 'top')

        self.f = f
        self.a = a
        return self

class NetstateNNplus1(BaseNetstate):
    """Analysis of amount of mutual information between N cells and the N+1th cell"""
    def calc(self, Nplus1s=None, maxN=15, maxnsamples=10):
        """Calculates Schneidman Figure 5b. Averages over as many as
        maxnsamples different groups of N cells for each N+1th cell in Nplus1s,
        all done for different values of N up to maxN"""
        if Nplus1s == None: # list of all indices of neurons that will be treated as the N+1th neuron
            Nplus1s = self.cs.nis
        else:
            Nplus1s = toiter(Nplus1s)
        nNplus1s = len(Nplus1s)
        dims = (maxN, nNplus1s, maxnsamples)
        mask = np.zeros(dims) # this will be converted to an array of Falses
        IdivS = np.ma.array(mask, mask=mask, fill_value=666) # masked array that holds the mutual info between N and N+1th cells, as a ratio of the N+1th cell's entropy. Index like: IdivS[ni, Nplus1i, samplei], ie group size, N+1th cell you're comparing to, and number of samples of size N taken from the possible combs
        self.N = range(1, maxN+1) # cell group size, excluding the N+1th neuron. This will be the x axis in the plot
        #self.N=[15]#.reverse() # for fun and pleasure
        nsamples = [ min(maxnsamples, nCr(nNplus1s-1, r)) for r in self.N ] # take up to maxnsamples of all the other neurons, if that many even exist (for the lower N values, the constraint may end up being the total number of possible combinations of cells), for each N+1th cell. Taken from nNplus1s-1 cuz you always have to exclude an N+1th neurons
        for ni, n in enumerate(self.N): # for all group sizes
            IdivS.mask[ni, :, nsamples[ni]::] = True # mask out the sampleis that are out of range for this value of N, if any
        maximum = nNplus1s*sum(nsamples)
        pd = wx.ProgressDialog(title='NNplus1 progress', message='', maximum=maximum, # create a progress dialog
                               style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME)

        # get the binary array for the whole population, then index into it appropriately in the sample loop, find the corresponding integer codes, and feed it to MIbinarray, so you don't have to unnecessarily re-generate it on every iteration
        nis = self.cs.nis
        nis2niis = self.cs.nis2niis
        counter = 0 # counts inner loop for the progress dialog
        for ni, n in enumerate(self.N):
            for Nplus1i, Nplus1 in enumerate(Nplus1s): # for each N+1th neuron to compare to
                mii = nis2niis(Nplus1)
                niscopy = copy(nis) # make a copy of neuron indices
                niscopy.remove(Nplus1) # keep just the indices of all the other neurons
                samples = nCrsamples(niscopy, n, nsamples[ni]) # returns nsamples random unique choices of n items from niscopy
                for samplei, sample in enumerate(samples): # collect nsamples different combinations of the N other cells
                    niis = np.array([ nis2niis(s) for s in toiter(sample) ]) # most of the time (for n>1), sample will be a sequence of nis. Build an array of niis out of it to use as indices into the binary code array. Sometimes (for n=1) sample will be a scalar, hence the need to push it through toiter()
                    IdivS[ni, Nplus1i, samplei] = MIbinarrays(Nbinarray=self.cs.c[niis], Mbinarray=self.cs.c[mii]).IdivS # do it
                    cont, skip = pd.Update(counter, newmsg='N: %d; N+1th neuron: %d; samplei: %d' % (n, Nplus1, samplei))
                    if not cont:
                        pd.Destroy()
                        return
                    counter += 1
        pd.Destroy()
        self.IdivS = IdivS.reshape(maxN, nNplus1s*maxnsamples) # reshape such that you collapse all Nplus1s and samples into a single dimension (columns). The N are still in the rows
        #logIdivS = log10(IdivS)
        self.IdivSmeans = self.IdivS.mean(axis=1) # average over all Nplus1s and all samples. Values that are masked are ignored
        self.IdivSstds = self.IdivS.std(axis=1) # find stdev for the same
        assert self.IdivSmeans.shape == (maxN,)
        assert self.IdivSstds.shape == (maxN,)
        self.IdivSsems = self.IdivSstds / sqrt(asarray(nsamples)*nNplus1s)

        return self

    def plot(self, maxN=15, maxnsamples=10, xlim=(10**log10(0.9), 1e3), ylim=(1e-3, 1e1)):
        """Plots the figure with error bars"""

        try: self.IdivS
        except AttributeError: self.calc(maxN=maxN, maxnsamples=maxnsamples)

        f = figure()
        gcfm().frame.SetTitle(lastcmd())
        a = f.add_subplot(111)
        a.hold(True)
        for n, row in zip(self.N, self.IdivS): # underplot the samples for each value of N
            a.plot([n]*len(row), row, '_', markersize=4, color='deepskyblue')
        a.errorbar(self.N, self.IdivSmeans, yerr=self.IdivSsems, fmt='b.') # now plot the means and sems
        # do some linear regression in log10 space
        m, b = sp.polyfit(log10(self.N), log10(self.IdivSmeans), 1) # returns slope and y intercept
        x = array([log10(0.9), 3]) # define x in log10 space, this is really [0.9, 1000]
        y = m*x + b
        xintersect = (0-b) / m # intersection point of regression line with y=1=10**0 line
        plot(10.0**x, 10.0**y, 'b-') # raise them to the power to make up for the fact that both the x and y scales will be log
        plot(10.0**x, [1e0]*2, 'r--') # plot horizontal line at y=1
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim(xlim)
        a.set_ylim(ylim)
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

        self.f = f
        self.a = a
        return self
        '''
        # plot the distributions of IdivS
        for ni, n in enumerate(self.N):
            f = figure()
            gcfm().frame.SetTitle('%s IdivS distrib for N=%d' % (lastcmd(), n))

            notmaskedis = self.IdivS[ni].mask==False # indexes the non-masked entries in IdivS, for this ni

            a1 = f.add_subplot(211) # axes with linear bins
            heights, bins = histogram(self.IdivS[ni, notmaskedis], bins=arange(0, 1, 0.02))
            barwidth = bins[1]-bins[0]
            a1.bar(left=bins, height=heights, width=barwidth, bottom=0, color='k')
            #a1.set_xlabel('mutualinfo(N, N+1th) / entropy(N+1th)')
            a1.set_ylabel('count')
            a1.set_title('IdivS distrib for N=%d' % n)

            a2 = f.add_subplot(212) # axes with log bins
            start = log10(0.001)
            stop = log10(1)
            bins = np.logspace(start=start, stop=stop, num=50, endpoint=True, base=10.0)
            heights, bins = histogram(self.IdivS[ni, notmaskedis], bins=bins)
            barwidth = list(diff(bins)) # each bar will have a different width, convert to list so you can append
            # need to add one more entry to barwidth to the end to get nbins of them:
            #barwidth.append(barwidth[-1]) # not exactly correct
            logbinwidth = (log10(stop)-log10(stop)) / float(len(bins))
            barwidth.append(10**(log10(stop)+logbinwidth) - stop) # should be exactly correct
            a2.bar(left=bins, height=heights, width=barwidth, bottom=0, color='k')
            a2.set_xscale('log')
            a2.set_xlabel('mutualinfo(N, N+1th) / entropy(N+1th)')
            a2.set_ylabel('count')
        '''

class NetstateCheckcells(BaseNetstate):
    """Analysis of how activity rates of each cell in the population vary with
    the overall amount of activity in the rest of the population"""
    def _calc(self, ni=None, othernis=None, shufflecodes=False):
        """Calculates the joint pdf of cell ni activity and the number of cells in
        othernis being active at the same time. ni should not be in othernis"""
        assert ni not in othernis
        nis2niis = self.cs.nis2niis
        nii = nis2niis(ni)
        otherniis = nis2niis(othernis)
        nothers = len(othernis)

        nicode = self.cs.c[nii] # 0s and 1s, this picks out the row in the binary code array that corresponds to ni
        othercodes = self.cs.c[otherniis]
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

    def calc(self, nis=None, othernis=None, nothers=None, nsamples=10, shufflecodes=False):
        """Calcs the probability of each cell (in nis) being active vs. the number of
        other active cells (in the Recording) at that time. For each ni, calcs an average over
        nsamples, each being a different sample of nothers from othernis.
        See Schneidman figure 5c"""
        if nis == None:
            self.nis = self.cs.nis
        else:
            self.nis = toiter(nis)
        if othernis == None:
            self.othernis = self.cs.nis
        else:
            self.othernis = toiter(othernis)
        if nothers == None:
            less = 0
            while nCr(len(self.othernis), len(self.othernis)-less) < nsamples:
                less += 1
            nothers = len(self.othernis) - 1 - less # -1 to remove ni, -less again to allow for at least nsamples combos of othernis
        self.N = arange(nothers+1)

        try: self.jpdfss
        except AttributeError:
            # init dicts to store jpdfs and other stuff in
            self.jpdfss = {}
            self.jpdfmeans = {}
            self.jpdfstds = {}
            self.jpdfsems = {}

        for ni in self.nis:
            try:
                self.jpdfss[ni]
            except KeyError:
                othernis = copy(self.othernis) # don't modify the original
                try:
                    othernis.remove(ni) # all the other possible nis, excluding the current ni
                except ValueError: # ni isn't in othernis, nothing to remove
                    pass
                otherniss = nCrsamples(objects=othernis, r=nothers, nsamples=nsamples) # get nsamples unique random samples of length nothers from othernis
                jpdfs = []
                for othernis in otherniss: # collect jpdfs across all random samples
                    jpdf = self._calc(ni=ni, othernis=othernis, shufflecodes=shufflecodes)
                    jpdfs.append(jpdf)
                jpdfs = asarray(jpdfs) # this is an nsamples x 2 x (1+nothers) matrix
                self.jpdfss[ni] = jpdfs
                self.jpdfmeans[ni] = jpdfs.mean(axis=0) # find the mean jpdf across all nsamples jpdfs
                self.jpdfstds[ni] = jpdfs.std(axis=0) # find the stdev across all nsamples jpdfs
                self.jpdfsems[ni] = self.jpdfstds[ni] / sqrt(nsamples)
        return self

    def plot(self, nis=None, nothers=None, nsamples=10):
        """Plots the desired neurons so you can see if they behave like check cells"""
        try: self.jpdfss
        except AttributeError: self.calc(nis=nis, nothers=nothers, nsamples=nsamples)
        try: self.f
        except AttributeError: self.f = {}
        try: self.a
        except AttributeError: self.a = {}
        if nis == None:
            nis = self.nis
        else:
            nis = toiter(nis)
        for ni in nis:
            f = figure()
            gcfm().frame.SetTitle('%s for ni=%d' % (lastcmd(), ni))
            a = f.add_subplot(111)
            a.hold(True)
            # plot all the samples first
            for jpdf in self.jpdfss[ni]: # iter over the hyperrows
                a.plot(self.N, jpdf[1], '_', markersize=4, color='grey') # marginal pdf of getting a 1 for the check cell
            # plot the stdevs, means, and sems of marginal pdf of getting a 1 for the check cell
            #a.errorbar(self.N, self.jpdfmeans[ni][1], yerr=self.jpdfstds[ni][1], fmt=None, capsize=0, ecolor='grey')
            a.errorbar(self.N, self.jpdfmeans[ni][1], yerr=self.jpdfsems[ni][1], fmt='k.-')
            a.set_ylim(ymin=0, ymax=1)

            titlestr = '%s\nni=%d' % (lastcmd(), ni)
            titlestr += ', nsamples=%d' % nsamples
            a.set_title(titlestr)
            a.set_xlabel('Number of other active cells')
            a.set_ylabel('Probability of cell ni being active')

            self.f[ni] = f
            self.a[ni] = a
        return self


class NetstateTriggeredAverage(BaseNetstate):
    """Analysis that reverse correlates the occurence of a specific netstate
    to the stimuli in the Experiments in this Recording to build up a netstate
    triggered average"""
    def cut(self, trange):
        """Cuts network state word times according to trange"""
        lo, hi = self.wordts.searchsorted([trange[0], trange[1]]) # returns indices where tstart and tend would fit in wordts
        if trange[1] == self.wordts[min(hi, len(self.wordts)-1)]: # if tend matches a word time (protect from going out of index bounds when checking)
            hi += 1 # inc to include a word time if it happens to exactly equal tend. This gives us end inclusion
            hi = min(hi, len(self.wordts)) # limit hi to max slice index (==max value index + 1)
        cutwordts = self.wordts[lo:hi] # slice it
        return cutwordts

    def calc(self, intcode=None, nt=9, ti0=-4):
        """Calculate the network state triggered average for word intcode, using nt revcorr timepoints,
        starting at revcorr timepoint index ti0.

        For now, this uses the Codes object created across the entire Recording
        """
        self.intcode = intcode
        self.nt = nt # number of revcorr timepoints
        self.ti0 = ti0
        self.tis = range(ti0, ti0+nt, 1) # revcorr timepoint indices, can be -ve. these will be multiplied by the movie frame time
        words = binarray2int(self.cs.c)
        i = (words == intcode)
        assert i.any(), 'netstate intcode %d never occured'
        self.wordts = self.cs.t[i] # netstate intcode word times

        if self.e == None:
            self.e = self.r.e # if no specific experiments were specified to revcorr to in __init__, revcorr to all of them
        self.frames = {} # accumulates movie frames separately for each timepoint across all Experiments' movies
        for ti in self.tis:
            self.frames[ti] = []
        tstart = time.clock()
        pd = wx.ProgressDialog(title='NSTA progress: loading movies, collecting frames', message='',
                               maximum=len(self.e), style=1) # create a progress dialog
        for ei, e in enumerate(self.e.values()):
            cont, skip = pd.Update(ei-1, newmsg='experiment %d\nelapsed: %.1fs' % (e.id,
                                   time.clock()-tstart))
            if not cont:
                pd.Destroy()
                return
            # for now, we're using the Codes object created across the entire Recording. It might be slightly more correct to generate a separate codes object for each Experiment. That way, the wordts for would be aligned to the start of each Experiment, as opposed to the start of the Recording, as they are now.
            cutwordts = self.cut(e.trange) # get wordts that were active during this experiment. This isn't really necessary is it??????????????????????????? Might speed things up for the next searchsorted call, since cutwordts is shorter than wordts
            rcdini = e.din[:, 0].searchsorted(cutwordts) - 1 # revcorr dini. Find where the cutwordts times fall in the din, dec so you get indices that point to the most recent din value for each cutwordt
            #self.din = e.din[rcdini, 1] # get the din (frame indices) at the rcdini
            # for now, only do revcorr if experiment.stims has only one entry. Stims no longer exists in dimstim >= 0.16 anyway
            assert len(e.stims) == 1
            movie = e.stims[0] # this will have to change for dimstim >= 0.16
            movie.load() # ensure movie is loaded
            data = movie.data
            try:
                self.width
            except AttributeError: # init stuff
                self.width = movie.data.shape[-1]
                self.height = movie.data.shape[-2] # dims are nframes, height, width
                self.sweeptimeMsec = movie.sweeptimeMsec # changes in dimstim >= 0.16
                self.ndinperframe = intround(movie.sweeptimeMsec / float(e.REFRESHTIME / 1000.)) # changes in dimstim >= 0.16
                self.regionwidthDeg = movie.regionwidthDeg
                self.regionheightDeg = movie.regionheightDeg
                self.origDeg = e.origDeg
            # assert that all movies are the same size and in the same spot. that way you can just accumulate frames in a single array
            assert self.width == movie.data.shape[-1] # dims are nframes, height, width
            assert self.height == movie.data.shape[-2]
            assert self.sweeptimeMsec == movie.sweeptimeMsec # changes in dimstim >= 0.16
            assert self.regionwidthDeg == movie.regionwidthDeg
            assert self.regionheightDeg == movie.regionheightDeg
            assert self.origDeg == e.origDeg
            for ti in self.tis:
                shiftedrcdini = rcdini - ti*self.ndinperframe # this can unintentionally introduce -ve valued indices at the left boundary, or out of range values at right boundary
                shiftedrcdini = shiftedrcdini[shiftedrcdini >= 0] # remove any -ve valued indices. Is this the most efficient way to do this?
                shiftedrcdini = shiftedrcdini[shiftedrcdini <= len(e.din)-1] # remove any out of range values
                frameis = e.din[shiftedrcdini, 1] # get the din values (frame indices) at the rcdini for this timepoint
                # in Cat 15, we erroneously duplicated the first frame of the mseq movies at the end, giving us one more frame (0 to 65535 for mseq32) than we should have had (0 to 65534 for mseq32). We're now using the correct movies, but the din for Cat 15 mseq experiments still have those erroneous frame indices (65535 and 16383 for mseq32 and mseq16 respectively), so we'll just ignore them for revcorr purposes.
                if movie.oname == 'mseq32':
                    frameis = frameis[frameis != 65535] # remove all occurences of 65535
                elif movie.oname == 'mseq16':
                    frameis = frameis[frameis != 16383] # remove all occurences of 16383
                frames = data.take(frameis.astype(np.int32), axis=0) # collect the relevant frames for this timepoint, take is much faster than direct indexing, but have to typecast indices to int32, maybe cuz this machine is 32bit?
                self.frames[ti].append(frames)
            pd.Destroy()


        # now that we're outside the experiment loop, get the mean for each timepoint across all experiments
        self.rf = zeros([self.nt, self.height, self.width], dtype=np.float64) # 3D matrix to store the NSTA at each timepoint. rf == 'receptive field'
        self.t = [ ti*intround(self.sweeptimeMsec) for ti in self.tis ]
        pd = wx.ProgressDialog(title='NSTA progress: averaging frames', message='',
                               maximum=self.nt, style=1) # create a progress dialog
        for tii, ti in enumerate(self.tis):
            cont, skip = pd.Update(tii-1, newmsg='timepoint: %dms\nelapsed: %.1fs' % (self.t[tii], time.clock()-tstart))
            if not cont:
                pd.Destroy()
                return
            self.frames[ti] = cat(tuple(self.frames[ti])) # need to concatenate all lists for this ti into a single array
            self.rf[tii] = mean_accum(self.frames[ti]) # this is much faster than frames.mean()
            #self.rf[ti] = mean_accum2(data, frameis)
        pd.Destroy()
        return self

    def plot(self, intcode=None, nt=10, ti0=-4, interp='nearest', normed=True, scale=2.0):
        """Plots the spatiotemporal RF as bitmaps in a wx.Frame"""
        try:
            self.rf
        except AttributeError:
            self.calc(intcode=intcode, nt=nt, ti0=ti0)
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
        self.rfframe = NetstateReceptiveFieldFrame(title=lastcmd(), rfs=[rf], intcodes=self.intcode, t=self.t, scale=scale)
        self.rfframe.Show()
        return self


class RecordingNetstate(BaseRecording):
    """Mix-in class that defines Netstate related Recording methods"""
    def ns_(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a BaseNetstate object"""
        return BaseNetstate(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_isinghist(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateIsingHist object"""
        return NetstateIsingHist(recording=self, experiments=experiments, kind=kind, tres=tres, phase=phase)
    def ns_nspikingpmf(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateNspikingPMF object"""
        return NetstateNspikingPMF(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_scatter(self, experiments=None, nis=None, radius=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateScatter object"""
        ns_so = NetstateScatter(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
        ns_so.calc(radius=radius)
        return ns_so
    def ns_i2vsin(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateI2vsIN object"""
        return NetstateI2vsIN(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_djshist(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateDJSHist object"""
        return NetstateDJSHist(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_s1invsn(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateS1INvsN object"""
        return NetstateS1INvsN(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_nnplus1(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateNNplus1 object"""
        return NetstateNNplus1(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_checkcells(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateCheckcells object"""
        return NetstateCheckcells(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)
    def ns_ta(self, experiments=None, nis=None, kind=CODEKIND, tres=CODETRES, phase=CODEPHASE):
        """Returns a NetstateTriggeredAverage object"""
        return NetstateTriggeredAverage(recording=self, experiments=experiments, nis=nis, kind=kind, tres=tres, phase=phase)


class Recording(RecordingRaster,
                RecordingCode,
                RecordingNetstate,
                BaseRecording):
    """Inherits all the Recording classes into a single Recording class"""
    pass
