"""Defines the Recording class and all of its support classes"""

# set the self.trange attribe in Base class for Recording, to be consistent with Neuron and Experiment

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
        firstexp = min(self.e.keys())
        lastexp = max(self.e.keys())
        self.trange = self.e[0].trange[firstexp], self.e[lastexp].trange[1] # start of the first experiment to end of the last one

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
    def __init__(self, recording, sortby='id'):
        self.r = recording
        self.t0 = self.r.trange[0]
        self.e = recording.e # dictionary
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
        gcfm().frame.SetTitle('r%d.raster(sortby=%s)' % (self.r.id, repr(self.sortby)))
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
        """Sorts self.neurons according to their attribute specified by self.sortby"""
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
                startlines = self.a.vlines(x=estart/1000.0, ymin=self.yrange[0], ymax=self.yrange[1], fmt='k--') # marks exp start, convert to ms
                startlines[0].set_color((0, 1, 0)) # set to bright green
            if left <= eend and eend <= left+width: # experiment end point is within view
                endlines = self.a.vlines(x=eend/1000.0, ymin=self.yrange[0], ymax=self.yrange[1], fmt='k--') # marks exp end, convert to ms
                endlines[0].set_color((1, 0, 0)) # set to bright red
        # plot the rasters
        for nii, neuron in enumerate(self.neurons):
            x = (neuron.cut((self.t0+left, self.t0+left+width)) - self.t0) / 1000.0 # make spike times always relative to t0, convert to ms
            self.a.vlines(x=x, ymin=nii, ymax=nii+1, fmt='k-')
        self.a.set_xlim(left/1000.0, (left+width)/1000.0) # convert from us to ms
    def panx(self, nsteps=None, left=None):
        """Pans the raster along the x axis"""
        self.a.lines=[] # first, clear all the vlines, this is easy but a bit innefficient, since we'll be redrawing most of the ones we just cleared
        if left != None: # use left
            self.plot(left=left, width=self.width)
        else: # use nsteps instead
            self.plot(left=self.left+self.width*nsteps, width=self.width)
        self.f.canvas.draw() # redraw the figure
    def zoomx(self, nsteps):
        """Zooms the raster along the x axis"""
        self.a.lines=[] # first, clear all the vlines, this is a bit innefficient, since we'll be redrawing most of the ones we just cleared
        self.plot(left=self.left+self.width*nsteps, width=self.width-2*self.width*nsteps)
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
                if estart < event.xdata and event.xdata < eend:
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
        #print key
        #import pdb; pdb.set_trace()
        if not event.guiEvent.ControlDown(): # wx dependent
            if key == wx.WXK_RIGHT:
                self.panx(+0.1)
            elif key == wx.WXK_LEFT:
                self.panx(-0.1)
            elif key == wx.WXK_UP:
                self.zoomx(+0.1)
            elif key == wx.WXK_DOWN:
                self.zoomx(-0.1)
            elif key == wx.WXK_NEXT: # PGDN (page right)
                self.panx(+1)
            elif key == wx.WXK_PRIOR: # PGUP (page left)
                self.panx(-1)
            elif key == wx.WXK_HOME: # go to start of first Experiment
                self.panx(left=0)
            elif key == wx.WXK_END: # go to end of last Experiment
                self.panx(left=self.r.trange[1]-self.t0-self.width)
        else: # Ctrl key is down, skip forward or back to next experiment marker
            if key == wx.WXK_LEFT:
                i = self.experimentmarkers.searchsorted(self.left, side='left') # current position of left edge of the window in experimentmarkers list
                i = max(0, i-1) # decrement by 1, do bounds checking
                self.panx(left=self.experimentmarkers[i])
            if key == wx.WXK_RIGHT:
                i = self.experimentmarkers.searchsorted(self.left, side='right') # current position of left edge of the window in experimentmarkers list
                i = min(i, len(self.experimentmarkers)-1) # bounds checking
                self.panx(left=self.experimentmarkers[i])

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


class Recording(RecordingRaster,
                BaseRecording):
    """Inherits all the Recording classes into a single Recording class"""
    pass
