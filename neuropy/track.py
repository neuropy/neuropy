"""Defines the Track class"""

from __future__ import division
from __future__ import print_function

import os
import StringIO

import numpy as np

import matplotlib as mpl
import pylab as pl
from pylab import get_current_fig_manager as gcfm

import core
from core import dictattr, TAB, td2usec, lastcmd, intround
from recording import Recording
from sort import TrackSort


class Track(object):
    """A track can have multiple recordings"""
    def __init__(self, path, animal=None):
        self.level = 2 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # string buffer to print tree hierarchy to
        self.path = path
        self.animal = animal
        if animal != None:
            # update parent animal's track dict, in case self wasn't loaded by its parent
            animal.tr[self.id] = self
        self.r = dictattr() # store recordings in a dictionary with attrib access

    def get_name(self):
        return os.path.split(self.path)[-1]
    
    name = property(get_name)

    def get_absname(self):
        """Absolute name including parent Animal, kind of like absolute path, but more
        abbreviated, as one would enter it at the IPython prompt"""
        if self.animal == None: # no parent animal
            return self.name
        else:
            return '.'.join((self.animal.name, self.name))

    absname = property(get_absname)

    def get_id(self):
        # replace first occurrence of 'tr' with nothing, keep the rest:
        id = self.name.lower().replace('tr', '', 1)
        if not id:
            raise NameError('Badly formatted track name: %s' % self.name)
        '''
        try:
            id = int(id) # convert string to int if possible
        except ValueError:
            pass # it's alphanumeric, as in '7c', leave it as a string
        '''
        return id
        
    id = property(get_id)

    def tree(self):
        """Print tree hierarchy"""
        print(self.treebuf.getvalue(), end='')

    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        if self.animal != None:
            self.animal.writetree(string)

    def load(self):
        treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)
        # collect recording names: 1st char of each name must be a digit, that's all:
        rnames = [ name for name in os.listdir(self.path)
                   if os.path.isdir(os.path.join(self.path, name))
                   and name[0].isdigit() ]
        rnames.sort() # alphabetical order
        dt = 0 # calculate total track duration by summing durations of all recordings
        for rname in rnames:
            path = os.path.join(self.path, rname)
            recording = Recording(path, track=self)
            recording.load()
            self.r[recording.id] = recording
            self.__setattr__('r' + str(recording.id), recording) # add shortcut attrib
            dt += recording.dt
        self.rnames = rnames # easy way to print out all recording names
        self.dt = dt
        self.dtsec = self.dt / 1e6
        self.dtmin = self.dtsec / 60
        self.dthour = self.dtmin / 60

        if len(rnames) == 0:
            return # no recordings in this track, nothing else to do

        if recording.sort == None or type(recording.sort.header) == core.SPKHeader:
            return # skip everything else due to impoverished .spk files

        # create a TrackSort with TrackNeurons:
        self.sort = TrackSort(self)
        self.sort.load()
        # load RF type for each cell, should be one big dict indexed by nid:
        rftypefname = os.path.join(self.path, self.absname + '.rftype')
        try:
            with open(rftypefname, 'r') as f:
                rftypestr = f.read()
            rftypes = eval(rftypestr)
            for nid, rftype in rftypes.items():
                assert rftype in ['simple', 'complex', 'LGN', None]
                self.alln[nid].rftype = rftype
        except IOError: # no absname.rftype file denoting RF type of each cell
            pass
        # load spike type for each cell, should be one big dict indexed by nid:
        spiketypefname = os.path.join(self.path, self.absname + '.spiketype')
        try:
            with open(spiketypefname, 'r') as f:
                spiketypestr = f.read()
            spiketypes = eval(spiketypestr)
            for nid, spiketype in spiketypes.items():
                assert spiketype in ['fast', 'slow', 'fastasym', 'slowasym']
                self.alln[nid].spiketype = spiketype
        except IOError: # no absname.spiketype file denoting RF type of each cell
            pass

        # one way of calculating self.trange:
        #tranges = np.asarray([ n.trange for n in self.alln.values() ])
        #self.trange = min(tranges[:, 0]), max(tranges[:, 1])
        # better way of calculating self.trange:
        rids = sorted(self.r.keys()) # all recording ids in self
        r0 = self.r[rids[0]]
        r1 = self.r[rids[-1]]
        assert r0.datetime == self.datetime
        self.trange = r0.td+r0.trange[0], r1.td+r1.trange[1]

        self.calc_meanrates()

        # pttype better be the same for all member recordings:
        pttype = self.r[rids[0]].pttype # init to pttype of first recording
        for rid in rids[1:]:
            r = self.r[rid]
            # if recording doesn't have a pttype, it's probably from an old .spk file,
            # so don't bother doing this test:
            if hasattr(r, 'pttype') and pttype != r.pttype:
                raise ValueError("inconsistent polytrode types %r and %r in track %s"
                                 % (pttype, r.pttype, self.id))

    def calc_meanrates(self):
        """Calculate mean firing rates of all TrackNeurons in this track"""
        TRACKNEURONPERIOD = get_ipython().user_ns['TRACKNEURONPERIOD']
        if TRACKNEURONPERIOD == 'track':
            # calc tn.meanrate using entire track duration:
            for tn in self.alln.values():
                tn.meanrate = tn.nspikes / self.dtsec
        elif TRACKNEURONPERIOD == 'trange':
            # calc tn.meanrate using duration between its first and last spike:
            for tn in self.alln.values():
                if tn.dtsec == 0:
                    tn.meanrate = 0.0
                else:
                    tn.meanrate = tn.nspikes / tn.dtsec
        else:
            raise ValueError("invalid value for TRACKNEURONPERIOD: %r" % TRACKNEURONPERIOD)

    def get_meanrates(self):
        """Return mean firing rates of all TrackNeurons in this track"""
        return np.asarray([ n.meanrate for n in self.alln.values() ])

    meanrates = property(get_meanrates)

    def meanratepdf(self, bins=None, figsize=(7.5, 6.5)):
        """Plot histogram of mean firing rates"""
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        if bins == None:
            bins = np.arange(0, 1, 0.05)
        n, mr = np.histogram(self.meanrates, bins=bins, density=False)
        binwidth = mr[1] - mr[0] # take width of first bin
        a.bar(left=mr[:-1], height=n, width=binwidth, bottom=0, color='k', ec='k')
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        a.set_xlabel('mean firing rate (Hz)')
        a.set_ylabel('neuron count')
        f.tight_layout(pad=0.3) # crop figure to contents

    # shortcuts to various attribs and properties in default sort:
    n = property(lambda self: self.sort.n)
    qn = property(lambda self: self.sort.qn)
    alln = property(lambda self: self.sort.alln)
    nspikes = property(lambda self: self.sort.nspikes)
    nneurons = property(lambda self: self.sort.nneurons)
    nqneurons = property(lambda self: self.sort.nqneurons)
    nallneurons = property(lambda self: self.sort.nallneurons)
    datetime = property(lambda self: self.sort.datetime)
    pttype = property(lambda self: self.sort.pttype)
    chanpos = property(lambda self: self.sort.chanpos)
    samplerate = property(lambda self: self.sort.samplerate)
    tres = property(lambda self: self.sort.tres)

    def get_nids(self, rids=None):
        """Return nids of active neurons common to all recordings specified in rids.
        Otherwise, return all active nids in all recordings. Active neurons in a recording
        are those with at least MINRATE mean spike rate during the recording"""
        if rids == None: # return all nids in all recordings
            rids = list(self.r.keys())
            return np.unique(np.hstack([ self.r[rid].n.keys() for rid in rids ]))
        else: # return intersection of nids of specified recordings
            nids = [ self.r[rid].n.keys() for rid in rids ]
            return core.intersect1d(nids, assume_unique=True)

    def get_allnids(self, rids=None):
        """Return nids of all neurons (active and quiet) common to all recordings
        specified in rids, ie return the intersection. If rids==None, return the union
        of all nids in the track instead"""
        if rids == None:
            rids = sorted(self.r.keys())
            allnids = np.hstack([ self.r[rid].alln.keys() for rid in rids ])
            return np.unique(allnids)
        else:
            allnids = [ self.r[rid].alln.keys() for rid in rids ]
            return core.intersect1d(allnids, assume_unique=True)

    def pospdfrec(self, neurons=None, rids=None, dim='y', nbins=10, a=None, figsize=(7.5, 6.5)):
        """Plot PDF of cell positions ('x' or 'y') along the polytrode
        for each recording to get an idea of how cells are distributed in space,
        and how that changes from one recording to the next"""
        if a == None:
            f = pl.figure(figsize=figsize)
            a = f.add_subplot(111)
        else: # add to existing axes
            a.hold(True)
            f = pl.gcf()

        if rids == None:
            rids = self.r.keys()
            rids.sort()
        for rid in rids:
            r = self.r[rid]
            r.pospdf(neurons=neurons, dim=dim, nbins=nbins, stats=False, a=a, figsize=figsize)

        return a

    def pospdf(self, neurons='all', dim='y', edges=None, nbins=10, stats=False, labels=True,
               a=None, figsize=(7.5, 6.5)):
        """Plot PDF of cell positions ('x' or 'y') along the polytrode
        to get an idea of how cells are distributed in space"""
        if neurons == 'all':
            neurons = self.alln.values()
        elif neurons == 'quiet':
            neurons = self.qn.values()
        elif neurons == 'active':
            neurons = self.n.values()
        dimi = {'x':0, 'y':1}[dim]
        p = [ n.pos[dimi] for n in neurons ] # all position values
        if edges != None:
            nbins = len(edges) - 1
            bins = edges # assume it includes rightmost bin edge
        else:
            nbins = max(nbins, 2*intround(np.sqrt(self.nneurons)))
            bins = nbins
        n, p = np.histogram(p, bins=bins) # p includes rightmost bin edge
        binwidth = p[1] - p[0] # take width of first bin in p

        if stats:
            mean = np.mean(p)
            median = np.median(p)
            argmode = n.argmax()
            mode = p[argmode] + binwidth / 2 # middle of tallest bin
            stdev = np.std(p)

        if a == None:
            f = pl.figure(figsize=figsize)
            a = f.add_subplot(111)
        else: # add to existing axes
            a.hold(True)
            f = pl.gcf()

        # use CCWHITEDICT1 for familiarity with len 10 1-based id to colour mapping
        #color = CCWHITEDICT1[int(self.id)]
        color = 'k'

        # exclude rightmost bin edge in p
        a.bar(left=p[:-1], height=n, width=binwidth, bottom=0, color=color, ec=color)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        if labels:
            a.set_title(titlestr)
            a.set_xlabel('neuron %s position (um)' % dim)
            a.set_ylabel('neuron count')

        if stats:
            # add stuff to top right of plot:
            uns = get_ipython().user_ns
            a.text(0.99, 0.99, 'mean = %.3f\n'
                               'median = %.3f\n'
                               'mode = %.3f\n'
                               'stdev = %.3f\n'
                               'minrate = %.2f Hz\n'
                               'nneurons = %d\n'
                               'dt = %d min'
                               % (mean, median, mode, stdev,
                                  uns['MINRATE'], self.nneurons, intround(self.dtmin)),
                               transform = a.transAxes,
                               horizontalalignment='right',
                               verticalalignment='top')
        f.tight_layout(pad=0.3) # crop figure to contents
        f.canvas.draw() # this is needed if a != None when passed as arg
        return a

    def templates(self, chans='max', cindex='nidi'):
        """Plot cell templates in their polytrode layout. chans can be 'max', 'nneigh', 'all'.
        cindex can be 'nidi' or 'nid', but best to colour cells by nidi to maximize
        alternation."""
        from colour import CCBLACKDICT0, CCBLACKDICT1
        from matplotlib.collections import LineCollection

        HUMPERINCH = 80 # for setting figure size in inches
        VUMPERINCH = 160 # for setting figure size in inches
        USPERUM = 15
        UVPERUM = 3
        HBORDERUS = 50 # us, horizontal border around chans
        VBORDERUV = 150 # uV, vertical border around plots
        HBORDER = HBORDERUS / USPERUM
        VBORDER = VBORDERUV / UVPERUM
        BG = 'black'
        SCALE = 500, 100 # scalebar size in (us, uV)
        SCALE = SCALE[0]/USPERUM, SCALE[1]/UVPERUM # um
        SCALEXOFFSET = 2 # um
        SCALEYOFFSET = 4 # um

        if chans not in ['max', 'nneigh', 'all',]:
            raise ValueError('unknown chans arg %r' % chans)
        if cindex == 'nidi':
            ccdict = CCBLACKDICT0 # use nidi to maximize colour alternation
        elif cindex == 'nid':
            ccdict = CCBLACKDICT1 # use nid to have colours that correspond to those in spyke
        else:
            raise ValueError('unknown cindex arg %r' % cindex)

        # for mpl, convert probe chanpos to center bottom origin instead of center top,
        # i.e. invert the y values:
        chanpos = self.sort.chanpos.copy()
        maxy = chanpos[:, 1].max()
        for chan, (x, y) in enumerate(chanpos):
            chanpos[chan, 1] = maxy - y

        if chans == 'nneigh': # generate dict of nearest neighbours indexed by maxchan
            dm = core.eucd(chanpos) # distance matrix
            minspace = dm[dm!=0].min()
            rincl = minspace * 1.1 # inclusion radius
            nneighs = {}
            for maxchan, pos in enumerate(chanpos):
                d = dm[maxchan]
                nnchans = np.where(d < rincl)[0]
                nneighs[maxchan] = nnchans

        colxs = np.unique(chanpos[:, 0]) # unique column x positions, sorted
        rowys = np.unique(chanpos[:, 1]) # unique row y positions, sorted
        ncols = len(colxs)
        nrows = len(rowys)
        hspace = (colxs[-1]-colxs[0]) / (ncols-1)
        vspace = (rowys[-1]-rowys[0]) / (nrows-1)

        # setting figure size actually sets window size, including toolbar and statusbar
        figwidth = (ncols*hspace + 2*HBORDER) / HUMPERINCH # inches
        figheight = (nrows*vspace + 2*VBORDER) / VUMPERINCH # inches
        dpi = mpl.rcParams['figure.dpi']
        #figwidth = (ncols*hspace) / HUMPERINCH # inches
        #figheight = (nrows*vspace) / VUMPERINCH # inches
        figwidth = intround(figwidth * dpi) / dpi # inches, rounded to nearest pixel
        figheight = intround(figheight * dpi) / dpi # inches, rounded to nearest pixel
        figsize = figwidth, figheight
        f = pl.figure(figsize=figsize, facecolor=BG, edgecolor=BG)
        a = f.add_subplot(111)

        # plot chan lines? maybe just the vertical lines?
        #for pos in chanpos:

        tres = self.sort.tres # time resolution, in us
        nts = np.unique([ neuron.nt for neuron in self.alln.values() ])
        if len(nts) != 1:
            raise RuntimeError("Not all neuron templates have the same number of timepoints. "
                               "That's probably bad.")
        nt = nts[0]
        ts = np.arange(0, neuron.nt*tres, tres) # time values in us

        nids = sorted(self.alln)
        for nidi, nid in enumerate(nids):
            colour = ccdict[eval(cindex)]
            neuron = self.alln[nid]
            # ncs (neuron channels) should be 0-based channel IDs:
            if chans == 'max':
                ncs = [neuron.maxchan]
            elif chans == 'nneigh':
                ncs = nneighs[neuron.maxchan]
            elif chans == 'all':
                ncs = neuron.chans
            # exclude channels of data within neigh that are missing from wavedata
            ncs = [ nc for nc in ncs if nc in neuron.chans ]
            # indices into neuron.chans, use to index into wavedata:
            ncis = np.hstack([ np.where(neuron.chans == nc)[0] for nc in ncs ])
            #import pdb; pdb.set_trace()
            wavedata = neuron.wavedata[ncis]
            # much less efficient, but much simpler than spyke code:
            for c, wd in zip(ncs, wavedata):
                x = chanpos[c, 0] + ts / USPERUM # um
                y = chanpos[c, 1] + wd / UVPERUM # um
                a.plot(x, y, ls='-', marker=None, lw=1, c=colour)

        a.set_axis_bgcolor(BG)
        a.set_xlabel('')
        a.set_ylabel('')
        a.xaxis.set_ticks([])
        a.yaxis.set_ticks([]) # if displayed, y ticks would be distance from bottom chan

        a.set_xlim(colxs[0]-HBORDER, colxs[-1]+nt*tres/USPERUM+HBORDER) # um
        a.set_ylim(rowys[0]-VBORDER, rowys[-1]+VBORDER) # um

        # add scale bars:
        r, b = a.get_xlim()[1]-SCALEXOFFSET, a.get_ylim()[0]+SCALEYOFFSET # um
        hbar = (r-SCALE[0], b), (r, b) # um
        vbar = (r, b), (r, b+SCALE[1]) # um
        scale = LineCollection([hbar, vbar], lw=1, colors='white', zorder=-1,
                               antialiased=True, visible=True)
        a.add_collection(scale) # add to axes' pool of LCs

        f.tight_layout(pad=0)
        #f.canvas.toolbar.hide()
        #f.canvas.window().statusBar().hide()
        f.canvas.set_window_title(lastcmd())


    def npos(self, inchespermicron=0.007, colour='active', legend=False):
        """Plot (x, y) cell positions over top of polytrode channel positions, to get an idea
        of how cells are distributed in space. Color cells by 'active', 'rftype' or
        'spiketype'."""
        uns = get_ipython().user_ns
        npos = np.asarray([ neuron.pos for neuron in self.alln.values() ])
        chanpos = self.chanpos
        chanxs, chanys = chanpos[:, 0], chanpos[:, 1]
        uchanxs = np.unique(chanxs)
        xspace = np.diff(uchanxs).max() # max spacing of consecutive unique x chan positions
        hsw = uns['PTSHANKWIDTHS'][self.pttype] / 2 # half shank width
        xs = np.hstack((npos[:, 0], chanxs, [-hsw, hsw]))
        ys = np.hstack((npos[:, 1], chanys))
        ymin = min(min(ys), 0)
        xlim = min(xs.min(), uchanxs[0]-xspace/2), max(xs.max(), uchanxs[-1]+xspace/2)
        ylim = ys.max()+xspace, ymin # inverted y axis
        
        figwidth = inchespermicron * np.ptp(xlim) * 2 + 3*legend # make space for y axis labels
        figheight = inchespermicron * np.ptp(ylim)
        f = pl.figure(figsize=(figwidth, figheight))
        a = f.add_subplot(111, aspect='equal')
        a.set_frame_on(False)
        # plot rectangle representing shank width and length, excluding the tip:
        sl = ylim[0]
        # starting from bottom left, going clockwise:
        shankxs = -hsw, -hsw, hsw, hsw
        shankys = sl, ymin, ymin, sl
        a.fill(shankxs, shankys, color='lightgrey', ec='none')
        # plot electrode sites:
        a.plot(chanpos[:, 0], chanpos[:, 1], 'k.', ms=5)
        if colour == 'active':
            # plot active and quiet cell positions in red and blue, respectively:
            anpos = np.asarray([ neuron.pos for neuron in self.n.values() ])
            qnpos = np.asarray([ neuron.pos for neuron in self.qn.values() ])
            na = len(anpos)
            nq = len(qnpos)
            # layer in inverse order of importance:
            if na: a.plot(qnpos[:, 0], qnpos[:, 1], 'b.', ms=10, alpha=0.6, label='quiet')
            if nq: a.plot(anpos[:, 0], anpos[:, 1], 'r.', ms=10, alpha=0.6, label='active')
        elif colour == 'rftype':
            # plot simple, complex, LGN afferent and None in red, blue, green and grey:
            spos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                if neuron.rftype == 'simple' ])
            cpos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                if neuron.rftype == 'complex' ])
            Lpos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                if neuron.rftype == 'LGN' ])
            Npos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                if neuron.rftype == None ])
            ns = len(spos)
            nc = len(cpos)
            nL = len(Lpos)
            nN = len(Npos)
            # layer in inverse order of importance:
            if nN: a.plot(Npos[:, 0], Npos[:, 1], 'e.', ms=10, alpha=0.6, label='unknown')
            if nL: a.plot(Lpos[:, 0], Lpos[:, 1], 'g.', ms=10, alpha=0.6, label='LGN afferent')
            if nc: a.plot(cpos[:, 0], cpos[:, 1], 'b.', ms=10, alpha=0.6, label='complex')
            if ns: a.plot(spos[:, 0], spos[:, 1], 'r.', ms=10, alpha=0.6, label='simple')
        elif colour == 'spiketype':
            # plot fast, slow, fastasym and slowasym in red, blue, green and grey:
            fpos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                if neuron.spiketype == 'fast' ])
            spos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                if neuron.spiketype == 'slow' ])
            fapos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                 if neuron.spiketype == 'fastasym' ])
            sapos = np.asarray([ neuron.pos for neuron in self.alln.values()
                                 if neuron.spiketype == 'slowasym' ])
            nf = len(fpos)
            ns = len(spos)
            nfa = len(fapos)
            nsa = len(sapos)
            # layer in inverse order of frequency:
            if nf: a.plot(fpos[:, 0], fpos[:, 1], 'r.', ms=10, alpha=0.6, label='fast')
            if ns: a.plot(spos[:, 0], spos[:, 1], 'b.', ms=10, alpha=0.6, label='slow')
            if nfa: a.plot(fapos[:, 0], fapos[:, 1], 'g.', ms=10, alpha=0.6,
                           label='fast asymmetric')
            if nsa: a.plot(sapos[:, 0], sapos[:, 1], 'e.', ms=10, alpha=0.6,
                           label='slow asymmetric')
        else:
            raise RuntimeError("unknown colour kwarg %r" % colour)
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xticks(uchanxs)
        a.set_yticks(np.arange(0, ylim[0], 200))
        #a.xaxis.set_ticks_position('bottom')
        #a.yaxis.set_ticks_position('left')
        # put legend to right of the axes:
        if legend:
            a.legend(loc='center left', bbox_to_anchor=(1.2, 0.5), frameon=False)
        bbox = a.get_position()
        wh = bbox.width / bbox.height # w:h ratio of axes, includes all ticks and labels?
        w, h = gcfm().canvas.get_width_height()
        gcfm().resize(w*wh, h)
        titlestr = lastcmd()
        gcfm().set_window_title(titlestr)
        a.set_title(self.absname)
        #a.set_xlabel('$\mu$m')
        #a.set_ylabel('$\mu$m')
        #f.tight_layout(pad=0.3) # resizes contents to figure (not crop figure to contents!)

    def scstim(self, method='mean', width=None, tres=None, figsize=(7.5, 6.5)):
        """Scatter plot some summary statistic of spike correlations of each recording,
        classified by the stimulus group each recording falls into. width and tres dictate
        tranges to split recordings up into, if any"""

        ## TODO: for each pair of recordings, find common subset of active neurons and calculate
        ## pairwise corrs for each recording in that pair using just those neurons

        ## TODO: maybe limit to visually responsive cells

        uns = get_ipython().user_ns
        if width == None:
            width = uns['SCWIDTH']
        if tres == None:
            tres = width
        blankmseqrids = uns['BSRIDS'][self.absname] + uns['MSRIDS'][self.absname]
        movdriftrids = uns['NSRIDS'][self.absname] + uns['DBRIDS'][self.absname]

        blankmseqcorrs = []
        movdriftcorrs = []
        for rid in (blankmseqrids + movdriftrids):
            r = self.r[rid]
            print('%s: %s' % (r.absname, r.name))
            spikecorr = r.sc(width=width, tres=tres)
            sc = spikecorr.sct(method=method)[0]
            sc = sc[0] # pull out the spike correlation values that span all laminae
            if rid in blankmseqrids:
                blankmseqcorrs.append(sc)
            else:
                movdriftcorrs.append(sc)
        blankmseqcorrs = np.hstack(blankmseqcorrs)
        movdriftcorrs = np.hstack(movdriftcorrs)
        # repeat each element in blankmseqcorrs len(movdriftcorrs) times:
        x = np.repeat(blankmseqcorrs, len(movdriftcorrs))
        # tile movdriftcorrs len(blankmseqcorrs) times:
        y = np.tile(movdriftcorrs, len(blankmseqcorrs))

        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        lim = min([x.min(), y.min(), 0]), max([x.max(), y.max()])
        a.plot(lim, lim, c='e', ls='--', marker=None) # y=x line
        a.plot(x, y, 'k.')
        #a.set_xlim(lim)
        #a.set_ylim(lim)
        a.set_xlabel('%s spike correlations: blankscreen and mseq' % method)
        a.set_ylabel('%s spike correlations: movie and drift bar' % method)
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        f.tight_layout(pad=0.3) # crop figure to contents
        f.show()

    def scsistim(self, method='mean', width=None, tres=None, timeaverage=False,
                 plottime=False, s=5, figsize=(7.5, 6.5)):
        """Scatter plot some summary statistic of spike correlations of each recording vs
        LFP synchrony index SI. Colour each point according to stimulus type. width and tres
        (sec) dictate tranges to split recordings up into. timeaverage averages across time
        values of both sc and si for each recording. s is point size"""
        ## TODO: maybe limit to visually responsive cells
        ## TODO: add linear regression of si vs log(sc)

        uns = get_ipython().user_ns
        if width == None:
            width = uns['LFPSIWIDTH']
        if tres == None:
            tres = width
        bsrids = uns['BSRIDS'][self.absname]
        msrids = uns['MSRIDS'][self.absname]
        mvrids = uns['NSRIDS'][self.absname]
        dbrids = uns['DBRIDS'][self.absname]
        rids = sorted(bsrids + msrids + mvrids + dbrids) # do everything in rid order
        print('blankscreen: %r' % [self.r[rid].name for rid in bsrids])
        print('mseq: %r' % [self.r[rid].name for rid in msrids])
        print('movie: %r' % [self.r[rid].name for rid in mvrids])
        print('driftbar: %r' % [self.r[rid].name for rid in dbrids])
        isect = core.intersect1d([msrids, bsrids, mvrids, dbrids])
        if len(isect) != 0:
            raise RuntimeError("some rids were classified into more than one type: %r" % isect)

        scs, sis, c = [], [], []
        for rid in rids:
            r = self.r[rid]
            print('%s: %s' % (r.absname, r.name))
            spikecorr = r.sc(width=width, tres=tres)
            sc, si = spikecorr.si(method=method, plot=False) # calls sc.sct() and sc.si()
            sc = sc[0] # pull out the spike correlation values that span all laminae
            if timeaverage:
                # average across all time values of sc and si to get a single coordinate
                # per recording
                sc = sc.mean()
                si = si.mean()
            scs.append(sc)
            sis.append(si)
            if rid in bsrids: color = 'e'
            elif rid in msrids: color = 'k'
            elif rid in mvrids: color = 'r'
            elif rid in dbrids: color = 'b'
            else: raise ValueError("unclassified recording: %r" % r.name)
            c.append(np.tile(color, len(sc)))
        scs = np.hstack(scs)
        sis = np.hstack(sis)
        c = np.hstack(c)
        
        f = pl.figure(figsize=figsize)
        a = f.add_subplot(111)
        if plottime: # underplot lines connecting points adjacent in time
            a.plot(scs, sis, 'e--')
        a.scatter(scs, sis, c=c, edgecolors='none', s=s)
        a.set_ylim(0, 1)
        a.set_xlabel('%s spike correlations' % method)
        a.set_ylabel('synchrony index')
        titlestr = lastcmd()
        gcfm().window.setWindowTitle(titlestr)
        a.set_title(titlestr)
        # make proxy line artists for legend:
        bs = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='e', mec='e')
        ms = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='k', mec='k')
        mv = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='r', mec='r')
        db = mpl.lines.Line2D([1], [1], color='white', marker='o', mfc='b', mec='b')
        # add legend:
        a.legend([bs, ms, mv, db],
                 ['blank screen', 'mseq', 'movie', 'drift bar'],
                 numpoints=1, loc='lower right',
                 handlelength=1, handletextpad=0.5, labelspacing=0.1)
        f.tight_layout(pad=0.3) # crop figure to contents
        return scs, sis, c
