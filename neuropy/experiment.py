"""Defines the Experiment class"""

from __future__ import division

import os
import StringIO

import numpy as np
import matplotlib as mpl

from core import getargstr, TAB, warn, rstrip, dictattr, intround, toiter
from core import _movies, MOVIEPATH, MSEQ16, MSEQ32, joinpath, lastcmd
from core import PopulationRaster, Codes, CodeCorrPDF, RevCorrWindow
import neuron

from dimstimskeletal import deg2pix, InternalParams, StaticParams, DynamicParams
from dimstimskeletal import Variable, Variables, Runs, BlankSweeps, Dimension, SweepTable
from dimstimskeletal import Movie, Grating, Bar, SparseNoise, BlankScreen


class BaseExperiment(object):
    """An experiment corresponds to a single contiguous stimulus session.
    It contains information about the stimulus during that session, including
    the DIN values and the text header. For data generated with dimstim >= 0.16,
    it includes the entire dimstim.Experiment object as an attribute (.e).
    A neuropy.Experiment is basically a container for a dimstim.Experiment"""
    def __init__(self, path, id=None, recording=None):
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        self.path = path
        self.id = id
        self.r = recording

    def get_name(self):
        fname = os.path.split(self.path)[-1]
        return rstrip(fname, '.din')

    name = property(get_name)

    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),

    def writetree(self, string):
        """Write to self's tree buffer and to parent's too"""
        self.treebuf.write(string)
        self.r.writetree(string)

    def load(self):
        f = open(self.path, 'rb')
        self.din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to nrows x 2 cols
        f.close()
        try:
            txthdrpath = rstrip(self.path, '.din') + '.textheader'
            f = open(txthdrpath, 'rU') # use universal newline support
            self.textheader = f.read() # read it all in
            f.close()
        except IOError:
            warn('Error reading: <%s>, text header not loaded' % f.name)
            self.textheader = '' # set to empty

        treestr = self.level*TAB + self.name + '/'
        # print string to tree hierarchy and screen
        self.writetree(treestr + '\n')
        print(treestr)

        if self.textheader != '':
            # comment out all lines starting with "from dimstim"
            self.textheader = self.textheader.replace('from dimstim', '#from dimstim')
            names1 = locals().copy() # namespace before execing the textheader
            exec(self.textheader)
            names2 = locals().copy() # namespace after
            # names that were added to the namespace, excluding the 'names1' name itself:
            newnames = [ n2 for n2 in names2 if n2 not in names1 and n2 != 'names1' ]
            try:
                # dimstim up to Cat 15 didn't have a version, neither did NVS display
                self.__version__ = eval('__version__')
            except NameError:
                self.__version__ = 0.0
            if self.__version__ >= 0.16: # after major refactoring of dimstim
                for newname in newnames:
                    # bind each variable in the textheader as an attrib of self
                    self.__setattr__(newname, eval(newname))
                self.sweeptable = SweepTable(experiment=self.e) # build the sweep table
                self.st = self.sweeptable.data # synonym, used a lot by experiment subclasses
                # this doesn't work for textheaders from dimstim 0.16, since xorigDeg and yorigDeg
                # were accidentally omitted from all the experiment scripts and hence the textheaders
                # too:
                '''
                self.e.xorig = deg2pix(self.e.static.xorigDeg, self.I) + self.I.SCREENWIDTH / 2
                self.e.yorig = deg2pix(self.e.static.yorigDeg, self.I) + self.I.SCREENHEIGHT / 2
                '''
                self.REFRESHTIME = intround(1 / float(self.I.REFRESHRATE) * 1000000) # us
                # prevent replication of movie frame data in memory
                if type(self.e) == Movie:
                    fname = os.path.split(self.e.static.fname)[-1] # pathless fname
                    if fname not in _movies:
                        # add movie experiment, indexed according to movie data file name,
                        # to prevent from ever loading its frames more than once
                        _movies[fname] = e
            else:
                self.oldparams = dictattr()
                for newname in newnames:
                    # bind each variable in the textheader to oldparams
                    self.oldparams[newname] = eval(newname)
                self.loadCat15exp()
        else:
            # use the time difference between the first two din instead
            self.REFRESHTIME = self.din[1, 0] - self.din[0, 0]

        # add an extra refresh time after last din, that's when screen actually turns off
        self.trange = (self.din[0, 0], self.din[-1, 0] + self.REFRESHTIME)

    def loadCat15exp(self):
        ## TODO: - fake a .e dimstim.Experiment object, to replace what used to be the
        ## .stims object for movie experiments
        '''           - self.movie = self.experiment.stims[0]
                - need to convert sweeptimeMsec to sweepSec
                   - assert len(self.experiment.stims) == 1
                   - self.movie = self.experiment.stims[0]
                   - self.movie.load() # ensure the movie's data is loaded

            if self.movie.oname == 'mseq32':
                frameis = frameis[frameis != 65535] # remove all occurences of 65535
            elif self.movie.oname == 'mseq16':
                frameis = frameis[frameis != 16383] # remove all occurences of 16383
        '''
        # Add .static and .dynamic params to fake dimstim experiment
        self.e = dictattr()
        self.e.I = dictattr() # fake InternalParams object
        self.e.static = dictattr() # fake StaticParams object
        self.e.dynamic = dictattr() # fake DynamicParams object
        # maps Cat 15 param names to dimstim 0.16 param types and names, wherever possible
        ## TODO: fill in params for experiment types other than Movie??
        _15to16 = {'EYE': ('I', 'EYE'),
                   'PIXPERCM': ('I', 'PIXPERCM'),
                   'REFRESHRATE': ('I', 'REFRESHRATE'),
                   'SCREENDISTANCECM': ('I', 'SCREENDISTANCECM'),
                   'SCREENHEIGHT': ('I', 'SCREENHEIGHT'),
                   'SCREENHEIGHTCM': ('I', 'SCREENHEIGHTCM'),
                   'SCREENWIDTH': ('I', 'SCREENWIDTH'),
                   'SCREENWIDTHCM': ('I', 'SCREENWIDTHCM'),

                   'fname': ('static', 'fname'),
                   'preexpSec': ('static', 'preexpSec'),
                   'postexpSec': ('static', 'postexpSec'),
                   'orioff': ('static', 'orioff'),
                   'regionwidthDeg': ('static', 'widthDeg'),
                   'regionheightDeg': ('static', 'heightDeg'),
                   'mask': ('static', 'mask'),
                   'diameterDeg': ('static', 'diameterDeg'),
                   'GAMMA': ('static', 'gamma'),

                   'framei': ('dynamic', 'framei'),
                   'ori': ('dynamic', 'ori'),
                   'polarity': ('dynamic', 'invert'),
                   'bgbrightness': ('dynamic', 'bgbrightness'),
                   'sweeptimeMsec': ('dynamic', 'sweepSec'),
                   'postsweepMsec': ('dynamic', 'postsweepSec'),
                   }

        # collect any Cat 15 movie attribs and add them to self.oldparams
        try:
            # can't really handle more than 1 movie, since dimstim 0.16 doesn't
            assert len(np.unique(self.oldparams.playlist)) == 1
            # bind it, movie was the only possible stim object anyway in Cat 15
            self.movie = self.oldparams.playlist[0]
            # returns dict of name:val pair attribs excluding __ and methods:
            movieparams = self.oldparams[self.movie.oname].__dict__
            self.oldparams.update(movieparams)
        except AttributeError:
            # no playlist, no movies, and therefore no movie attribs to deal with
            pass

        # convert Cat 15 params to dimstim 0.16
        for oldname, val in self.oldparams.items():
            if 'msec' in oldname.lower():
                val = val / 1000. # convert to sec
            elif oldname == 'polarity':
                val = bool(val) # convert from 0/1 to boolean
            if oldname == 'origDeg': # split old origDeg into new separate xposDeg and yposDeg
                self.e.dynamic.xposDeg = val[0]
                self.e.dynamic.yposDeg = val[1]
            else:
                try:
                    paramtype, newname = _15to16[oldname]
                    self.e[paramtype][newname] = val
                except KeyError: # oldname doesn't have a newname equivalent
                    pass

        try:
            m = self.movie
        except AttributeError:
            m = None

        if m:
            # make fake dimstim experiment a Cat15Movie object, bind all of the attribs of
            # the existing fake dimstim experiment
            old_e = self.e
            self.e = m
            for name, val in old_e.__dict__.items():
                # bind each variable in the textheader as an attrib of self
                self.e.__setattr__(name, val)
            # deal with movie filename:
            # didn't have a chance to pass this exp as the parent in the movie init,
            # so just set the attribute manually:
            m.e = self
            # if fname refers to a movie whose local name is different, rename it to match
            # the local movie name
            _old2new = {'mseq16.m': MSEQ16, 'mseq32.m': MSEQ32}
            try:
                m.fname = _old2new[m.fname]
            except KeyError:
                pass # old name not in _old2new, leave it be
            self.e.static.fname = m.fname # update fake dimstim experiment's fname too
            # extensionless fname, fname should've been defined in the textheader
            m.name = os.path.splitext(m.fname)[0]
            if m.name not in _movies:
                # and it very well may not be, cuz the textheader inits movies with no args,
                # leaving fname==None at first, which prevents it from being added to
                # _movies
                _movies[m.name] = m # add m to _movies dictattr
            # Search self.e.moviepath string (from textheader) for 'Movies' word. Everything
            # after that is the relative path to your base movies folder. Eg, if
            # self.e.moviepath = 'C:\\Desktop\\Movies\\reliability\\e\\single\\', then set
            # self.e.relpath = '\\reliability\\e\\single\\'
            spath = self.oldparams.moviepath.split('\\') # Cat15 has purely windows seperators
            matchi = spath.index('Movies')
            relpath = joinpath(spath[matchi+1 ::])
            path = os.path.join(MOVIEPATH, relpath)
            m.fname = os.path.join(path, m.fname)
            self.e.static.fname = m.fname # update

            # Generate the sweeptable:
            # dict of lists, ie sweeptable={'ori':[0,45,90], 'sfreq':[1,1,1]}, so you index
            # into it with self.sweeptable['var'][sweepi]
            #self.sweeptable = {[]}
            #vars = self.sweeptable.keys()
            # need to check if varlist exists, if so use it (we're dealing with Cat 15),
            # if not, use revamped dimstim.SweepTable class
            varvals = {} # init a dictionary that will contain variable values
            for var in m.varlist:
                # generate dict with var:val entries, to pass to buildSweepTable
                varvals[var] = eval('m.' + var)
            # pass varlist by reference, dim indices end up being modified:
            m.sweepTable = self.buildCat15SweepTable(m.varlist, varvals, m.nruns,
                                                     m.shuffleRuns, m.blankSweep,
                                                     m.shuffleBlankSweeps,
                                                     makeSweepTableText=0)[0]
        else: # this is a simple stim (not object oriented movie)
            varvals = {} # init a dictionary that will contain variable values
            for var in self.oldparams.varlist:
                # generate dict with var:val entries, to pass to buildSweepTable
                varvals[var] = eval('self.oldparams.' + var)
            # pass varlist by reference, dim indices end up being modified:
            self.sweepTable = self.buildCat15SweepTable(self.oldparams.varlist, varvals,
                                                        self.oldparams.nruns,
                                                        self.oldparams.shuffleRuns,
                                                        self.oldparams.blankSweep,
                                                        self.oldparams.shuffleBlankSweeps,
                                                        makeSweepTableText=0)[0]
        try:
            self.REFRESHTIME = intround(1 / float(self.oldparams.REFRESHRATE) * 1000000) # us
        except AttributeError:
            pass

    def buildCat15SweepTable(self, varlist, varvals, nruns=1, shuffleRuns=0, blankSweep=(0,0),
                             shuffleBlankSweeps=0, makeSweepTableText=0):
        """Deprecated: kept only for backward compatibility with Cat 15
        Builds a sweep table
        Returns: 'sweeptable': a dictionary where each entry (key) is a variable name, followed
                               by its values, one value per unique sweep (ie permutation of dims)
                 'dimlist': a dictionary of renumbered dimensions
                 'sweeplist': stores the # of times and in what order we'll be stepping through
                              sweeptable during the experiment
                 'sweeptabletext': a formatted tab-separated string printout of 'sweeptable'
        Passed:  'varlist': a dictionary (see below)
                 'varvals': a dictionary with var:val(s) key:value(s) entries
                 'nruns': number of times to run the whole combination of parameters
                 'blankSweep': (n,duration) tuple, do a blank sweep every n'th sweep for
                               duration in seconds, n must be in: {0=off ,2,3,4....}
                 'shuffleBlankSweeps': shuffle the position of the blank sweeps? (0=no, 1=yes)
                 'makeSweepTableText': generate text printout of sweeptable? (0=no, 1=yes)

        Each entry in 'varlist' represents one variable. Each entry is itself made up of a 2-entry dictionary with 'shuffle' and 'dim' fields. The 'shuffle' field is a flag (0=leave ordered, 1=shuffle, 2=randomize). shuffle==1 shuffles the variable's dimension during the experiment, shuffle==2 randomly samples it instead. Variable position in the varlist doesn't matter. However, the variable's associated 'dim' value relative to all the other 'dim' values in the variable list determines its order in the nested for loops that generate the combinations of values for each sweep: the variable with the lowest 'dim' value becomes the outermost for loop and changes least often; the variable with the highest 'dim' value becomes the innermost for loop and changes on every sweep. 'dim' must be an integer (+ve, -ve, or 0). Variables with the same 'dim' value are part of the same dimension, are shuffled together, and must therefore be assigned the same number of values and the same shuffle flag"""

        # Do some error checking
        if type(nruns) != int or nruns < 0:
            raise ValueError, 'nruns must be a non-negative integer'
        if type(shuffleRuns) != int or shuffleRuns not in (0, 1):
            raise ValueError, 'shuffleRuns must be 0 or 1'
        if type(blankSweep) != tuple or type(blankSweep[0]) != int or blankSweep[0] == 1 or blankSweep[0] < 0:
            raise ValueError, 'blankSweep must be a tuple, and blankSweep[0] must be a non-negative integer and cannot equal 1'
        if type(shuffleBlankSweeps) != int or shuffleBlankSweeps not in (0, 1):
            raise ValueError, 'shuffleBlankSweeps must be 0 or 1'
        # check that shuffle flag is valid for each variable
        for var in varlist:
            if type(varlist[var]['shuffle']) != int or varlist[var]['shuffle'] not in (0, 1, 2):
                raise ValueError, 'Variable %s has a shuffle value outside of (0, 1, 2)' % var

        for (key, val) in varvals.items():
            if val == None:
                raise ValueError, 'Variable %s was passed with value None. Can\'t generate a sweeptable with None values' % key
            if type(val) != list:
                val = [val] # turn any single value non-list vals into single value list vals
            exec(key + '=' + str(val)) # locally allocate the vals for all vars in varlist - bad!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nvars = len(varlist)

        # Print out values for vars in varlist, as well as the shuffle and dim fields for each var
        #print 'varlist with vals:'
        #for var in varlist:
        #    print var + ' =', eval(var), ',', varlist[var]
        #print

        # First find total number of distinct dimensions in varlist. The dimension numbers entered in the 'dim' field in varlist don't have to be consecutive. Eg, they could be: 0,5,5,5,2,666,-74,...
        dimindexlist = []
        for key in varlist:
            dim = varlist[key]['dim']
            if dim not in dimindexlist:
                dimindexlist.append(dim)
        ndims = len(dimindexlist)
        dimindexlist.sort()

        # Now, renumber all the dimensions 0,1,2,3... according to their numerical order. Also, copy info stored in the varlist dictionary into simpler lists that can then be easily used to generate the sweeptable
        dimlist = [] # Every n'th entry in dimlist is a dictionary that corresponds to dimension n. Keys in each dictionary include 'vars' (list of variable names that belong to the dim), 'shuffle' (shuffle flag), and 'dim' (n)
        dimlengths = [] # stores the length of each dimension (the # of conds of each of its vars, all vars in a dim should have the same # of conds)
        dimshuffles = [] # stores the shuffle flag for each dimensionas
        dimvars = [] # stores the variable names that belong to each dimension
        for dim in range(ndims): # for all dimensions
            varlengths=[] # this will temporarily store the lengths of each var in this dim
            varshuffles=[] # this will temporarily store the shuffle flags of each var in this dim
            dimvars.append([]) # init dimvars for this dim, ready to add the variable names that belong to this dim
            dimi = dimindexlist[dim] # dimi is the original dimension # assigned in varlist, dim is the corresponding 0 based dimension #
            # go through varlist, grabbing dimi number from each member in varlist
            for (var, vardict) in varlist.items(): # varlist is a dictionary of dictionaries, var is the variable name
                if vardict['dim']==dimi:
                    vardict['dim']=dim # overwrite dimi with dim, so dim nums are 0 based and consecutive, note that this overwrites the dimi value in varlist (OK) cuz it was passed by reference and we're stepping through it by reference!
                    if type(eval(var)) == list:
                        varlengths.append(len(eval(var))) # store the variable's length
                    else:
                        varlengths.append(1) # store length of single numbers (non-lists) as 1
                    varshuffles.append(vardict['shuffle']) # store the shuffle flag for this var in this dim
                    dimvars[dim].append(var) # store the variable names that are in this dimension
            # check if each var in this dim has the same length as the first var, and if they don't all have the same length, raise error
            for varlength in varlengths:
                if varlength != varlengths[0]:
                    raise ValueError, 'Dimension %d contains variables with an unequal number of conditions' % dimi
            dimlengths.append(varlength) # store the length of this dimension
            # check if each var in this dim has the same shuffle flag as the first var
            for varshuffle in varshuffles:
                if varshuffle != varshuffles[0]:
                    raise ValueError, 'Dimension %d contains variables with different shuffle flags' % dimi
            dimshuffles.append(varshuffle) # store the shuffle flag for this dimension

            dimlist.append({}) # add a dictionary for this dim in dimlist
            dimlist[dim]={'vars':dimvars[dim], 'shuffle':dimshuffles[dim], 'dim':dim} # write the entry for this dim into dimlist

        # Print dimlist
        #print 'Dimension list:'
        #for dim in range(ndims):
        #    print dimlist[dim]
        #print

        # Build the ordered (unshuffled/unrandomized) indextable
        vali = [None]*ndims # stores the index we're on in each dimension
        indextable = [] # stores the ordered table of indices for all dimensions over all sweeps
        # generate code with the right number of nested for loops
        code = ''
        tabs = ''
        for dim in range(ndims): # here come the nested for loops...

            code += tabs+'for vali['+str(dim)+'] in xrange(dimlengths['+str(dim)+']):\n' # val index for this dim, used for all vars in this dim
            tabs +='\t'

        code += tabs+'indextable.append(vali[:])\n' # here's the innermost part of the nested for loops, val[:] returns a copy (important!)
        #print code
        #print
        exec(code) # run the generated code, this builds the ordered indextable with all the permutations

        # example of what the generated code looks like for 3 dimensions:
        #for vali[0] in range(dimlengths[0]):
        #    for vali[1] in range(dimlengths[1]):
        #        for vali[2] in range(dimlengths[2]):
        #            indextable.append(vali[:])

        nsweeps = len(indextable)

        # Now use indextable to build the sweeptable
        sweeptable = {}
        for var in varlist:
            sweeptable[var] = [None]*nsweeps # init with all the var names as the key names
        for dimi in range(ndims):
            for var in dimvars[dimi]:
                evalvar = eval(var)
                for sweepi in range(nsweeps):
                    sweeptable[var][sweepi] = evalvar[ indextable[sweepi][dimi] ] # eg sweeptable['ori'][sweepi]=ori[ indextable[sweepi][dimi] ]

        #sweeptable['sweepi'] = range(nsweeps) # add ordered sweep indices to a new key in sweeptable called 'sweepi'

        # The 'sweeplist' will store the # of times and in what order we'll be stepping through the sweeptable during the experiment
        if shuffleRuns:
            nshuffles = nruns # re-shuffle sweeplist for each run
            ncopies = 1 # make sweeplist ncopies of itself long after shuffling
        else:
            nshuffles = 1 # shuffle sweeplist only once, use same sweeplist for each run
            ncopies = nruns # make sweeplist ncopies of itself long after shuffling

        #t1 = time.clock()
        sweeplistlist = [None]*nsweeps*nshuffles # init list of concatenated sweeplists for speed
        for shufflei in range(nshuffles):
            sweeplist = range(nsweeps) # init with an ordered sweeplist for each shuffle run

            # Do the appropriate shuffling by overwriting the appropriate entries in sweeplist in the appropriate way

            # check if all dims are set to shuffle/randomize, if so, do it the fast way
            allshuffled = 1
            allrandomized = 1
            for dimi in range(ndims):
                allshuffled *= (dimshuffles[dimi]==1)
                allrandomized *= (dimshuffles[dimi]==2)
            if allshuffled:
                sweeplist = shuffle(sweeplist)
            elif allrandomized:
                sweeplist = randomize(sweeplist)
            else: # shuffle/randomize each dim individually (slower)
                for dimi in xrange(ndims):
                    if dimshuffles[dimi] in (1, 2): # if flag is set to shuffle or randomize
                        col = [None]*nsweeps # init for speed
                        for sweepi in sweeplist: # extract column from indextable for dimi
                            col[sweepi] = indextable[sweepi][dimi]
                        #col = [ indextable[sweepi][dimi] for sweepi in sweeplist ]
                        # use numpy's argsort instead: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        colvi = [ (v, i) for i, v in enumerate(col) ] # load up col values and indices into a list of tuples
                        colvi.sort() # sort the pair according to v (first entry in (v,i) tuple), the indices will go along for the ride
                        sortindices = [ i for v, i in colvi ] # unload the sorted indices, these are the indices into sweeplist that would get you a sorted sweeplist
                        sortedsweeplist = [None]*nsweeps # init for speed
                        for i, si in enumerate(sortindices): # create sweeplist sorted in order of dimi
                            sortedsweeplist[i] = sweeplist[si]
                        #sortedsweeplist = [ sweeplist[si] for si in sortindices ]

                        offset = 1 # offset is the product of the lengths of all dims other than dimi
                        for dimj in xrange(ndims):
                            if dimj != dimi: # for all dims other than dimi
                                offset *= dimlengths[dimj]

                        ldimi = dimlengths[dimi] # length of dimi, including pre-shuffle repeats
                        if nsweeps % ldimi != 0:
                            raise ValueError, 'Somehow, nsweeps isn\'t an integer multiple of length (including pre-shuffle repeats) of dim %d' % dimi
                        ncollections = nsweeps // ldimi # can safely divide now

                        for colli in xrange(ncollections):
                            # collis is a collection of indices to shuffle over, made up of every offset'th index, starting from colli
                            collis = [None]*ldimi # init for speed
                            for i, j in enumerate(xrange(colli, offset*ldimi, offset)):
                                collis[i] = j
                            #collis = [ j for j in xrange(colli,offset*ldimi,offset) ]
                            if dimshuffles[dimi] == 1: # shuffle this dim
                                shuffcollis = shuffle(collis)
                            elif dimshuffles[dimi] == 2: # randomize this dim
                                shuffcollis = randomize(collis)
                            for i in xrange(len(collis)):
                                sweeplist[sortindices[collis[i]]] = sortedsweeplist[shuffcollis[i]] # update sweeplist appropriately, this is the trickiest bit
            for i, j in enumerate(range(nsweeps*shufflei, nsweeps*(shufflei+1))):
                sweeplistlist[j] = sweeplist[i] # enters a copy of sweeplist for this shuffle into the concatenated list of sweeplists

        sweeplist = sweeplistlist # rename it now that we're done the shuffle loop
        #t2 = time.clock()
        #print 'shuffling took', t2-t1

        sweeplist *= ncopies # does nothing if shuffleRuns==1 (ncopies is 1), make nruns copies of sweeplist if shuffleRuns==0 (ncopies is nruns)

        nuniquesweeps = nsweeps
        nsweeps = len(sweeplist) # now nsweeps includes repeats (if any)

        # If blank sweeps are requested, then build the stimulus state list (stimOn[sweepi]: 0=stimulus off, 1=stimulus on) and use it to insert blank sweeps (None values) into sweeplist
        if blankSweep[0]:
            # increment nsweeps to have correct additional amount to account for added blank sweeps
            addedsweeps = 0 # number of sweeps added to nsweeps
            addsweeps = nsweeps // blankSweep[0] # both nsweeps and blankSweep[0] are integers
            while addsweeps > 0:
                nsweeps += addsweeps
                addedsweeps += addsweeps
                addsweeps = nsweeps // blankSweep[0] - addedsweeps
            stimOn = [1]*nsweeps # init stimOn to all 1s (stim is on for all sweeps)
            for sweepi in xrange(nsweeps):
                if (sweepi+1) % blankSweep[0] == 0: # if 1 based sweepi is multiple of blankSweep[0]
                    stimOn[sweepi] = 0 # stimulus will be off on sweepi
            if shuffleBlankSweeps == 1:
                stimOn = shuffle(stimOn) # shuffle the stimulus state list
            #print 'stimOn is', stimOn
            # insert blank sweeps into sweeplist, according to stimulus state list
            for sweepi in xrange(nsweeps): # for all sweeps (values) in stimulus table
                if stimOn[sweepi] == 0: # if we're on a blank sweep
                    insertioni = min(sweepi, len(sweeplist)-1) # insert at the lesser of sweepi or what's currently the last index in sweeplist. This prevents the case where many of the zeros in a shuffled stimOn are clustered at the end, and you end up trying to insert a repeat at an index value that's out of range of sweeplist
                    sweeplist.insert(insertioni, None) # insert a 'None' value into sweeplist at insertioni, this indicates a blank sweep
        #print 'sweeplist is', sweeplist

        # Print sweeptable as formatted text to string
        if makeSweepTableText:
            f = cStringIO.StringIO() # create a string file-like object, implemented in C, fast
            f.write('sweepi\t') # 'sweepi' column header
            for dim in xrange(ndims): # print out column headers in dim order
                f.write( str(dimvars[dim]) + '\t' )
            f.seek(-1, 2) # move the file pointer back one relative to end
            f.truncate() # get rid of the last tab
            f.write('\n')
            for sweepi in xrange(nuniquesweeps):
                f.write('%.6s\t' % sweepi) # write the sweep index in the first column
                for dimi in xrange(ndims): # print out columns in dim order
                    for var in dimvars[dimi]:
                        f.write('%.6s\t' % sweeptable[var][sweepi])
                f.seek(-1, 2) # move the file pointer back one relative to end
                f.truncate() # get rid of the last tab
                f.write('\n')
            f.seek(-1, 2) # move the file pointer back one relative to end
            f.truncate() # get rid of the last newline
            sweeptabletext = f.getvalue()
        else:
            sweeptabletext = None

        return sweeptable, dimlist, sweeplist, sweeptabletext


class ExperimentCode(BaseExperiment):
    """Mix-in class that defines the spike code related experiment methods"""
    def code(self, neuron=None, **kwargs):
        """Returns a Neuron.Code object, constraining it to the time range of this experiment. Takes either a Neuron object or just a Neuron id"""
        try:
            return neuron.code(tranges=[self.trange], **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].code(tranges=[self.trange], **kwargs) # neuron is probably a Neuron id
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nNeuron.code: '+getargstr(neuron.Neuron.code)
    code.__doc__ += '\nbinary: '+getargstr(neuron.BinaryCode.__init__)

    def codes(self, neurons=None, **kwargs):
        """Returns a 2D array where each row is a neuron code constrained to the time range of this experiment"""
        if neurons == None:
            neurons = self.r.n
        codeso = self.r.codes(neurons=neurons, experiments=[self], **kwargs)
        codeso.calc()
        return codeso
    codes.__doc__ += '\n\n**kwargs:'
    codes.__doc__ += '\nCodes: '+getargstr(Codes.__init__)
    codes.__doc__ += '\nNeuron.code: '+getargstr(neuron.Neuron.code)
    codes.__doc__ += '\nbinary: '+getargstr(neuron.BinaryCode.__init__)

    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code objects"""
        code1 = self.code(neuron1, **kwargs)
        code2 = self.code(neuron2, **kwargs)
        return corrcoef(code1.c, code2.c)
    codecorr.__doc__ += '\n\n**kwargs:'
    codecorr.__doc__ += '\nNeuron.code: '+getargstr(neuron.Neuron.code)
    codecorr.__doc__ += '\nbinary: '+getargstr(neuron.BinaryCode.__init__)

    def codecorrpdf(self, **kwargs):
        """Returns an existing CodeCorrPDF object, or creates a new one if necessary"""
        try:
            self._codecorrpdfs
        except AttributeError: # doesn't exist yet
            self._codecorrpdfs = [] # create a list that'll hold CodeCorrPDF objects
        cco = self.r.codecorrpdf(experiments={self.id: self}, **kwargs) # init a new one
        for ccpdf in self._codecorrpdfs:
            if cco == ccpdf: # need to define special == method for class CodeCorrPDF()
                return ccpdf # returns the first object whose attributes match what's desired. This saves on calc() time and avoids duplicates in self._codecorrpdfs
        cco.calc() # no matching object was found, calculate it
        self._codecorrpdfs.append(cco) # add it to the object list
        return cco
    codecorrpdf.__doc__ += '\n\n**kwargs:'
    codecorrpdf.__doc__ += '\nCodeCorrPDF: '+getargstr(CodeCorrPDF.__init__)
    codecorrpdf.__doc__ += '\nNeuron.code: '+getargstr(neuron.Neuron.code)
    codecorrpdf.__doc__ += '\nbinary: '+getargstr(neuron.BinaryCode.__init__)
    '''
    def netstate(self):
        """Returns a Netstate object"""
        nso = Netstate(experiment=self)
        return nso
    def codewords(self, **kwargs):
        cw = CodeWords(trange=self.trange)
        cw.calc()
        return cw
    '''

class ExperimentRate(BaseExperiment):
    """Mix-in class that defines the spike rate related experiment methods"""
    def rate(self, neuron, **kwargs):
        """Returns a Neuron.Rate object, constraining it to the time range of this experiment. Takes either a neuron object or just a neuron id"""
        try:
            return neuron.rate(trange=self.trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].rate(trange=self.trange, **kwargs) # neuron is probably a Neuron id
    rate.__doc__ += '\n\n**kwargs:'
    rate.__doc__ += neuron.Neuron._rateargs

    def ratepdf(self, neuron, **kwargs):
        """Returns a Neuron.RatePDF object, constraining it to the time range of this experiment. Takes either a Neuron object or just a Neuron id"""
        try:
            return neuron.ratepdf(trange=self.trange, **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].ratepdf(trange=self.trange, **kwargs) # neuron is probably a Neuron id
    ratepdf.__doc__ += '\n\n**kwargs:'
    ratepdf.__doc__ += '\nNeuron.RatePDF: '+getargstr(neuron.RatePDF.__init__)
    ratepdf.__doc__ += '\nNeuron.rate: '+getargstr(neuron.Neuron.rate)
    ratepdf.__doc__ += neuron.Neuron._rateargs


class RevCorrs(object):
    """Base class for doing reverse correlation of multiple neurons to a simulus"""
    def __init__(self, neurons=None, experiment=None, trange=None, nt=10):
        self.neurons = neurons
        self.experiment = experiment
        if trange == None:
            self.trange = self.experiment.trange
        else:
            self.trange = trange
        self.nt = nt # number of revcorr timepoints
        self.tis = range(0, nt, 1) # revcorr timepoint indices
        self.t = [ intround(ti * self.experiment.e.dynamic.sweepSec * 1000) for ti in self.tis ] # revcorr timepoint values, in ms

    def plot(self, interp='nearest', normed=True, title='ReceptiveFieldFrame', scale=2.0):
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
            rf = cmap(rf)[::, ::, ::, 0:3] # convert normalized luminance to RGB via the colormap, throw away alpha channel (not used for now in ReceptiveFieldFrame)
            rf = rf * 255 # scale up to 8 bit values
            rf = rf.round().astype(np.uint8) # downcast from float to uint8 for feeding to ReceptiveFieldFrame
            rfs.append(rf)
        frame = ReceptiveFieldFrame(title=title, rfs=rfs, neurons=self.neurons, t=self.t, scale=scale)
        frame.Show()

class STAs(RevCorrs):
    """Just a container class for multiple Neuron.STA objects. The plot() method is unique
    though: it plots all the Neuron.STA objects in a single figure"""
    def calc(self):
        self.stas = [] # store STAs in a list
        for neuron in self.neurons:
            stao = neuron.sta(experiment=self.experiment, trange=self.trange, nt=self.nt)
            self.stas.append(stao)

    def plot(self, interp='nearest', normed=True, scale=2.0):
        super(STAs, self).plot(interp=interp, normed=normed,
                               title=lastcmd(),
                               scale=scale)
    plot.__doc__ = RevCorrs.plot.__doc__


class STCs(RevCorrs):
    """Just a container class for multiple Neuron.STC objects. The plot() method is unique
    though: it plots all the Neuron.STC objects in a single figure"""
    def calc(self):
        self.stcs = [] # store STCs in a list
        for neuron in self.neurons:
            stco = neuron.stc(experiment=self.experiment, trange=self.trange, nt=self.nt)
            self.stcs.append(stco)

    def plot(self, interp='nearest', normed=True, scale=2.0):
        super(STCs, self).plot(interp=interp, normed=normed,
                               title=lastcmd(),
                               scale=scale)
    plot.__doc__ = RevCorrs.plot.__doc__


class ExperimentRevCorr(BaseExperiment):
    """Mix-in class that defines the reverse correlation related experiment methods"""
    def sta(self, neurons=None, **kwargs):
        """Returns an STAs RevCorrs object"""
        if neurons == None: # no Neurons were passed, use all the Neurons from the default Sort for this experiment's Recording
            keyvals = self.r.n.items() # get key val pairs in a list of tuples
            keyvals.sort() # make sure they're sorted by key
            neurons = []
            for key, val in keyvals:
                neurons.append(val)
        else:
            try: # assume neurons is a Neuron id or list of Neuron ids, get the associated Neuron objects from the default Sort for this experiment's Recording
                neurons = [ self.r.n[ni] for ni in toiter(neurons) ]
            except KeyError: # neurons is probably a list of Neuron objects
                pass
        staso = STAs(neurons=neurons, experiment=self, **kwargs) # init a new STAs object
        staso.calc()
        return staso
    sta.__doc__ += '\n\n**kwargs:\n'
    sta.__doc__ += getargstr(STAs.__init__)

    def stc(self, neurons=None, **kwargs):
        """Returns an STCs RevCorrs object"""
        if neurons == None: # no neurons were passed, use all the neurons from the default sort for this experiment's recording
            keyvals = self.r.n.items() # get key val pairs in a list of tuples
            keyvals.sort() # make sure they're sorted by key
            neurons = []
            for key, val in keyvals:
                neurons.append(val)
        else:
            try: # assume neurons is a neuron id or list of neuron ids, get the associated neuron objects from the default sort for this experiment's recording
                neurons = [ self.r.n[ni] for ni in tolist(neurons) ]
            except KeyError: # neurons is probably a list of Neuron objects
                pass
        stcso = STCs(neurons=neurons, experiment=self, **kwargs) # init a new STCs object
        stcso.calc()
        return stcso
    stc.__doc__ += '\n\n**kwargs:\n'
    stc.__doc__ += getargstr(STCs.__init__)


class ExperimentPopulationRaster(PopulationRaster):
    """A population raster limited to a single experiment"""
    def __init__(self, experiment):
        super(ExperimentPopulationRaster, self).__init__(recording=experiment.r, experiments={experiment.id: experiment})


class ExperimentRaster(BaseExperiment):
    """Mix-in class that defines the raster related experiment methods"""
    def raster(self, **kwargs):
        """Creates a population spike raster plot"""
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
    """Inherits all the experiment classes into a single one"""
    pass
