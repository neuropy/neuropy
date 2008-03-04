"""Defines the Experiment class and all of its support classes."""

#print 'importing Experiment'

import dimstim.Experiment

from Core import *
from Core import _data
import Neuron
from Recording import PopulationRaster, Codes, CodeCorrPDF

class BaseExperiment(object):
    """An Experiment corresponds to a single contiguous VisionEgg stimulus session.
    It contains information about the stimulus during that session, including
    the DIN values and the text header. For data generated with dimstim >= 0.16,
    it includes the entire dimstim.Experiment object as an attribute (.e)
    a neuropy.Experiment is basically a container for a dimstim.Experiment"""

    from Recording import Recording

    def __init__(self, id=None, name=None, parent=None):
        self.level = 4 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        if parent:
            self.r = parent # save parent Recording object
        else:
            self.r = Recording() # init parent Recording object
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

    # doesn't need an id2name or name2id method, neither can really be derived from the other in an easy way (although could use re), the id is just chronological (which is also alphabetical) order, at least for now

    def load(self):

        from Recording import Recording

        f = file(os.path.join(self.path, self.name) + '.din', 'rb') # open the din file for reading in binary mode
        self.din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to nrows x 2 columns
        f.close()
        try:
            f = file(os.path.join(self.path, self.name) + '.textheader', 'r') # open the textheader file for reading
            self.textheader = f.read() # read it all in
            f.close()
        except IOError:
            warn('Error reading: <%s>, text header not loaded' % f.name)
            self.textheader = '' # set to empty

        treestr = self.level*TAB + self.name + '/'
        self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

        if self.textheader: # if it isn't empty
            names1 = locals().copy() # namespace before execing the textheader
            # should this instead be done line by line, using some regexps, as in dimstim?????????????
            exec(self.textheader) # execute the textheader as Python code. TODO: maybe add some checks here to prevent changes to filesystem from accidental code, unload the os or sys modules or something? But that wouldn't prevent 'accidental' code from 'accidentally' re-importing such modules
            names2 = locals().copy() # namespace after
            newnames = [ n2 for n2 in names2 if n2 not in names1 and n2 != 'names1' ] # names that were added to the namespace, excluding the 'names1' name itself
            for newname in newnames:
                self.__setattr__(newname, eval(newname)) # for each variable that was defined in the textheader, bind it as an attribute of this Experiment
            try:
                self.__version__ # dimtim up to ptc15 didn't have a version, neither did NVS display
            except AttributeError:
                self.__version__ = 0.0
            if self.__version__ >= 0.16: # after major refactoring of dimstim
                self.sweeptable = Core.SweepTable(experiment=self) # build the sweep table
                self.st = self.sweeptable.data # synonym, used a lot by Experiment subclasses
                self.xorig = deg2pix(self.static.xorigDeg) + I.SCREENWIDTH / 2
                self.yorig = deg2pix(self.static.yorigDeg) + I.SCREENHEIGHT / 2

                self.REFRESHTIME = intround(1 / float(self.I.REFRESHRATE) * 1000000) # in us, keep 'em integers
            else:
                self.loadCat15exp()
        else:
            self.REFRESHTIME = self.din[1, 0] - self.din[0, 0] # use the time difference between the first two din instead

        self.trange = (self.din[0, 0], self.din[-1, 0]+self.REFRESHTIME) # add an extra refresh time after last din, that's when screen actually turns off

    def loadCat15exp(self):

        ## TODO: - fake a .e dimstim.Experiment object, to replace what used to be the .stims object for Movie experiments
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
        from Movie import Movie

        try:
            self.stims = unique(self.playlist) # self.stims is a non-repeating list of object oriented stim objects (Movie is the only possible type) in this Experiment
        except AttributeError: # this was a simple non object-oriented stim, has no playlist
            self.stims = []
        for s in self.stims:
            s.e = self # If you inited stim object(s) (like a movie) while execing the textheader, you didn't have a chance to pass this exp as the parent in the init. So just set the attribute manually
            try: # this'll probably only apply to Movies stim, cuz others won't have fnames
                s.name = os.path.splitext(s.fname)[0] # extensionless fname, fname should've been defined in the textheader
                if s.name not in _data.movies: # and it very well may not be, cuz the textheader inits movies with no args, leaving fname==None at first, which prevents it from being added to _data.movies
                    _data.movies[s.name] = s # add s to _data.movies dictattr
            except AttributeError: # s.fname doesn't exist? probably not a Movie stim
                pass
            # Search self.moviepath string (from textheader) for 'Movies' word (preferably case insensitive). Everything after that is the relative path to your base movies folder. Eg, if self.moviepath = 'C:\\Desktop\\Movies\\reliability\\e\\single\\', then set self.relpath = '\\reliability\\e\\single\\'
            spath = splitpath(self.moviepath)
            try:
                matchi = spath.index('Movies')
            except ValueError:
                matchi = spath.index('movies')
            s.relpath = joinpath(spath[matchi+1 ::])
            s.path = os.path.join(MOVIEPATH, s.relpath)

        # Generate the sweeptable here, no need to load it from files anymore...
        # self.sweeptable = {[]} # dictionary of lists, ie sweeptable={'ori':[0,45,90], 'sfreq':[1,1,1]}
        # so you index into it with self.sweeptable['var'][sweepi]
        # vars = self.sweeptable.keys()
        # need to check if varlist exists, if so use it (we're dealing with Cat 15), if not, use revamped dimstim.SweepTable class
        if len(self.stims) > 0: # this Experiment has object-oriented stim(s)
            for s in self.stims:
                varvals = {} # init a dictionary that will contain variable values
                for var in s.varlist:
                    varvals[var] = eval('s.' + var) # generate a dictionary with var:val entries, to pass to buildSweepTable
                s.sweepTable = dimstim.Core.buildSweepTable(s.varlist, varvals, s.nruns, s.shuffleRuns, s.blankSweep, s.shuffleBlankSweeps, makeSweepTableText=0)[0] # passing varlist by reference, dim indices end up being modified
        else: # this is a simple stim (not object oriented)
            varvals = {} # init a dictionary that will contain variable values
            for var in self.varlist:
                varvals[var] = eval('self.' + var) # generate a dictionary with var:val entries, to pass to buildSweepTable
            self.sweepTable = self.buildCat15SweepTable(self.varlist, varvals, self.nruns, self.shuffleRuns, self.blankSweep, self.shuffleBlankSweeps, makeSweepTableText=0)[0] # passing varlist by reference, dim indices end up being modified

        # for all (Movie) stims inited by the textheader, enter each Movie into _data.movies
        for s in self.stims:
            assert s.__class__ == Movie
            if s.name not in _data.movies:
                _data.movies[s.name] = s # add the Movie stim to the movies dictattr, with the extensionless fname as the key
                #treestr = s.level*TAB + os.path.join(s.path, s.fname)
                #self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
        try:
            self.REFRESHTIME = intround(1 / float(self.REFRESHRATE) * 1000000) # in us, keep 'em integers
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
        if nruns.__class__ != int or nruns < 0:
            raise ValueError, 'nruns must be a non-negative integer'
        if shuffleRuns.__class__ != int or shuffleRuns not in (0, 1):
            raise ValueError, 'shuffleRuns must be 0 or 1'
        if blankSweep.__class__ != tuple or blankSweep[0].__class__ != int or blankSweep[0] == 1 or blankSweep[0] < 0:
            raise ValueError, 'blankSweep must be a tuple, and blankSweep[0] must be a non-negative integer and cannot equal 1'
        if shuffleBlankSweeps.__class__ != int or shuffleBlankSweeps not in (0, 1):
            raise ValueError, 'shuffleBlankSweeps must be 0 or 1'
        # check that shuffle flag is valid for each variable
        for var in varlist:
            if varlist[var]['shuffle'].__class__ != int or varlist[var]['shuffle'] not in (0, 1, 2):
                raise ValueError, 'Variable %s has a shuffle value outside of (0, 1, 2)' % var

        for (key, val) in varvals.items():
            if val == None:
                raise ValueError, 'Variable %s was passed with value None. Can\'t generate a sweeptable with None values' % key
            if val.__class__ != list:
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
                    if eval(var).__class__ == list:
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
    """Mix-in class that defines the spike code related Experiment methods"""
    def code(self, neuron=None, **kwargs):
        """Returns a Neuron.Code object, constraining it to the time range of this Experiment. Takes either a Neuron object or just a Neuron id"""
        try:
            return neuron.code(tranges=[self.trange], **kwargs) # see if neuron is a Neuron
        except AttributeError:
            return self.r.n[neuron].code(tranges=[self.trange], **kwargs) # neuron is probably a Neuron id
    code.__doc__ += '\n\n**kwargs:'
    code.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code)
    code.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__)

    def codes(self, neurons=None, **kwargs):
        """Returns a 2D array where each row is a neuron code constrained to the time range of this Experiment"""
        if neurons == None:
            neurons = self.r.n
        codeso = self.r.codes(neurons=neurons, experiments=[self], **kwargs)
        codeso.calc()
        return codeso
    codes.__doc__ += '\n\n**kwargs:'
    codes.__doc__ += '\nCodes: '+getargstr(Codes.__init__)
    codes.__doc__ += '\nNeuron.code: '+getargstr(Neuron.Neuron.code)
    codes.__doc__ += '\nbinary: '+getargstr(Neuron.BinaryCode.__init__)

    def codecorr(self, neuron1, neuron2, **kwargs):
        """Calculates the correlation of two Neuron.Code objects"""
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
        cco = self.r.codecorrpdf(experiments={self.id: self}, **kwargs) # init a new one
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
    """Mix-in class that defines the reverse correlation related Experiment methods"""
    def sta(self, neurons=None, **kwargs):
        """Returns an STAs RevCorrs object"""
        if neurons == None: # no Neurons were passed, use all the Neurons from the default Rip for this experiment's Recording
            keyvals = self.r.n.items() # get key val pairs in a list of tuples
            keyvals.sort() # make sure they're sorted by key
            neurons = []
            for key, val in keyvals:
                neurons.append(val)
        else:
            try: # assume neurons is a Neuron id or list of Neuron ids, get the associated Neuron objects from the default Rip for this experiment's Recording
                neurons = [ self.r.n[ni] for ni in tolist(neurons) ]
            except KeyError: # neurons is probably a list of Neuron objects
                pass
        staso = STAs(neurons=neurons, experiment=self, **kwargs) # init a new STAs object
        staso.calc()
        return staso
    sta.__doc__ += '\n\n**kwargs:\n'
    sta.__doc__ += getargstr(STAs.__init__)

    def stc(self, neurons=None, **kwargs):
        """Returns an STCs RevCorrs object"""
        if neurons == None: # no Neurons were passed, use all the Neurons from the default Rip for this experiment's Recording
            keyvals = self.r.n.items() # get key val pairs in a list of tuples
            keyvals.sort() # make sure they're sorted by key
            neurons = []
            for key, val in keyvals:
                neurons.append(val)
        else:
            try: # assume neurons is a Neuron id or list of Neuron ids, get the associated Neuron objects from the default Rip for this experiment's Recording
                neurons = [ self.r.n[ni] for ni in tolist(neurons) ]
            except KeyError: # neurons is probably a list of Neuron objects
                pass
        stcso = STCs(neurons=neurons, experiment=self, **kwargs) # init a new STCs object
        stcso.calc()
        return stcso
    stc.__doc__ += '\n\n**kwargs:\n'
    stc.__doc__ += getargstr(STCs.__init__)


class ExperimentPopulationRaster(PopulationRaster):
    """A population raster limited to a single Experiment"""
    def __init__(self, experiment):
        super(ExperimentPopulationRaster, self).__init__(recording=experiment.r, experiments={experiment.id: experiment})


class ExperimentRaster(BaseExperiment):
    """Mix-in class that defines the raster related Experiment methods"""
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
    """Inherits all the Experiment classes into a single Experiment class"""
    pass
