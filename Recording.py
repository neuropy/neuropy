"""Defines the Recording class. If need be,  This could be split up into multiple base classes
(one for basic init and loading RecordingBase), one for rates (RecordingRate), one for codes (RecordingCodes), etc...), and then combined
using multiple inheritance into a single child Class called Recording"""

print 'importing Recording'

from Core import *

class Recording(object):
    """A Recording corresponds to a single SURF file, ie everything recorded between when
    the user hits record and when the user hits stop and closes the SURF file, including any
    pauses in between Experiments within that Recording. A Recording can have multiple Experiments,
    and multiple spike extractions, called Rips"""
    def __init__(self, id=None, name=None, parent=Track):
        self.level = 3 # level in the hierarchy
        self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
        try:
            self.t = parent() # init parent Track object
        except TypeError: # parent is an instance, not a class
            self.t = parent # save parent Track object
        if id is not None:
            name = self.id2name(self.t.path, id) # use the id to get the name
        elif name is not None:
            id = self.name2id(name) # use the name to get the id
        else:
            raise ValueError, 'recording id and name can\'t both be None'
        self.id = id
        self.name = name
        self.path = self.t.path + self.name + SLASH
    def tree(self):
        """Print tree hierarchy"""
        print self.treebuf.getvalue(),
    def writetree(self,string):
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
        self.e = {} # store Experiments in a dictionary
        experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.din') ] # returns din filenames without their .din extension
        for (experimentid, experimentName) in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
            experiment = Experiment(id=experimentid, name=experimentName, parent=self) # pass both the id and the name
            experiment.load() # load the Experiment
            self.e[experiment.id] = experiment # save it
        #if len(self.e) == 1:
        #   self.e = self.e.values[0] # pull it out of the dictionary
        self.rip = {} # store Rips in a dictionary
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
        #if len(self.rip) == 1:
        #   self.rip = self.rip.values[0] # pull it out of the dictionary

    def plot_codeProbScatter(self, nbits=DEFAULTCODEBITLENGTH, randomneurons=False, **kwargs):
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
