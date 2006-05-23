'''Core neuropy functions and classes

          object hierarchy:     level:

                Data              0
                  |
                 Cat              1
                  |
                Track             2
                  |
              Recording           3
             /         \\
       Experiment      Rip        4
            |           |
          Movie       Neuron      5
'''

DEFAULTDATAPATH = 'C:/data/' # the convention in neuropy will be that all 'path' var names have a trailing slash
DEFAULTCATID    = 15
DEFAULTTRACKID  = '7c'
RIPKEYWORDS = ['best', 'liberal spikes'] # a Rip with one of these keywords (listed in decreasing priority) will be loaded as the default Rip for its Recording
SLASH = '/' # use forward slashes instead of having to use double backslashes
TAB = '    ' # 4 spaces

DEFAULTMOVIEPATH = 'C:/pub/Movies/'
DEFAULTMOVIENAME = 'mseq32.m'

import os, types, pprint, numpy, numarray, struct, re, StringIO, sys
import Dimstim.Movies
from Dimstim.Core import buildSweepTable

pp = pprint.pprint
INF = numpy.inf

# generate sweeptables on the fly in Experiment.load()
# need to delete the extra lines in all the textheaders, and uncomment some lines too!
# Rips should really have ids to make them easier to reference to: r[83].rip[0] instead of r[83].rip['conservative spikes'] - this means adding id prefixes to rip folder names (or maybe suffixes: 'conservative spikes.0.rip', 'liberal spikes.1.rip', etc...). Prefixes would be better cuz they'd force sorting by id in explorer (which uses alphabetical order) - ids should be 0-based of course
# worry about conversion of ids to strings: some may be only 1 digit and may have a leading zero!
# maybe make two load() f'ns for Experiment and Neuron: one from files, and a future one from a database
# make a save() f'n that pickles the object (including any of its results, like its STA, tuning curve points, etc)?

"""
def str2(data):
	if type(data) is types.IntTypes:
		s = str(data)
		if len(s) == 1:
			s = '0'+s # add a leading zero for single digits
"""

def txtdin2binarydin(fin, fout):
	'''Converts a csv text .din file to a int64 binary .din file'''
	fi = file(fin, 'r') # open the din file for reading in text mode
	fo = file(fout,'wb') # for writing in binary mode
	for line in fi:
		line = line.split(',')
		'''
		# for old NVS display, converts from NVS condition numbers (which increment with repeats) to Dimstim sweepis (which don't)
		nruns = 18
		line[1] = int(line[1]) % nruns
		'''
		fo.write( struct.pack('@qq',int(line[0]),int(line[1])) ) # read both values in as a C long longs, using the system's native ('@') byte order
	fi.close()
	fo.close()
	print 'Converted ascii din: ', fin, ' to binary din: ', fout
def warn(msg,level=2,exit_val=1):
	'''Standard warning printer. Gives formatting consistency. Stolen from IPython.genutils'''
	if level>0:
		header = ['','','WARNING: ','ERROR: ','FATAL ERROR: ']
		print >> sys.stderr, '%s%s' % (header[level],msg)
		if level == 4:
			print >> sys.stderr,'Exiting.\n'
			sys.exit(exit_val)
"""
def warn(msg):
	import warnings
	warnings.warn(msg, category=RuntimeWarning, stacklevel=2)
"""
"""
def unique(objlist):
	'''Does in-place removal of non-unique objects in a list of objects'''
	for (i,obj1) in enumerate(objlist):
		for (j,obj2) in enumerate(objlist):
			if i != j and obj1 == obj2:
				del objlist[j]
"""

def unique(inseq):
	'''Return unique items from a 1-dimensional sequence. Stolen from numpy.unique(), modified to return list instead of array'''
	# Dictionary setting is quite fast.
	outseq = {}
	for item in inseq:
		outseq[item] = None
	return list(outseq.keys())

def iterable(y):
	'''Check if the input is iterable, stolen from numpy.iterable()'''
	try: iter(y)
	except: return 0
	return 1

def histogramSorted(sorteda, bins=10, range=None, normed=False):
	'''Builds a histogram, stolen from numpy.histogram(), modified to assume sorted input'''
	a = numpy.asarray(sorteda).ravel()
	if not iterable(bins):
		if range is None:
			range = (a.min(), a.max())
		mn, mx = [mi+0.0 for mi in range]
		if mn == mx:
			mn -= 0.5
			mx += 0.5
		bins = numpy.linspace(mn, mx, bins, endpoint=False)
	#n = numpy.sort(a).searchsorted(bins)
	n = a.searchsorted(bins)
	n = numpy.concatenate([n, [len(a)]]) # don't understand what this does
	n = n[1:]-n[:-1]
	if normed:
		db = bins[1] - bins[0]
		return 1.0/(a.size*db) * n, bins
	else:
		return n, bins

"""
def unique(objlist):
	'''Returns the input list minus any repeated objects it may have had. Also defined in Dimstim'''
	return list(set(objlist)) # this requires Python >= 2.4
"""
"""
def tolist(obj):
	'''Takes either scalar or sequence input and returns a list,
	useful when you want to iterate over an object (like in a for loop),
	and you don't want to have to do type checking or handle exceptions
	when the object isn't a sequence'''
	try: # assume obj is a sequence
		return list(obj) # converts any sequence to a list
	except TypeError: # obj is probably a scalar
		return [obj] # converts any scalar to a list
"""
###########################

class Data(object): # use 'new-style' classes
	'''Data can have multiple Cats'''
	def __init__(self, dataPath=DEFAULTDATAPATH):
		self.level = 0 # level in the hierarchy
		self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		self.name = 'Data'
		self.path = dataPath
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		# Data has no parent to write to
	def load(self):
		treestr = self.level*TAB + self.name + '/'
		self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		self.c = {} # store Cats in a dictionary
		catNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Cat ') ] # os.listdir() returns all dirs AND files
		for catName in catNames:
			cat = Cat(id=None, name=catName, parent=self) # make an instance using just the catName (let it figure out the cat id)
			cat.load() # load the Cat
			self.c[cat.id] = cat # save it
		#if len(self.c) == 1:
		#	self.c = self.c.values[0] # pull it out of the dictionary

class Cat(object):
	'''A Cat can have multiple Tracks'''
	def __init__(self, id=DEFAULTCATID, name=None, parent=Data):
		self.level = 1 # level in the hierarchy
		self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		try:
			self.d = parent() # init parent Data object
		except TypeError: # parent is an instance, not a class
			self.d = parent # save parent Data object
		if id is not None:
			name = self.id2name(self.d.path,id) # use the id to get the name
		elif name is not None:
			id = self.name2id(name) # use the name to get the id
		else:
			raise ValueError, 'cat id and name can\'t both be None'
		self.id = id
		self.name = name
		self.path = self.d.path + self.name + SLASH
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.d.writetree(string)
	def id2name(self, path, id):
		name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith('Cat '+str(id)) ]
		if len(name) != 1:
			raise NameError, 'Ambiguous or non-existent Cat id: %s' % id
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		id = name.replace('Cat ','',1) # replace first occurrence of 'Cat ' with nothing, keep the rest
		if not id:
			raise NameError, 'Badly formatted Cat name: %s' % name
		try:
			id = int(id) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return id
	def load(self):
		treestr = self.level*TAB + self.name + '/'
		self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		self.t = {} # store Tracks in a dictionary
		trackNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Track ') ]
		for trackName in trackNames:
			track = Track(id=None, name=trackName, parent=self) # make an instance using just the track name (let it figure out the track id)
			track.load() # load the Track
			self.t[track.id] = track # save it
		#if len(self.t) == 1:
		#	self.t = self.t.values[0] # pull it out of the dictionary

class Track(object):
	'''A Track can have multiple Recordings'''
	def __init__(self, id=DEFAULTTRACKID, name=None, parent=Cat):
		self.level = 2 # level in the hierarchy
		try:
			self.c = parent() # init parent Cat object
		except TypeError: # parent is an instance, not a class
			self.c = parent # save parent Cat object
		if id is not None:
			name = self.id2name(self.c.path,id) # use the id to get the name
		elif name is not None:
			id = self.name2id(name) # use the name to get the id
		else:
			raise ValueError, 'track id and name can\'t both be None'
		self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		self.id = id
		self.name = name
		self.path = self.c.path + self.name + SLASH
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.c.writetree(string)
	def id2name(self, path, id):
		name = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith('Track '+str(id)) ]
		if len(name) != 1:
			raise NameError, 'Ambiguous or non-existent Track id: %s' % id
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		id = name.replace('Track ','',1) # replace first occurrence of 'Track ' with nothing, keep the rest
		if not id:
			raise NameError, 'Badly formatted Track name: %s' % name
		try:
			id = int(id) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return id
	def load(self):
		treestr = self.level*TAB + self.name + '/'
		self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		self.r = {} # store Recordings in a dictionary
		recordingNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname[0:2].isdigit() and dirname.count(' - ') == 1 ] # 1st 2 chars in dirname must be digits, must contain exactly 1 occurrence of ' - '
		for recordingName in recordingNames:
			recording = Recording(id=None, name=recordingName, parent=self) # make an instance using just the recording name (let it figure out the recording id)
			recording.load() # load the Recording
			self.r[recording.id] = recording # save it
		#if len(self.r) == 1:
		#	self.r = self.r.values[0] # pull it out of the dictionary

class Recording(object):
	'''A Recording corresponds to a single SURF file, ie everything recorded between when
	the user hits record and when the user hits stop and closes the SURF file, including any
	pauses in between Experiments within that Recording. A Recording can have multiple Experiments,
	and multiple spike extractions, called Rips'''
	def __init__(self, id=None, name=None, parent=Track):
		self.level = 3 # level in the hierarchy
		self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		try:
			self.t = parent() # init parent Track object
		except TypeError: # parent is an instance, not a class
			self.t = parent # save parent Track object
		if id is not None:
			name = self.id2name(self.t.path,id) # use the id to get the name
		elif name is not None:
			id = self.name2id(name) # use the name to get the id
		else:
			raise ValueError, 'recording id and name can\'t both be None'
		self.id = id
		self.name = name
		self.path = self.t.path + self.name + SLASH
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.t.writetree(string)
	def id2name(self, path, id):
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
		treestr = self.level*TAB + self.name + '/'
		self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		self.e = {} # store Experiments in a dictionary
		experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.din') ] # returns din filenames without their .din extension
		for (experimentid, experimentName) in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
			experiment = Experiment(id=experimentid, name=experimentName, parent=self) # pass both the id and the name
			experiment.load() # load the Experiment
			self.e[experiment.id] = experiment # save it
		#if len(self.e) == 1:
		#	self.e = self.e.values[0] # pull it out of the dictionary
		self.rip = {} # store Rips in a dictionary
		ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
		defaultRipNames = [ ripName for ripName in ripNames for ripkeyword in RIPKEYWORDS if ripName.count(ripkeyword) ]
		if len(defaultRipNames) < 1:
			warn('Couldn''t find a default Rip for Recording(%s)' % self.id)
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
		#	self.rip = self.rip.values[0] # pull it out of the dictionary

class Experiment(object):
	'''An Experiment corresponds to a single contiguous VisionEgg stimulus session.
	It contains information about the stimulus during that session, including
	the DIN values, the text header, and any Movies that were involved'''
	def __init__(self, id=None, name=None, parent=Recording): # Experiment IDs are 1-based in the .din filenames, at least for now. They should be renamed to 0-based. Here, they're treated as 0-based
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
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.r.writetree(string)
	# doesn't need a id2name or name2id method, neither can really be derived from the other in an easy way (although could use re), the id is just alphabetical order, at least for now
	def load(self):
		f = file(self.path + self.name + '.din', 'rb') # open the din file for reading in binary mode
		self.din = numpy.fromfile(f, dtype=numpy.int64).reshape(-1,2) # reshape to nrows x 2 columns
		f.close()
		try:
			f = file(self.path + self.name + '.textheader', 'r') # open the textheader file for reading
			self.textheader = f.read() # read it all in
			f.close()
		except IOError:
			warn('Error reading: <%s.textheader>, text header not loaded' % self.name)
			self.textheader = '' # set to empty

		treestr = self.level*TAB + self.name + '/'
		self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

		if self.textheader: # if it isn't empty
			names1 = locals() # namespace before execing the textheader should this be locals() or globals()????????????
			exec(self.textheader) # execute the textheader as Python code, maybe add some checks here to prevent changes to filesystem from accidental code, unload the os or sys modules or something?
			names2 = locals() # namespace after
			newnames = [ n2 for n2 in names2 for n1 in names1 if n2 != n1 ] # names that were added to the namespace
			for newname in newnames:
				self.__setattr__(newname, eval(newname)) # for each variable that was defined in the textheader, bind it as an attribute of this Experiment
			'''
			thlines = self.textheader.splitlines()
			for line in thlines:
				if line and not line.startswith('#'): # if it ain't blank and it ain't commented out
					exec(line) # just exec the line first so everyting on the rhs is defined in the unbound namespace
					exec('self.'+line) # exec the statement on each line in the text header and save it as an attribute of this Experiment
			'''
			# finding the movie object name:
			# self.textheader.index(' = Movie()')
			# or look for self.textheader.index('.oname = ')
			# once you have the name in a string, run eval('self.'+oname)?
			'''
			# or, run something that returns all objects of type Movie in the attributes of self:
			self.movies = [] # stores all the movies inited by the textheader
			for objname in self.__dict__:
				obj = eval('self.'+objname)
				if isinstance(obj, Movie):
					self.movies.append(obj)
			'''
			try:
				self.stims = unique(self.playlist) # self.stims is a non-repeating list of object oriented stim objects (Movie is the only possible kind right now) in this Experiment
			except AttributeError: # this was a simple non object-oriented stim, has no playlist
				self.stims = []
			#if len(self.stims) == 1:
			#	self.stims = self.stims[0] # get rid of the list

			# if you inited a stim object(s) (like a movie) while execing the textheader, you didn't have a chance to pass this exp as the parent in the init. So just set the attribute manually:

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
				# also, if you initd a stim that needs to be loaded (like a movie), maybe you should also load it now (this wasn't done when execing the textheader)
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
					s.sweepTable = buildSweepTable(s.varlist, varvals, s.nruns, s.shuffleRuns, s.blankSweep, s.shuffleBlankSweeps, makeSweepTableText=0) # passing varlist by reference, dim indices end up being modified
			else: # this is a simple stim (not object oriented)
				varvals={} # init a dictionary that will contain variable values
				for var in self.varlist:
					varvals[var]=eval('self.'+var) # generate a dictionary with var:val entries, to pass to buildSweepTable
				self.sweepTable = buildSweepTable(self.varlist, varvals, self.nruns, self.shuffleRuns, self.blankSweep, self.shuffleBlankSweeps, makeSweepTableText=0) # passing varlist by reference, dim indices end up being modified

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
					if s.name == defaultm.name and isinstance(s,Movie):
						if defaultm.data == None: # see if this default movie has yet to be loaded
							defaultm.load() # load this default movie
						s.data = defaultm.data # point this movie's data to default movie data
						treestr = s.level*TAB + s.name
						self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen

		self.buildsweepranges()
	def buildsweepranges(self):
		pass
		self.sweepranges = []

class Rip(object):
	'''A Rip is a single spike extraction. Generally, Rips of the same name within the same Track
	were generated with the same spike template, though of course Rips in different Tracks must
	be generated from different templates, even if the Rips have the same name'''
	def __init__(self, id=None, name=None, parent=Recording):
		self.level = 4 # level in the hierarchy
		#self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		try:
			self.r = parent() # init parent Recording object
		except TypeError: # parent is an instance, not a class
			self.r = parent # save parent Recording object
		if name is None:
			raise ValueError, 'rip name can\'t be None'
		# rips don't have ids, at least for now. Just names
		self.id = id # not really used by the Rip class, just there for user's info
		self.name = name
		self.path = self.r.path + self.name + '.rip' + SLASH # have to add .rip extension to rip name to get its actual folder name
	"""
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.r.writetree(string)
	"""
	def load(self):
		#treestr = self.level*TAB + self.name + '/'
		#self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		self.n = {} # store Neurons in a dictionary
		neuronNames = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.spk') ] # returns spike filenames without their .spk extension
		for neuronName in neuronNames:
			neuron = Neuron(id=None, name=neuronName, parent=self) # make an instance using just the neuron name (let it figure out the neuron id)
			neuron.load() # load the neuron
			self.n[neuron.id] = neuron # save it
		# then, maybe add something that loads info about the rip, say from some file describing the template used, and all the thresholds, exported to the same folder by SURF
		# maybe also load the template used for the rip, perhaps also stored in the same folder

class Neuron(object):
	'''A Neuron object''s spike data spans all the Experiments within a Recording.
	If different Recordings have Rips with the same name, you can assume that the
	same spike template was used for all of those Recordings, and that therefore
	the neuron ids are the same'''
	def __init__(self, id=None, name=None, parent=Rip): # neuron names don't include the '.spk' ending, although neuron filenames do
		self.level = 5 # level in the hierarchy
		#self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		try:
			self.rip = parent() # init parent Rip object
		except TypeError: # parent is an instance, not a class
			self.rip = parent # save parent Rip object
		if id is not None:
			name = self.id2name(self.rip.path,id) # use the id to get the name
		elif name is not None:
			id = self.name2id(name) # use the name to get the id
		else:
			raise ValueError, 'neuron id and name can\'t both be None'
		self.id = id
		self.name = name
		self.path = self.rip.path
	"""
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.rip.writetree(string)
	"""
	def id2name(self, path, id):
		name = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(path) if os.path.isfile(path+fname) and \
		       ( fname.find('_t'+str(id)+'.spk')!=-1 or fname.find('_t0'+str(id)+'.spk')!=-1 or fname.find('_t00'+str(id)+'.spk')!=-1 ) ] # have to deal with leading zero ids, go up to 3 digit ids, should really use a re to do this properly...
		if len(name) != 1:
			raise NameError, 'Ambiguous or non-existent Neuron id: %s' % id
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		try:
			id = name[name.rindex('_t')+2::] # everything from just after the last '_t' to the end of the neuron name, index() raises ValueError if it can't be found
		except ValueError:
			raise ValueError, 'Badly formatted Neuron name: %s' % name
		try:
			id = int(id) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return id
	def load(self): # or loadspikes()?
		f = file(self.path + self.name + '.spk', 'rb') # open the spike file for reading in binary mode
		self.spikes = numpy.fromfile(f, dtype=numpy.int64) # read in all spike times in us
		f.close()
		self.results = {} # a dictionary to store results in
		#treestr = self.level*TAB + self.name + '/'
		#self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
	def cut(self, tstart=None, tend=None): # maybe use a masked array instead
		'''Returns all of the Neuron's spike times that fall within tstart and tend'''
		#return self.spikes[ numpy.greater_equal(self.spikes, tstart) & numpy.less_equal(self.spikes, tend) ]
		if tstart is None:
			tstart = self.spikes[0]
		if tend is None:
			tend = self.spikes[-1]
		'''
		# this is what we're trying to do:
		return self.spikes[ (self.spikes >= tstart) & (self.spikes <= tend) ]

		self.searchsorted(values) method does it faster. It returns an index where the value would fit in self. The index is such that self[index-1] < value <= self[index]. In this formula self[self.size]=inf and self[-1]= -inf
		'''
		lo, hi = self.spikes.searchsorted([tstart, tend]) # returns slice indices
		return self.spikes[ lo : hi ] # slice it
	def cutrel(self, tstart=None, tend=None):
		'''Cuts Neuron spike times and returns them relative to tstart'''
		if tstart is None:
			tstart = self.spikes[0]
		if tend is None:
			tend = self.spikes[-1]
		if tstart.__class__ is types.FloatType:
			warn('Converting tstart to int: '+str(int(tstart)))
		return self.cut(tstart, tend) - int(tstart)
	def rate(self, tres=50000, trange=None, method='bin'):
		if method == 'bin':
			'''bins sequence demarcates left bin edges'''
			if trange==None:
				trange = (self.spikes[0], self.spikes[1])
			bins = arange( trange[0], trange[1], tres )
			n, bins = numpy.histogram(self.spikes, bins=bins, normed=1)
			n = n / float(tres)
			self.results['rate'] = (n, bins)
			return n
	def raster(self):
		pass

class Movie(Dimstim.Movies.Movie): # inherit from Dimstim Movie() class (assumes it's new-style)
	'''A Movie stimulus object'''
	def __init__(self, name=None, path=DEFAULTMOVIEPATH, parent=None):
		super(Movie, self).__init__() # first run __init__() of inherited Dimstim Movie class
		self.level = 5 # level in the hierarchy
		try:
			self.e = parent() # init parent Experiment object
		except TypeError: # parent is an instance, not a class
			self.e = parent # save parent Experiment object
		self.name = name
		self.path = path
	def load(self):
		# Load movie data
		f = file(self.path + self.name, 'rb') # open the movie file for reading in binary format
		headerstring = f.read(5)
		if headerstring == 'movie': # a header has been added to the start of the file
			(self.ncellswide,) = struct.unpack('H', f.read(2)) # 'H'== unsigned short int == 2 bytes on this PC
			(self.ncellshigh,) = struct.unpack('H', f.read(2))
			(self.nframes,) = struct.unpack('H', f.read(2))
			if self.nframes == 0: # this was used in Cat 15 mseq movies to indicate 2**16 frames
				self.nframes = 2**16
			offset = 11 # header is this long
		else: # there's no header at the start of the file, set the file pointer back to the beginning and use these hard coded values:
			f.seek(0)
			self.ncellswide = self.ncellshigh = 64
			self.nframes = 6000
			offset = 0 # header is this long
		# read in all of the frames
		self.data = numpy.fromfile(f, numpy.uint8, count=self.nframes*self.ncellshigh*self.ncellswide).reshape(self.nframes,self.ncellshigh,self.ncellswide) # read it all in
		#self.data = numarray.fromfile(f, numpy.UInt8, (self.nframes,self.ncellshigh,self.ncellswide)) # read it all in
		leftover = f.read() # check if there are any leftover bytes in the file
		if leftover != '':
			pp(leftover)
			print self.nframes,self.ncellshigh,self.ncellswide
			raise RuntimeError('There are unread bytes in movie file \'%s\'. Width, height, or nframes is incorrect in the movie file header.' % self.name)
		#self.data = self.data[::,::-1,::] # flip the movie frames vertically for OpenGL's bottom left origin
		f.close() # close the movie file

################

# init some typical movies (but don't load 'em til needed), then just point to them within the appropriate Experiments
MSEQ32 = Movie(name='mseq32.m', parent=None)
MSEQ16 = Movie(name='mseq16.m', parent=None)
# shouldn't use sparse bar movies anymore, can access VisionEgg directly now, get the framebuffers to directly do STA
#sparsebars = Movie(path='C:/data/Cat 15/Track 7c/72 - track 7c sparseexps/', name='72 - track 7c sparseexps.sparsebars.movie');
'''
# init and load a Recording to play around with
print 'Initing and loading Recording(92):'
r=Recording(92)
r.load()
'''
print 'Initing and loading Recording(71):'
r=Recording(71)
r.load()
