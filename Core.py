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
DEFAULTRIPNAME  = 'liberal spikes' # a rip is the name of a spike sorting extraction; rips with the same name should have been done with the same templates, and possibly the same ripping thresholds
SLASH = '/' # use forward slashes instead of having to use double backslashes
TAB = '    ' # 4 spaces

DEFAULTMOVIEPATH = 'C:/pub/Movies/'
DEFAULTMOVIENAME = 'mseq32.m'

import os, types, pprint, numpy, numarray, struct, re, StringIO, sys, warnings

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
		#print line[0], line[1]
		fo.write( struct.pack('@qq',int(line[0]),int(line[1])) ) # read both values in as a C long longs, using the system's native ('@') byte order
	fi.close()
	fo.close()
	print 'Converted ascii din: ', fin, ' to binary din: ', fout

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
			raise NameError, 'Ambiguous or non-existent Cat id: '+str(id)
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		id = name.replace('Cat ','',1) # replace first occurrence of 'Cat ' with nothing, keep the rest
		if not id:
			raise NameError, 'Badly formatted Cat name: '+name
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
			raise NameError, 'Ambiguous or non-existent Track id: '+str(id)
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		id = name.replace('Track ','',1) # replace first occurrence of 'Track ' with nothing, keep the rest
		if not id:
			raise NameError, 'Badly formatted Track name: '+name
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
			raise NameError, 'Ambiguous or non-existent Recording id: '+str(id)
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		try:
			id = name[0:name.index(' - ')] # everything before the first ' - ', index() raises ValueError if it can't be found
		except ValueError:
			raise ValueError, 'Badly formatted Recording name: '+name
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
		self.rip = {} # store Rips in a dictionary
		ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
		for ripName in ripNames:
			rip = Rip(name=ripName, parent=self) # pass just the name, ain't no such thing as a ripid, at least for now
			rip.load() # load the Rip
			self.rip[rip.name] = rip # save it
		try: # make the Neurons from the default rip (if it exists in the Recording path) available in the Recording, so you can access them via r.n[nid] instead of having to do r.rip[name].n[nid]. Make them just another pointer to the data in r.rip[DEFAULTRIPNAME].n
			self.n = self.rip[DEFAULTRIPNAME].n
			self.defaultRipPath = self.path + DEFAULTRIPNAME + SLASH
		except:
			pass

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
		f = file(self.path + self.name + '.textheader', 'r') # open the textheader file for reading
		self.textheader = f.read() # read it all in
		f.close()
		# then, for each line in the textheader, exec() it so you get self.varname saved directly within in the Experiment object - watch out, will try and make Movie() objects and Bar() objects, etc?
		# also need to generate sweeptable here
		treestr = self.level*TAB + self.name + '/'
		self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		self.m = [] # Experiments can potentially have multiple movies
		for m in [mseq32, mseq16]: # check if this Experiment uses specific movies
			if self.textheader.count(m.name): # search the textheader for the movie's filename
				try:
					m.data # see if this movie has already been loaded
				except AttributeError: # load this movie
					m.load()
				self.m.append(Movie(parent=self)) # add a new movie to this Experiment
				self.m[-1].data = m.data # point data to already loaded data in specific movie
				treestr = m.level*TAB + m.name
				self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
		if len(self.m) == 1:
			self.m = self.m[0] # get rid of the list if there's only a single movie in this Experiment (usual case)

class Rip(object):
	'''A Rip is a single spike extraction. Generally, Rips of the same name within the same Track
	were generated with the same spike template, though of course Rips in different Tracks must
	be generated from different templates, even if the Rips have the same name'''
	def __init__(self, name=DEFAULTRIPNAME, parent=Recording):
		self.level = 4 # level in the hierarchy
		self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
		try:
			self.r = parent() # init parent Recording object
		except TypeError: # parent is an instance, not a class
			self.r = parent # save parent Recording object
		if name is None:
			raise ValueError, 'rip name can\'t be None'
		# rips don't have ids, at least for now. Just names
		self.name = name
		self.path = self.r.path + self.name + '.rip' + SLASH # have to add .rip extension to rip name to get its actual folder name
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.r.writetree(string)
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
	'''A Neuron object's spike data spans all the Experiments within a Recording.
	If different Recordings have Rips with the same name, you can assume that the
	same spike template was used for all of those Recordings, and that therefore
	the neuron ids are the same'''
	def __init__(self, id=None, name=None, parent=Rip): # neuron names don't include the '.spk' ending, although neuron filenames do
		self.level = 5 # level in the hierarchy
		self.treebuf = StringIO.StringIO() # create a string buffer to print tree hierarchy to
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
	def tree(self):
		'''Print tree hierarchy'''
		print self.treebuf.getvalue(),
	def writetree(self,string):
		'''Write to self's tree buffer and to parent's too'''
		self.treebuf.write(string)
		self.rip.writetree(string)
	def id2name(self, path, id):
		name = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(path) if os.path.isfile(path+fname) and \
		             ( fname.find('_t'+str(id)+'.spk')!=-1 or fname.find('_t0'+str(id)+'.spk')!=-1 or fname.find('_t00'+str(id)+'.spk')!=-1 ) ] # have to deal with leading zero ids, go up to 3 digit ids, should really use a re to do this properly...
		if len(name) != 1:
			raise NameError, 'Ambiguous or non-existent Neuron id: '+str(id)
		else:
			name = name[0] # pull the string out of the list
		return name
	def name2id(self, name):
		try:
			id = name[name.rindex('_t')+2::] # everything from just after the last '_t' to the end of the neuron name, index() raises ValueError if it can't be found
		except ValueError:
			raise ValueError, 'Badly formatted Neuron name: '+name
		try:
			id = int(id) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return id
	def load(self): # or loadspikes()?
		f = file(self.path + self.name + '.spk', 'rb') # open the spike file for reading in binary mode
		self.spikes = numpy.fromfile(f, dtype=numpy.int64) # read in all spike times in us
		f.close()
		#treestr = self.level*TAB + self.name + '/'
		#self.writetree(treestr+'\n'); print treestr # print string to tree hierarchy and screen
	def cut(self, tstart=None, tend=None): # maybe use a masked array instead
		'''Returns all of the Neuron's spike times that fall within tstart and tend'''
		#return self.spikes[ numpy.greater_equal(self.spikes, tstart) & numpy.less_equal(self.spikes, tend) ]
		if tstart is None:
			tstart = self.spikes[0]
		if tend is None:
			tend = self.spikes[-1]
		return self.spikes[ (self.spikes >= tstart) & (self.spikes <= tend) ]
	def cutrel(self, tstart=None, tend=None):
		'''Cuts Neuron spike times and returns them relative to tstart'''
		if tstart is None:
			tstart = self.spikes[0]
		if tend is None:
			tend = self.spikes[-1]
		if tstart.__class__ is types.FloatType:
			warnings.warn('Converting tstart to int: '+str(int(tstart)))
		return self.cut(tstart, tend) - int(tstart)
	def rate(self, method='bin'):
		pass
	def raster(self):
		pass

class Movie(object): # careful with potential name conflict with Movies() in Dimstim
	'''A Movie stimulus object'''
	def __init__(self, name=DEFAULTMOVIENAME, path=DEFAULTMOVIEPATH, parent=None):
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
			raise RuntimeError('There are unread bytes in movie file \'%s\'. Width, height, or nframes is incorrect in the movie file header.' %self.name)
		#self.data = self.data[::,::-1,::] # flip the movie frames vertically for OpenGL's bottom left origin
		f.close() # close the movie file

################

# init some typical movies (but don't load 'em til needed), then just point to them within the appropriate Experiments
mseq32 = Movie(name='mseq32.m', parent=None)
mseq16 = Movie(name='mseq16.m', parent=None)
# shouldn't use sparse bar movies anymore, can access VisionEgg directly now, get the framebuffers to directly do STA
#sparsebars = Movie(path='C:/data/Cat 15/Track 7c/72 - track 7c sparseexps/', name='72 - track 7c sparseexps.sparsebars.movie');

# init and load a Recording to play around with
print 'Initing and loading Recording(92):'
r=Recording(92)
r.load()
