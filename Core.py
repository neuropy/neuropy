'''Core neuropy functions and classes'''

DEFAULTDATAPATH = 'C:/data/' # the convention in neuropy will be that all 'path' var names have a trailing slash
DEFAULTCATID    = 15
DEFAULTTRACKID  = '7c'
DEFAULTRIPNAME  = 'liberal spikes' # a rip is the name of a spike sorting extraction; rips with the same name should have been done with the same templates, and possibly the same ripping thresholds
SLASH = '/' # use forward slashes instead of having to use double backslashes

import os, types

# rename all the class args (like catid, catname) to be simply id and name, since the class itself already defines what kind of id or name it is
# worry about conversion of ids to strings: some may be only 1 digit and may have a leading zero!
# maybe make two load() f'ns for Experiment and Neuron: one from files, and a future one from a database
# make a save() f'n that pickles the object (including any of its results, like its STA, tuning curve points, etc)?
# maybe rename to fopen()? what's the standard f'n for opening a db?

"""
def str2(data):
	if type(data) is types.IntTypes:
		s = str(data)
		if len(s) == 1:
			s = '0'+s # add a leading zero for single digits
"""

class Data(object): # use 'new-style' classes
	'''Data can have multiple Cats'''
	def __init__(self, dataPath=DEFAULTDATAPATH):
		self.path = dataPath
		print 'Data() has been init\'d!'
	def load(self):
		self.c = {} # store Cats in a dictionary
		catNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Cat ') ] # os.listdir() returns all dirs AND files
		for catName in catNames:
			cat = Cat(catid=None, catName=catName, parent=self) # make an instance using just the catName (let it figure out the catid)
			cat.load() # load the Cat
			self.c[cat.id] = cat # save it

class Cat(object):
	'''A Cat can have multiple Tracks'''
	def __init__(self, catid=DEFAULTCATID, catName=None, parent=Data):
		try:
			self.d = parent() # init parent Data object
		except TypeError: # parent is an instance, not a class
			self.d = parent # save parent Data object
		if catid is not None:
			catName = self.id2name(self.d.path,catid) # use the id to get the name
		elif catName is not None:
			catid = self.name2id(catName) # use the name to get the id
		else:
			raise ValueError, 'catid and catName can\'t both be None'
		self.id = catid
		self.name = catName
		self.path = self.d.path + self.name + SLASH
	def id2name(self, path, catid):
		catName = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith('Cat '+str(catid)) ]
		if len(catName) != 1:
			raise NameError, 'Ambiguous or non-existent Cat id: '+str(catid)
		else:
			catName = catName[0] # pull the string out of the list
		return catName
	def name2id(self, catName):
		catid = catName.replace('Cat ','',1) # replace first occurrence of 'Cat ' with nothing
		if not catid:
			raise NameError, 'Badly formatted Cat name: '+catName
		try:
			catid = int(catid) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return catid
	def load(self):
		self.t = {} # store Tracks in a dictionary
		trackNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.startswith('Track ') ]
		for trackName in trackNames:
			track = Track(trackid=None, trackName=trackName, parent=self) # make an instance using just the trackName (let it figure out the trackid)
			track.load() # load the Track
			self.t[track.id] = track # save it

class Track(object):
	'''A Track can have multiple Recordings'''
	def __init__(self, trackid=DEFAULTTRACKID, trackName=None, parent=Cat):
		try:
			self.c = parent() # init parent Cat object
		except TypeError: # parent is an instance, not a class
			self.c = parent # save parent Cat object
		if trackid is not None:
			trackName = self.id2name(self.c.path,trackid) # use the id to get the name
		elif trackName is not None:
			trackid = self.name2id(trackName) # use the name to get the id
		else:
			raise ValueError, 'trackid and trackName can\'t both be None'
		self.id = trackid
		self.name = trackName
		self.path = self.c.path + self.name + SLASH
	def id2name(self, path, trackid):
		trackName = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith('Track '+str(trackid)) ]
		if len(trackName) != 1:
			raise NameError, 'Ambiguous or non-existent Track id: '+str(trackid)
		else:
			trackName = trackName[0] # pull the string out of the list
		return trackName
	def name2id(self, trackName):
		trackid = trackName.replace('Track ','',1) # replace first occurrence of 'Track ' with nothing
		if not trackid:
			raise NameError, 'Badly formatted Track name: '+trackName
		try:
			trackid = int(trackid) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return trackid
	def load(self):
		self.r = {} # store Recordings in a dictionary
		recordingNames = [ dirname for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname[0:2].isdigit() and dirname.count(' - ') == 1 ] # 1st 2 chars in dirname must be digits, must contain exactly 1 occurrence of ' - '
		for recordingName in recordingNames:
			recording = Recording(recordingid=None, recordingName=recordingName, parent=self) # make an instance using just the recordingName (let it figure out the recordingid)
			recording.load() # load the Recording
			self.r[recording.id] = recording # save it

class Recording(object):
	'''A Recording corresponds to a single SURF file, ie everything recorded between when
	the user hits record and when the user hits stop and closes the SURF file, including any
	pauses in between Experiments within that Recording. A Recording can have multiple Experiments,
	and multiple spike extractions, called Rips'''
	def __init__(self, recordingid=None, recordingName=None, parent=Track):
		try:
			self.t = parent() # init parent Track object
		except TypeError: # parent is an instance, not a class
			self.t = parent # save parent Track object
		if recordingid is not None:
			recordingName = self.id2name(self.t.path,recordingid) # use the id to get the name
		elif recordingName is not None:
			recordingid = self.name2id(recordingName) # use the name to get the id
		else:
			raise ValueError, 'recordingid and recordingName can\'t both be None'
		self.id = recordingid
		self.name = recordingName
		self.path = self.t.path + self.name + SLASH
		self.defaultripname = DEFAULTRIPNAME
	def id2name(self, path, recordingid):
		recordingName = [ dirname for dirname in os.listdir(path) if os.path.isdir(path+dirname) and dirname.startswith(str(recordingid)+' - ') ]
		if len(recordingName) != 1:
			raise NameError, 'Ambiguous or non-existent Recording id: '+str(recordingid)
		else:
			recordingName = recordingName[0] # pull the string out of the list
		return recordingName
	def name2id(self, recordingName):
		try:
			recordingid = recordingName[0:recordingName.index(' - ')] # everything before the first ' - ', index() raises ValueError if it can't be found
		except ValueError:
			raise ValueError, 'Badly formatted Recording name: '+recordingName
		try:
			recordingid = int(recordingid) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return recordingid
	def load(self):
		self.e = {} # store Experiments in a dictionary
		experimentNames = [ fname[0:fname.rfind('.din')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.din') ] # returns din filenames without their .din extension
		print experimentNames
		for (experimentid, experimentName) in enumerate(experimentNames): # experimentids will be according to alphabetical order of experimentNames
			experiment = Experiment(experimentid=experimentid, experimentName=experimentName, parent=self) # pass both the id and the name
			experiment.load() # load the Experiment
			self.e[experiment.id] = experiment # save it
		self.rip = {} # store Rips in a dictionary
		ripNames = [ dirname[0:dirname.rfind('.rip')] for dirname in os.listdir(self.path) if os.path.isdir(self.path+dirname) and dirname.endswith('.rip') ] # returns rip folder names without their .rip extension
		for ripName in ripNames:
			rip = Rip(ripName=ripName, parent=self) # pass just the name, ain't no such thing as a ripid, at least for now
			rip.load() # load the Rip
			self.rip[rip.name] = rip # save it
		try: # make the Neurons from the default rip (if it exists in the Recording path) available in the Recording, so you can access them via r.n[nid] instead of having to do r.rip[ripName].n[nid]. Make them just another pointer to the data in r.rip[defaultripname].n
			self.n = self.rip[self.defaultripname].n
			self.defaultrippath = self.path + self.defaultripname + SLASH
		except:
			pass

class Experiment(object):
	'''An Experiment corresponds to a single contiguous VisionEgg stimulus session.
	It contains information about the stimulus during that session, including
	the DIN values and the text header'''
	def __init__(self, experimentid=None, experimentName=None, parent=Recording): # Experiment IDs are 1-based in the .din filenames, at least for now. They should be renamed to 0-based. Here, they're treated as 0-based
		try:
			self.r = parent() # init parent Recording object
		except TypeError: # parent is an instance, not a class
			self.r = parent # save parent Recording object
		if experimentName is None:
			raise ValueError, 'experimentName can\'t be None'
		self.id = experimentid # not really used by the Experiment class, just there for user's info
		self.name = experimentName
		self.path = self.r.path
	# doesn't need a id2name or name2id method, neither can really be derived from the other in an easy, the id is just alphabetical order
	def load(self):
		self.din = [] # a list of tuples?
		self.textheader = ''
		# then, for each line in the textheader, exec() it so you get self.varname saved directly within in the Experiment object - watch out, will try and make Movie() objects and Bar() objects, etc?

class Rip(object):
	def __init__(self, ripName=DEFAULTRIPNAME, parent=Recording):
		try:
			self.r = parent() # init parent Recording object
		except TypeError: # parent is an instance, not a class
			self.r = parent # save parent Recording object
		if ripName is None:
			raise ValueError, 'ripName can\'t be None'
		# rips don't have ids, at least for now. Just names
		self.name = ripName
		self.path = self.r.path + self.name + '.rip' + SLASH # have to add .rip extension to ripName to get its actual folder name
	def load(self):
		self.n = {} # store Neurons in a dictionary
		neuronNames = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(self.path) if os.path.isfile(self.path+fname) and fname.endswith('.spk') ] # returns spike filenames without their .spk extension
		for neuronName in neuronNames:
			neuron = Neuron(neuronid=None, neuronName=neuronName, parent=self) # make an instance using just the neuronName (let it figure out the neuronid)
			neuron.load() # load the neuron
			self.n[neuron.id] = neuron # save it
		# then, maybe add something that loads info about the rip, say from some file describing the template used, and all the thresholds, exported to the same folder by SURF
		# maybe also load the template used for the rip, perhaps also stored in the same folder

class Neuron(object):
	'''A Neuron object's spike data spans all the Experiments within a Recording.
	If different Recordings have Rips with the same name, you can assume that the
	same spike template was used for all of them, and that therefore the neuronids
	are the same'''
	def __init__(self, neuronid=None, neuronName=None, parent=Rip): # neuron names don't include the '.spk' ending, although neuron filenames do
		try:
			self.rip = parent() # init parent Rip object
		except TypeError: # parent is an instance, not a class
			self.rip = parent # save parent Rip object
		if neuronid is not None:
			neuronName = self.id2name(self.rip.path,neuronid) # use the id to get the name
		elif neuronName is not None:
			neuronid = self.name2id(neuronName) # use the name to get the id
		else:
			raise ValueError, 'neuronid and neuronName can\'t both be None'
		self.id = neuronid
		self.name = neuronName
		self.path = self.rip.path
	def id2name(self, path, neuronid):
		neuronName = [ fname[0:fname.rfind('.spk')] for fname in os.listdir(path) if os.path.isfile(path+fname) and \
		             ( fname.find('_t'+str(neuronid)+'.spk')!=-1 or fname.find('_t0'+str(neuronid)+'.spk')!=-1 or fname.find('_t00'+str(neuronid)+'.spk')!=-1 ) ] # have to deal with leading zero ids, go up to 3 digit ids, should really use a re to do this properly...
		if len(neuronName) != 1:
			raise NameError, 'Ambiguous or non-existent Neuron id: '+str(neuronid)
		else:
			neuronName = neuronName[0] # pull the string out of the list
		return neuronName
	def name2id(self, neuronName):
		try:
			neuronid = neuronName[neuronName.rindex('_t')+2::] # everything from just after the last '_t' to the end of the neuron name, index() raises ValueError if it can't be found
		except ValueError:
			raise ValueError, 'Badly formatted Neuron name: '+neuronName
		try:
			neuronid = int(neuronid) # convert string to int if possible
		except ValueError:
			pass # it's alphanumeric, leave it as a string
		return neuronid
	def load(self): # or loadspikes()?
		self.spikes = [] # use numpy array instead of a list? faster?

"""
class Stim(Experiment):
	'''A Stim contains the visual stimulus information for a single Experiment. A Stim corresponds
	to a single contiguous VisionEgg stimulus session'''
	def __init__(self,fname):
		pass
	def load(self,fname):
		pass
"""
