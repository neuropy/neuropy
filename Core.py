'''Core neuropy functions and classes'''

DEFAULTBASEPATH  = 'C:/data'
DEFAULTCATPATH   = 'Cat 15'
DEFAULTTRACKPATH = 'Track 7c'

# maybe make two load() f'ns: one from files, and a future one from a database
# make a save() f'n that pickles the object (including any of its results, like its STA, tuning curve points, etc)?
# maybe rename to fopen()? what's the standard f'n for opening a db?

def isRelativePath(path):
	if path.count('/') != 0:
		raise Exception, path + ' is not a relative path'
	return path


class Data(object):
	'''Data can have multiple Cats'''
	def __init__(self,basePath=DEFAULTBASEPATH):
		self.basePath = basePath
	def load(self):
		self.cat = {} # store Cats in a dictionary

class Cat(Data):
	'''A Cat can have multiple Tracks'''
	def __init__(self,catPath=DEFAULTCATPATH):
		super(Cat, self).__init__() # run constructor from Cat's superclass
		self.catPath = isRelativePath(catPath)
	def load(self):
		self.track = {} # store Tracks in a dictionary

class Track(Cat):
	'''A Track can have multiple Recordings'''
	def __init__(self,trackPath=DEFAULTTRACKPATH):
		super(Cat, self).__init__() # run constructor from Track's superclass
		self.trackPath = isRelativePath(trackPath)
	def load(self):
		self.recording = {} # store Recordings in a dictionary

class Recording(Track):
	'''A Recording corresponds to a single SURF file, ie everything recorded between when
	the user hits record and when the user hits stop and closes the SURF file, including any
	pauses in between Experiments within that Recording. A Recording can have multiple Experiments'''
	def __init__(self,recordingPath):
		super(Cat, self).__init__() # run constructor from Recording's superclass
		self.recordingPath = isRelativePath(recordingPath)

class Experiment(Recording):
	'''An Experiment corresponds to a single contiguous VisionEgg stimulus session.
	An Experiment can have multiple Neurons, but only one Stim'''
	def __init__(self,path='path/to/Experiment'):
		self.neurons = []; # neuron list
		self.s = Stim()
		self.load(path) # can this be done here???????????????
	def load(self,path):
		for neuron in neurons:
			self.n.append(neuron)

class Neuron(Recording):
	'''Although a Neuron's id is the same across all Experiments in all Recordings in a Track, due to
	the way SurfBawd exports spike data, a Neuron object's spike data only spans all the Experiments within a Recording'''
	def __init__(self,fname):
		pass
	def load(self,fname): # or loadspikes()?
		self.spikes = [] # use numpy array instead of a list? faster?

class Stim(Experiment):
	'''A Stim contains the visual stimulus information for a single Experiment. A Stim corresponds
	to a single contiguous VisionEgg stimulus session'''
	def __init__(self,fname):
		pass
	def load(self,fname):
		pass

