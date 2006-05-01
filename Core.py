'''Core neuropy functions and classes'''

DEFAULTDATAPATH  = 'C:/data' # without trailing slash
DEFAULTCATID   = '15'
DEFAULTTRACKID = '7c'
SLASH = '/'

import os

# maybe make two load() f'ns: one from files, and a future one from a database
# make a save() f'n that pickles the object (including any of its results, like its STA, tuning curve points, etc)?
# maybe rename to fopen()? what's the standard f'n for opening a db?
"""
def isRelativePath(path):
	if path.count(SLASH) != 0:
		raise Exception, path + ' is not a relative path'
	return path

def parseCatName(catName):
	'''return everything after 'Cat ' '''

def parseTrackName(trackName):
	'''return everything after 'Track ' '''
"""
def parseRecordingName(recordingName):
	'''return everything before the the first occurence of ' -' '''




class Data(object):
	'''Data can have multiple Cats'''
	def __init__(self, dataPath=DEFAULTDATAPATH):
		self.dataPath = dataPath
		self.path = self.dataPath + SLASH
	def load(self):
		self.c = {} # store Cats in a dictionary


class Cat(Data): # subclass of Data
	'''A Cat can have multiple Tracks'''
	def __init__(self, catid=DEFAULTCATID, catName=None):
		super(Cat, self).__init__() # run constructor from Cat's superclass
		if catid is None: # if it was intentionally set to None during the call
			catid = self.name2id(self.path,catName) # use the name instead to get the id
		if catName is None:
			catName = self.id2name(self.path,catid)
		self.id = catid
		self.name = catName
		self.path += self.name + SLASH
	def name2id(self, catName):
		catid = catName.replace('Cat ','',1) # replace first instance of 'Cat ' with nothing
		if not catid:
			raise Exception, 'Badly formatted Cat name: '+catName
		try:
			catid = int(catid) # convert string to int if possible
		except:
			pass # it's alphanumeric, leave it as a string
		return catid
	def id2name(self, path, catid):
		catName = [ dirname for dirname in os.listdir(path) if ( os.path.isdir(path+SLASH+dirname) and dirname.startswith('Cat '+str(catid)) ) ] # os.listdir() actually returns all dirs AND files
		if len(catName) != 1:
			raise Exception, 'Ambiguous or non-existent Cat id: '+str(catid)
		else:
			catName = catName[0] # pull the string out of the list
		return catName
	def load(self):
		self.t = {} # store Tracks in a dictionary
		trackNames = [dirname for dirname in os.dirlisting(self.path) if dirname.startswith('Track')]
		for trackName in trackNames:
			track = Track(trackid=None, trackName=trackName) # make an instance using just the trackname (let it figure out the trackid)
			track.load() # load the track
			self.t[track.id] = track

class Track(Cat): # subclass of Cat
	'''A Track can have multiple Recordings'''
	def __init__(self, trackid=DEFAULTTRACKID, trackName=None):
		super(Track, self).__init__() # run constructor from Track's superclass
		if trackid is None: # if it was intentionally set to None during the call
			trackid = self.name2id(self.path,trackName) # use the name instead to get the id
		if trackName is None:
			trackName = self.id2name(self.path,trackid)
		self.path += self.trackName + SLASH
		self.catid = self.id # save inherited id
		self.catName = self.name # save inherited name
		self.id = trackid # overwrite
		self.name = trackName # overwrite
		self.path += self.name + SLASH
	def load(self):
		self.id = parseTrackName(self.trackName)
		self.r = {} # store Recordings in a dictionary
		recordingNames = os.dirlisting(self.path)
		for recordingName in recordingNames:
			recid = parseRecordingName(recordingName)
			self.r[recid] = Recording(recordingName)
			self.r[recid].load()

class Recording(Track):
	'''A Recording corresponds to a single SURF file, ie everything recorded between when
	the user hits record and when the user hits stop and closes the SURF file, including any
	pauses in between Experiments within that Recording. A Recording can have multiple Experiments'''
	def __init__(self,recordingName):
		super(Recording, self).__init__() # run constructor from Recording's superclass
		self.recordingName = isRelativePath(recordingName)
		self.path += self.recordingName + SLASH
	def load(self):
		self.id = parseRecordingName(self.recordingName)
		self.e = {} # store Recordings in a dictionary
		recordingNames = [fname for fname in os.dirlisting(self.path) if fname.endswith('.din')]
		for recordingName in recordingNames:
			recid = parseRecordingName(recordingName)
			self.r[recid] = Recording(recordingName)
			self.r[recid].load()

class Experiment(Recording):
	'''An Experiment corresponds to a single contiguous VisionEgg stimulus session.
	An Experiment can have multiple Neurons, but only one Stim'''
	def __init__(self,expi): # Experiment IDs are 1-based in the data, at least for now
		super(Recording, self).__init__() # run constructor from Recording's superclass
		self.neurons = []; # neuron list

class Neuron(Recording):
	'''Although a Neuron's id is the same across all Experiments in all Recordings in a Track, due to
	the way SurfBawd exports spike data, a Neuron object's spike data only spans all the Experiments within a Recording'''
	def __init__(self,fname):
	 	pass
	def load(self,fname): # or loadspikes()?
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
