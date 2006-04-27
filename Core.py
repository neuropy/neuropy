'''Core neuropy functions and classes'''

# maybe make two load() f'ns: one from files, and a future one from a database
# maybe rename to fopen()? what's the standard f'n for opening a db?

class Cat:
	'''A Cat can have multiple Tracks'''
	def __init__(self):
		pass
	def load(self):
		pass

class Track:
	'''A Track can have multiple Recordings'''
	def __init__(self):
		pass

class Recording:
	'''A Recording corresponds to a single SURF file, ie everything recorded between when the user hits record and when the user hits stop, including any pauses in between Experiments within that Recording. A Recording can have multiple Experiments'''

class Experiment:
	'''An Experiment can have multiple Neurons and Stims'''
	def __init__(self,folder='path/to/Experiment'):
		self.n = []; # neuron list
		self.s = []; # stim list
		self.load(folder) # can this be done here???????????????
	def load(self,folder):
		for neuron in neurons:
			self.n.append(neuron)
		for stim in stims:
			self.s.append(stim)

class Neuron:
	'''A neuron persists across all Experiments in all Recordings in a Track'''
	def __init__(self,fname):
		pass
	def load(self,fname):
		pass

class Stim:
	def __init__(self,fname):
		pass
	def load(self,fname):
		pass

