"""Defines the Movie class. This is really just a place to hold movie data,
don't confuse with dimstim.Movie subclass of dimstim.Experiment"""

#print 'importing Movie'

import dimstim
from Core import *

class Movie(object): # DON'T inherit from dimstim.Movie class!
    def __init__(self, fname=None, name=None, parent=None):
        """Movies don't need parents, they can just exist on their own and be used by anyone"""
        super(Movie, self).__init__()
        self.level = 5 # level in the hierarchy
        self.e = parent # save parent Experiment object
        self.fname = fname
        self.name = name
        if self.name == None and self.fname != None:
            self.name = os.path.split(self.fname)[-1] # get the name from the fname
    def load(self):
        """Load movie data"""
        f = file(self.fname, 'rb') # open the movie file for reading in binary format
        headerstring = f.read(5)
        if headerstring == 'movie': # a header has been added to the start of the file
            self.ncellswide, = struct.unpack('H', f.read(2)) # 'H'== unsigned short int
            self.ncellshigh, = struct.unpack('H', f.read(2))
            self.nframes, = struct.unpack('H', f.read(2))
            if self.nframes == 0: # this was used in Cat 15 mseq movies to indicate 2**16 frames, shouldn't really worry about this, cuz we're using slightly modified mseq movies now that we don't have the extra frame at the end that the Cat 15 movies had (see comment in Experiment module), and therefore never have a need to indicate 2**16 frames
                self.nframes = 2**16
            self.offset = 11 # header is this long
        else: # there's no header at the start of the file, set the file pointer back to the beginning and use these hard coded values:
            f.seek(0)
            self.ncellswide = self.ncellshigh = 64
            self.nframes = 6000
            self.offset = 0 # header is this long
        # read in all of the frames
        self.data = np.fromfile(f, np.uint8, count=self.nframes*self.ncellshigh*self.ncellswide)
        self.data = self.data.reshape(self.nframes, self.ncellshigh, self.ncellswide)
        #self.data = numarray.fromfile(f, np.UInt8, (self.nframes,self.ncellshigh,self.ncellswide))
        leftover = f.read() # check if there are any leftover bytes in the file
        if leftover != '':
            pprint(leftover)
            print self.ncellswide, self.ncellshigh, self.nframes
            raise RuntimeError, 'There are unread bytes in movie file %r. Width, height, or nframes is incorrect in the movie file header.' % self.fname
        #self.data = self.data[::, ::-1, ::] # don't need to flip the movie frames vertically for OpenGL's bottom left origin
        f.close() # close the movie file

# init some typical movies (but don't load 'em til needed). Then, just point to them within the appropriate Experiments
MSEQ32 = Movie(fname=os.path.join(DEFAULTMOVIEPATH, 'mseq32.m'), parent=None)
MSEQ16 = Movie(fname=os.path.join(DEFAULTMOVIEPATH, 'mseq16.m'), parent=None)

# shouldn't use sparse bar movies anymore, can access VisionEgg directly now, get the framebuffers to directly do STA
#sparsebars = Movie(path='C:/data/Cat 15/Track 7c/72 - track 7c sparseexps/', name='72 - track 7c sparseexps.sparsebars.movie');
