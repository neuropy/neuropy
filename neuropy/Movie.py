"""Defines the Movie class. This is really just a place to hold movie data,
don't confuse with dimstim.Movie subclass of dimstim.Experiment

FIXME: what's the point of this class now? Why not just let dimstim init a Movie experiment?
    - can understand that neuropy.Movie is needed for dimstim < 0.16
"""

#print 'importing Movie'

import dimstim
from Core import *
from Core import _data # ensure it's imported, in spite of leading _

class Movie(object): # DON'T inherit from dimstim.Movie class!
    def __init__(self, fname=None, name=None, parent=None):
        """Movies don't need parents, they can just exist on their own and be used by anyone"""
        self.level = 5 # level in the hierarchy
        self.fname = fname
        self.name = name
        self.parent = parent # save parent object, might be an Experiment, might not
        if self.name == None and self.fname != None:
            self.path, self.fname = os.path.split(self.fname) # separate path from fname
            self.name = os.path.splitext(self.fname)[0] # extentionless fname
            if self.name not in _data.movies:
                _data.movies[self.name] = self # add self to _data.movies dictattr
        else:
            pass # both self.name and self.fname are None, this happens when executing Cat 15 textheaders, where you init a movie with m = Movie(), and only later assign its fname field. In this case, the .loadprecat16exp() method handles adding movies init'd from textheader to the _data.movies dictattr

    def load(self):
        """Load movie data"""
        try:
            self.data # movie's already been loaded, don't do anything
        except AttributeError:
            try:
                self.data = _data.movies[self.name].data # if a Movie init'd with the same name already has its data loaded, use it
            except AttributeError:
                f = file(os.path.join(self.path, self.fname), 'rb') # open the movie file for reading in binary format
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
                treestr = self.level*TAB + os.path.join(self.path, self.fname)
                print treestr
