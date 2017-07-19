"""Export movie data from NVS movie format to .avi file.
Run from within neuropy using `run -i scripts/write_movie_avconv.py`
This doesn't need to be run in neuropy at all though. Plain Python works fine.

------------------
NOTE: To export from .avi to .jpg and then back to .avi (e.g., for creating shorter clips),
you can do this completely outside of Python with avconv at the command line.

export an .avi to a sequence of (1-based numbered) .jpgs:

avconv -i MVI_1400.AVI -f image2 -vcodec copy 201-500/%00d.jpg

Then, move the frames you want keep (say, 201-500) to another folder. To then build a
new .avi using just the remaining subset of (1-based numbered) frames:

avconv -f image2 -r 60 -start_number 201 -i 201-500/%00d.jpg -vcodec copy MAS_1400_CLR.avi

Or, you can rename them starting from 0 or 1 using Thunar, then you don't need to use the
'start_number' arg.

Note that the above command keeps the original colour. Maybe it's possible to do some basic
image manipulation with avconv, like conversion to grayscale, rotation, or contrast inversion.
------------------
"""

import numpy as np
import Image
import subprocess as sp
import os
import shutil
import struct

CONTRASTINVERT = False
REVERSE = False # in time
SHUFFLE = 'pink' # 'white': shuffle pixels of each frame; 'pink': shuffle spatial phase
invstr, revstr, shfstr = '', '', ''
if CONTRASTINVERT:
    invstr = '_INV'
if REVERSE:
    revstr = '_REV'
if SHUFFLE == 'white':
    shfstr = '_WHT'
elif SHUFFLE == 'pink':
    shfstr = '_PNK'
FPS = 60 # set the frame rate in the output .avi
SCALESPACE = 1 #16 # resize the movie by this factor in both x and y
SCALETIME = 1 #6 # repeat each frame this many times
# extracted jpg file quality, from 1 to 100, or None for PIL's default, only really necessary
# for mseq:
JPGQUALITY = 100

# choose your movie and frame range:

mvifname = 'MVI_1400'
path = os.path.expanduser('~/data/slab/mov/2007-11-24')
framei0, framei1  = 200, 500 # aka, MAS_1400

'''
mvifname = 'MVI_1400'
path = os.path.expanduser('~/data/slab/mov/2007-11-24')
framei0, framei1  = 3300, 3600 # aka, MAS_1400_B
'''
'''
mvifname = 'MVI_1403'
path = os.path.expanduser('~/data/slab/mov/2007-11-24')
framei0, framei1 = 0, 300 # aka, MAS_1403
'''
'''
mvifname = 'MVI_1419'
path = os.path.expanduser('~/data/slab/mov/2007-11-25')
framei0, framei1 = 3000, 3300
'''
'''
mvifname = 'MSEQ16'
path = os.path.expanduser('~/data/slab/mov/mseq')
framei0, framei1 = 0, 16383
'''


class Movie(object):
    def __init__(self, fname):
        """Load NVS movie data"""
        f = file(fname, 'rb') # open the movie file for reading in binary format
        self.f = f
        headerstring = f.read(5)
        assert headerstring == 'movie'
        self.ncellswide, = struct.unpack('H', f.read(2)) # 'H'== unsigned short int
        self.ncellshigh, = struct.unpack('H', f.read(2))
        self.nframes, = struct.unpack('H', f.read(2))
        print('width, height, nframes:', self.ncellswide, self.ncellshigh, self.nframes)
        self.framesize = self.ncellshigh*self.ncellswide
        self.frames = np.fromfile(f, dtype=np.uint8, count=self.nframes*self.framesize)
        self.frames.shape = self.nframes, self.ncellshigh, self.ncellswide
        leftover = f.read() # check if there are any leftover bytes in the file
        if leftover != '':
            #pprint(leftover)
            raise RuntimeError('There are unread bytes in movie file %r. Width, height, or'
                               'nframes are incorrect in the movie file header.' % fname)

def shuffle_phase_separate_halves(mvi):
    """Shuffle phases of movie while maintaining amplitudes. For a natural scene
    movie, this should result in pink noise"""
    # make sure there are an even number of pixels, for indexing the FFT components:
    n = mvi.size
    assert n % 2 == 0
    MVI = np.fft.fftn(mvi) # complex numbers, same shape as mvi
    MVI = MVI.ravel() # flatten temporarily
    DC = MVI[0] # DC term
    #PMVI = MVI[1:n/2] # increasingly positive freq terms, excluding DC and Nyquist terms
    #NQ = MVI[n/2] # Nyquist term
    #NMVI = MVI[n/2+1:] # decreasingly negative freq terms, excluding DC and Nyquist terms
    MVI = MVI[1:n/2+1] # PMVI, NQ
    r = np.abs(MVI) # amplitudes
    #nr = np.abs(NMVI) # negative freq amplitudes
    # get uniform random distribution of angles in radians from -pi to pi
    theta = np.pi * (np.random.random(MVI.shape) * 2 - 1)
    #ntheta = ptheta[::-1]
    # convert polar back to rect:
    #PPNK = r * np.cos(ptheta) +  1j * r * np.sin(ptheta)
    PNK = r * np.exp(1j*theta) # this is equivalent to above, but slightly faster
    NPNK = np.conj(PNK[:-1][::-1])
    #NPNK = nr * np.exp(-1j*ntheta)
    #NPNK = np.conj(PPNK[::-1])
    PNK = np.concatenate([[DC], PNK, NPNK])
    PNK.shape = mvi.shape # restore shape
    pnk = np.fft.ifftn(PNK)
    pnk = np.real(pnk)
    # normalize back to exactly the same pixel value range as input mvi:
    pnk = pnk - pnk.min()
    pnk = pnk / pnk.max() # now ranges from 0 to 1
    pnk = pnk * mvi.ptp() # same dynamic range as mvi
    pnk = pnk + mvi.min() # same offset as mvi
    pnk = np.int64(np.round(pnk))
    assert pnk.min() == mvi.min()
    assert pnk.max() == mvi.max()
    assert pnk.min() >= 0
    assert pnk.max() <= 255
    mvi = np.uint8(pnk)
    return mvi

def shuffle_phase(mvi):
    """Shuffle phases of movie while maintaining amplitudes. For a natural scene
    movie, this should result in pink noise"""
    ## TODO: confirm that output spectrogram is very similar to input spectrogram
    ## TODO: generate pink noise movie explicitly, see how it compares to phase-shuffled
    ## natscene movie
    # make sure there are an even number of pixels, for indexing the FFT components:
    n = mvi.size
    assert n % 2 == 0
    MVI = np.fft.fftn(mvi) # complex numbers, same shape as mvi
    #MVI = MVI.ravel() # flatten temporarily
    r = np.abs(MVI) # amplitudes
    # get uniform random distribution of angles in radians from -pi to pi
    theta = np.pi * (np.random.random(MVI.shape) * 2 - 1)
    # convert polar back to rect:
    #PPNK = r * np.cos(ptheta) +  1j * r * np.sin(ptheta)
    PNK = r * np.exp(1j*theta) # this is equivalent to above, but slightly faster
    #PNK.shape = mvi.shape # restore shape
    pnk = np.fft.ifftn(PNK)
    pnk = np.real(pnk)
    #pnk = np.abs(pnk)
    # normalize back to exactly the same pixel value range as input mvi:
    pnk = pnk - pnk.min()
    pnk = pnk / pnk.max() # now ranges from 0 to 1
    pnk = pnk * mvi.ptp() # same dynamic range as mvi
    pnk = pnk + mvi.min() # same offset as mvi
    pnk = np.int64(np.round(pnk))
    assert pnk.min() == mvi.min()
    assert pnk.max() == mvi.max()
    assert pnk.min() >= 0
    assert pnk.max() <= 255
    mvi = np.uint8(pnk)
    return mvi

fname = os.path.join(path, mvifname)
m = Movie(fname)
mvi = m.frames[framei0:framei1] # 3D numpy array
assert mvi.dtype == np.uint8

if CONTRASTINVERT:
    mvi = 255 - mvi # invert contrast of all pixels, assumes 8 bit pixels

if REVERSE:
    mvi = mvi[::-1] # reverse frame order

if SHUFFLE == 'white':
    for framei, frame in enumerate(mvi):
        np.random.shuffle(frame.ravel())
        mvi[framei] = frame
elif SHUFFLE == 'pink':
    mvi = shuffle_phase(mvi)

if SCALESPACE > 1: # scale it up in both spatial dimensions
    mvi = np.repeat(np.repeat(mvi, SCALESPACE, axis=1), SCALESPACE, axis=2)

if SCALETIME > 1: # scale it up in time (number of frames)
    mvi = np.repeat(mvi, SCALETIME, axis=0)

basename = mvifname + '_' + str(framei0) + '-' + str(framei1) # e.g. MVI_1403_0-300
fullname = ('%s%s%s%s_%sfps_%sxy_%st' % (os.path.join(path, basename),
            invstr, revstr, shfstr, FPS, SCALESPACE, SCALETIME))
if JPGQUALITY != None:
    fullname += '_Q%d' % JPGQUALITY
fnameavi = fullname + '.avi'
framespath = fullname + '_frames'
try:
    os.mkdir(framespath)
except OSError:
    pass # frames folder already exists

print('writing frames to %s' % framespath)
for i, frame in enumerate(mvi):
    im = Image.fromarray(frame)
    # save a sequence of .jpg files to disk
    if JPGQUALITY == None: # use PIL's default
        im.save(os.path.join(framespath, "%d.jpg" % i))
    else:
        im.save(os.path.join(framespath, "%d.jpg" % i), format='JPEG', subsampling=0,
                quality=JPGQUALITY)

# convert .jpg files to .avi using external program:
FFMPEG_BIN = 'avconv'
command = [ FFMPEG_BIN,
            #'-y', # overwrite output file if it exists, otherwise prompt
            '-f', 'image2',
            '-r', '%d' % FPS, # frames per second
            '-i', os.path.join(framespath, r'%00d.jpg'), # input
            '-vcodec', 'copy',
            #'-qscale',  '1',
            #'-vcodec', 'rawvideo',
            #'-vcodec', 'mjpeg',
            #'-vcodec', 'mpeg4',
            #'-s', '320x240', # size of one frame
            #'-pix_fmt', 'gray',
            fnameavi ]
sp.call(command)
print('saved .avi movie to %s' % fnameavi)

#shutil.rmtree(framespath) # recursive delete
#print('removed %s' % framespath)


'''
# Export movie data from neuropy:

fnamenpy = os.path.join(path, basename) + '.npy'
np.save(fnamenpy, mvi)
print('saved .npy movie to %s' % fnamenpy)

# load movie data exported from neuropy:
mvi = np.load(fnamenpy)
os.remove(fnamenpy)
print('removed %s' % fnamenpy)
'''


"""
To potentially skip the intermediate .jpg files written to disk, perhaps this method using pipes would work, from http://stackoverflow.com/a/13298538:

Ok I got it working. thanks to LordNeckbeard suggestion to use image2pipe. I had to use jpg encoding instead of png because image2pipe with png doesn't work on my verision of ffmpeg. The first script is essentially the same as your question's code except I implemented a simple image creation that just creates images going from black to red. I also added some code to time the execution.

serial execution

import subprocess, Image

fps, duration = 24, 100
for i in range(fps * duration):
    im = Image.new("RGB", (300, 300), (i, 1, 1))
    im.save("%07d.jpg" % i)
subprocess.call(["ffmpeg","-y","-r",str(fps),"-i", "%07d.jpg","-vcodec","mpeg4", "-qscale","5", "-r", str(fps), "video.avi"])
parallel execution (with no images saved to disk)

import Image
from subprocess import Popen, PIPE

fps, duration = 24, 100
p = Popen(['ffmpeg', '-y', '-f', 'image2pipe', '-vcodec', 'mjpeg', '-r', '24', '-i', '-', '-vcodec', 'mpeg4', '-qscale', '5', '-r', '24', 'video.avi'], stdin=PIPE)
for i in range(fps * duration):
    im = Image.new("RGB", (300, 300), (i, 1, 1))
    im.save(p.stdin, 'JPEG')
p.stdin.close()
p.wait()
the results are interesting, I ran each script 3 times to compare performance: serial:

12.9062321186
12.8965060711
12.9360799789
parallel:

8.67797684669
8.57139396667
8.38926696777
So it seems the parallel version is faster about 1.5 times faster.

"""
