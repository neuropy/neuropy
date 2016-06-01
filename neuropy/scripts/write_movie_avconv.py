"""Export movie data to .avi file.
Run from within neuropy using `run -i scripts/write_movie_avconv.py`"""

import numpy as np
import Image
import subprocess as sp
import os
import shutil

CONTRASTINVERT = False
REVERSE = True
invstr, revstr = '', ''
if CONTRASTINVERT:

path = os.path.expanduser('~/data/mov/2007-11-24')
    invstr = '_INV'
if REVERSE:
    revstr = '_REV'
'''
mvifname = 'MVI_1400'
framei0, framei1  = 200, 500
e = ptc17.tr2b.r58.e0.e
'''
mvifname = 'MVI_1403'
framei0, framei1 = 0, 300
e = ptc22.tr1.r08.e0.e

basename = mvifname + '_' + str(framei0) + '-' + str(framei1) # e.g. MVI_1403_0-300
fnameavi = os.path.join(path, basename) + invstr + revstr + '.avi'
e.load()
assert e.f.name == os.path.join(path, mvifname)
mvi = np.asarray(e.frames[framei0:framei1]) # a 3D numpy array
if CONTRASTINVERT:
    assert mvi.dtype == np.uint8
    mvi = 255 - mvi # invert contrast of all pixels, assumes 8 bit pixels
if REVERSE:
    mvi = mvi[::-1] # reverse frame order

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
mvi = mvi[:, ::-1, :] # invert vertically for PIL
nframes = len(mvi)
assert nframes <= 1000 # otherwise the file numbering scheme below needs more leading 0s

framespath = os.path.join(path, 'frames')
try:
    os.mkdir(framespath)
except OSError:
    pass # frames folder already exists

print('writing frames to %s' % framespath)
for i, frame in enumerate(mvi):
    im = Image.fromarray(frame)
    im.save(os.path.join(path, "frames/%03d.jpg" % i)) # save a sequence of .jpg files to disk

# convert .jpg files to .avi using external program:
FFMPEG_BIN = 'avconv'
command = [ FFMPEG_BIN,
            '-y', # (optional) overwrite output file if it exists
            '-f', 'image2',
            '-r', '60', # frames per second
            '-i', os.path.join(path, r'frames/%03d.jpg'), # input
            '-vcodec', 'copy',
            #'-vcodec', 'rawvideo',
            #'-vcodec', 'mjpeg',
            #'-vcodec', 'mpeg4',
            #'-s', '320x240', # size of one frame
            #'-pix_fmt', 'gray',
            fnameavi ]
sp.call(command)
print('saved .avi movie to %s' % fnameavi)

shutil.rmtree(framespath) # recursive delete
print('removed %s' % framespath)



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
