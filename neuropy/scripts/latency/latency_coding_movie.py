"""Creates a latency coding movie from an image, to try and change perception
of the image by manipulating its latency-contrast profile"""

from __future__ import division

import Image
import numpy as np
import sys
import os
import struct

nframes = 4
moviepath = '/mov/latency/'

assert len(sys.argv) <= 2
if len(sys.argv) == 2:
    fname = sys.argv[1]
else:
    fname = 'lincoln_head_small.jpg'
print 'Processing image file', fname

im = Image.open(fname)
# could just convert to luminance, but warn instead for now
assert im.mode == 'L', 'image file should be luminance only (greyscale)'
imdata = np.asarray(im)
height, width = imdata.shape
maxval = imdata.max()
minval = imdata.min()
meanval = int(round(imdata.mean()))
span = maxval - minval
step = span / nframes
threshes = [minval + (i+1)*step for i in range(nframes)] # pixel threshold values
print 'max pixel value is', maxval
print 'min pixel value is', minval
print 'mean pixel value is', meanval
print 'span is', span
print 'step is', step

# init movie
frames = meanval * np.ones((nframes, height, width), dtype=np.uint8) # all pixels same as mean image pix val

for framei, thresh in enumerate(threshes):
    if framei == 0:
        condition = (imdata <= thresh)
    else:
        condition = (threshes[framei-1] < imdata) * (imdata <= thresh)
    xy = condition.nonzero() # tuple of x and y indices where condition is true
    frames[framei][xy] = imdata[xy]

# save to a movie file
outfname, ext = os.path.splitext(fname)
outfname = outfname.upper() # same as image name, no ext, all caps
print 'Writing to movie file', moviepath+outfname
f = file(moviepath+outfname, 'wb') # open a movie file for writing in binary format
# write movie header
f.write('movie')
f.write(struct.pack('H', width)) # 'H'== unsigned short int == 2 bytes on this PC
f.write(struct.pack('H', height))
f.write(struct.pack('H', nframes))
for frame in frames:
    f.write(frame)
f.close()
