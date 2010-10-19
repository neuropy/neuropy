"""Script for generating mpegs from raw Movie data"""

import sys
from neuropy.Movie import Movie
from neuropy.core import enlarge

#import pymedia
#import pymedia.video
import pymedia.video.vcodec as vcodec
import numpy as np

#out = 'C:/home/mspacek/Desktop/tracking.mpeg'
out = 'C:/home/mspacek/Desktop/reliability.mpeg'
enlargex = 2 # factor to enlarge each frame by, set to 1 to not enlarge

#m = Movie(name='ns1-64p-50h-2m-mnNAT-ctNAT', path='C:/pub/Movies/tracking/', parent=None)
m = Movie(name='nre-64p-50h-2m-mn125-ct035-single.m', path='C:/pub/Movies/reliability/e/single/', parent=None)
m.load()

#width = 800
#height = 600
width = m.ncellswide * enlargex
height = m.ncellshigh * enlargex
nruns = 10
iri = 50 # inter-run interval, in frames, during which to have greyscreen

def prepframe(frame, enlargex=1):
    """Takes raw luminance data frame and converts it to RGB yuvFrame suitable for encoding with pymedia"""
    height, width = np.asarray(frame.shape)*enlargex
    frame = frame.T
    frame = enlarge(frame, enlargex)
    frame = np.asarray([frame, frame, frame]).T
    bmpFrame = vcodec.VFrame( vcodec.formats.PIX_FMT_RGB24, (width, height), (frame.tostring(), None, None) )
    yuvFrame = bmpFrame.convert(vcodec.formats.PIX_FMT_YUV420P)
    return yuvFrame

blankframe = np.zeros((height, width), dtype=np.int8) + 127
blankyuvFrame = prepframe(blankframe, enlargex=enlargex)

f = file(out, 'wb')

params= { 'type': 0,
          'gop_size': 12,
          'max_b_frames': 0,
          'frame_rate_base': 1,
          'frame_rate': 50,
          'width': width,
          'height': height,
          'deinterlace': 0,
          'bitrate': 3000000, # normally 2700000 for a framerate of 2997/125 Hz
          'id': vcodec.getCodecID( 'mpeg1video' )
        }

e = vcodec.Encoder(params)

for runi in range(nruns):
    for frame in m.data:
        yuvFrame = prepframe(frame, enlargex=enlargex)
        d = e.encode(yuvFrame)
        f.write(d.data)
        sys.stdout.write('.')
    if iri:
        for framei in range(iri):
            d = e.encode(blankyuvFrame)
            f.write(d.data)
            sys.stdout.write('.')

f.close()
