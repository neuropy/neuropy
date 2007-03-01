#mpeg = 'C:/home/mspacek/Desktop/output.mpeg'
mp3name  = 'C:/data/Cat 15/Track 7c/76 - track 7c tracking/76 - track 7c tracking.2.ns1-64p-50h-2m-mnNAT-ctNAT.mp3'
enlargex = 2 # factor to enlarge each frame by, set to 1 to not enlarge

import pymedia.muxer as muxer
import pymedia.audio.acodec as acodec
import pymedia.video.vcodec as vcodec

from neuropy.Movie import Movie
from neuropy.Core import enlarge

m = Movie(name='ns1-64p-50h-2m-mnNAT-ctNAT', path='C:/pub/Movies/tracking/', parent=None)
m.load()

mp3file = file(mp3name, 'rb')
mp3data = mp3file.read() # read the whole file in
mp3file.close()

f = file('output.avi', 'wb') # output file

# decode a frame of the mp3 just to get its params, which are then used during muxing with the video
dm = muxer.Demuxer('mp3')
frames = dm.parse(mp3data) # demux the data into sub streams, in this case, there's only the audio stream
frame0 = frames[0] # get the first frame of the sub streams. This is a (stream id, data) tuple
dec = acodec.Decoder(dm.streams[ frame0[0] ]) # open a decoder
decodedframe0 = dec.decode(frame0[1]) # decode it
aparams = {'channels': decodedframe0.channels,
           'sample_rate': decodedframe0.sample_rate,
           'bitrate': decodedframe0.bitrate,
           'id': acodec.getCodecID('mp3')
          }



width = m.ncellswide * enlargex
height = m.ncellshigh * enlargex

vparams = {'type': 0,
           'gop_size': 12,
           'max_b_frames': 0,
           'frame_rate_base': 1,
           'frame_rate': 50,
           'width': width,
           'height': height,
           'deinterlace': 0,
           'bitrate': 300000, # normally 2700000 for a framerate of 2997/125 Hz
           'id': vcodec.getCodecID('mpeg1video')
          }

ve = vcodec.Encoder(vparams) # encodes raw video data into the desired format


mux = muxer.Muxer('avi') # create muxer, one of: ['avi', 'wmv', 'mov', 'mp4', 'wma', 'mp2', 'mp3', 'ac3', 'aac', 'flac', 'mpeg', 'mpg', 'mpeg', 'dat', 'vob', 'm1v', 'ogg']
vstreamid = mux.addStream(muxer.CODEC_TYPE_VIDEO, vparams) # add a video stream
#astreamid = mux.addStream(muxer.CODEC_TYPE_AUDIO, aparams) # add an audio stream
#import pdb; pdb.set_trace()
# write header to disk
header = mux.start()
print 'hello'
f.write(header)

# write mp3 to disk
#f.write(mp3data)
s = mux.write(astreamid, mp3data) # returns a stream
f.write(s) # the stream can then be written to disk



#s = mux.write(astreamid, encodedframe) # returns a stream
#f.write(s) # the stream can then be written to disk

for frame in m.data:
    frame = frame.T
    frame = enlarge(frame, enlargex)
    frame = np.asarray([frame, frame, frame]).T
    bmpFrame = vcodec.VFrame( vcodec.formats.PIX_FMT_RGB24, (width, height), (frame.tostring(), None, None) )
    yuvFrame = bmpFrame.convert(vcodec.formats.PIX_FMT_YUV420P)
    encodedframe = ve.encode(yuvFrame)
    s = mux.write(vstreamid, encodedframe) # returns a stream
    f.write(s) # the stream can then be written to disk

    #f.write(d.data)

f.close()



#create muxer
#mux = muxer.Muxer( 'avi' ) # or mpg or something
#s1= vc1.encode( vfr )
#m_data.append( ( v_id1, s1 ))

