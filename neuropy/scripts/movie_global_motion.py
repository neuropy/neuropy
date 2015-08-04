"""Extract optic flow vector fields from natural scene movies, average them to calculate
global motion within each movie for each specified recording. Adapted from
opencv/samples/python2/opt_flow.py. Run from within neuropy using `run -i
scripts/movie_global_motion.py`"""

from __future__ import division, print_function

import cv2
import numpy as np

recs = [ptc17.tr2b.r58, ptc18.tr1.r38, ptc18.tr2c.r58, ptc22.tr1.r08, ptc22.tr1.r10,
        ptc22.tr4b.r49]

pyr_scale = 0.5
levels = 3
winsize = 15
iterations = 3
poly_n = 5
poly_sigma = 1.2
flags = 0

FIGSIZE = (8, 4)

# calculate optic flow vector field between neighbouring pairs of frames, average their
# magnitudes to get global motion:
motion = {}
for rec in recs:
    name = rec.absname
    print(name)
    motion[name] = [] # optic flow magnitudes, in deg/sec, one per frame interval
    e0 = rec.e0
    movie = e0.e
    movie.load() # load movie data for this recording, flips frames vertically by default
    degpermoviepix = e0.s.widthDeg / movie.ncellswide
    dt = e0.d.sweepSec # frame duration in seconds
    frames = np.asarray(e0.e.frames)
    frameis = np.asarray(e0.d.framei) # movie frame indices used by this recording
    frames = frames[frameis] # dereference
    frame0 = frames[0] # init
    for frame1 in frames[1:]:
        ## TODO: if interested in flow direction, double-check vertical order of movie frames
        ## vs. what's expected by cv2.calcOpticalFlowFarneback. Should frames be flipped?
        flow = cv2.calcOpticalFlowFarneback(frame0, frame1, pyr_scale, levels, winsize,
                                            iterations, poly_n, poly_sigma, flags)
        mag, ang = cv2.cartToPolar(flow[:, :, 0], flow[:, :, 1]) # mag is in pix/frame
        # average over entire vector flow field in space, convert from pix/frame to deg/sec:
        motion[name].append(mag.mean() * degpermoviepix / dt)
        frame0 = frame1 # update for next iteration

    motion[name] = np.asarray(motion[name])
    figure(figsize=FIGSIZE)
    plot(frameis[1:], motion[name], 'k-', lw=1.5)
    xlabel('frame index')
    ylabel('motion amplitude (deg/s)')
    text(0.99, 0.98, '%s' % os.path.basename(e0.s.fname), # movie file name
                     horizontalalignment='right', verticalalignment='top',
                     transform=gca().transAxes, color='k')
    gcfm().window.setWindowTitle('movie_global_motion_%s' % name)
    tight_layout(pad=0.3)

pl.show()
