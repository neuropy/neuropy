"""Plot x and y positions of all spikes from each neuron in specified tracks,
as a function of time. Copy and paste into neuropy console"""

## TODO: make figure width proportional to track duration

from colour import CLUSTERCOLOURDICT # for plotting on white

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]

spikefnames = ['/home/mspacek/data/ptc15/tr7c/track7c.track_2012-08-07_08.59.10.spike',
               '/home/mspacek/data/ptc22/tr1/track1.track_2012-03-02_15.56.59.spike',
               '/home/mspacek/data/ptc22/tr2/track2.track_2012-10-04_11.00.00.spike']

for track, spikefname in zip(tracks, spikefnames):
    titlestr = track.absname
    spikes = np.load(spikefname) # record array
    nids, ts = sorted(track.alln), spikes['t']
    x0s, y0s = spikes['x0'], spikes['y0']
    nid2txyc = {}
    # collect t, x, y and c for all nids:
    for nid in nids:
        sids, = np.where(spikes['nid'] == nid)
        t, x, y = ts[sids], x0s[sids], y0s[sids]
        t = t / 1e6 / 3600 # convert from us to hours
        c = CLUSTERCOLOURDICT[nid]
        nid2txyc[nid] = t, x, y, c
    fx = figure()
    ax = gca()
    fy = figure()
    ay = gca()
    # plot data for each nid, one at a time:
    for nid in nids:
        t, x, y, c = nid2txyc[nid]
        ax.plot(t, x, '.', ms=1, c=c)
        ay.plot(t, y, '.', ms=1, c=c)
    ax.set_ylim((-100, 100)) # um
    ay.invert_yaxis()
    ax.set_xlabel('time (hours)')
    ax.set_ylabel('x position ($\mu$m)')
    ay.set_xlabel('time (hours)')
    ay.set_ylabel('y position ($\mu$m)')
    ax.set_title(titlestr)
    ay.set_title(titlestr)
    fx.canvas.manager.set_window_title(titlestr+'_spike_xpos')
    fy.canvas.manager.set_window_title(titlestr+'_spike_ypos')
    fx.tight_layout(pad=0.3) # resize contents to figure
    fy.tight_layout(pad=0.3) # resize contents to figure
