"""Plot histogram of total cell active durations (between first and last spike) normalized
by track duration, for all cells. Run with 'run -i scripts/active_duration.py' or copy and
paste into neuropy console"""

from __future__ import division
from __future__ import print_function

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
edges = np.arange(0, 1.1, 0.1)

adfs = [] # active duration fractions
for track in tracks:
    # total track duration, different from track.dt which excludes recording gaps:
    totaltrackdt = track.trange[1] - track.trange[0]
    for n in track.alln.values():
        adfs.append(n.dt / totaltrackdt)

adfs = np.hstack(adfs)

figure(figsize=(2.5, 2.92))
count = hist(adfs, bins=edges, color='k')[0]
xlim(0, 1)
ylim(ymax=count.max()) # normalize to max
xticks([0, 0.5, 1])
yticks([0, count.max()])
xlabel('active duration')
ylabel('neuron count')
gcfm().set_window_title('active_duration')
gcf().tight_layout(pad=0.3) # crop figure to contents

show()
