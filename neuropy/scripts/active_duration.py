"""Plot histogram of total cell active durations (between first and last spike) normalized
by track duration, for all cells. Run with 'run -i scripts/active_duration.py' or copy and
paste into neuropy console"""

from __future__ import division
from __future__ import print_function

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
edges = np.arange(0, 1.1, 0.1)

adfs = [] # active duration fractions
fsfs = [] # first spike fractions
lsfs = [] # last spike fractions
for track in tracks:
    # total track duration, different from track.dt which excludes recording gaps
    # the following potentially excludes periods before the start and after then end of the
    # first and last experiment:
    #totaltrackdt = track.trange[1] - track.trange[0]
    # define instead as time between first and last spike:
    fs = min([ n.spikes[0] for n in track.alln.values() ]) # first spike
    ls = max([ n.spikes[-1] for n in track.alln.values() ]) # last spike
    totaltrackdt = ls - fs
    for n in track.alln.values():
        adfs.append(n.dt / totaltrackdt)
        fsfs.append((n.spikes[0] - fs) / totaltrackdt)
        lsfs.append((n.spikes[-1] - fs) / totaltrackdt)

adfs = np.hstack(adfs)
fsfs = np.hstack(fsfs)
lsfs = np.hstack(lsfs)

figure(figsize=(2.5, 2.92))
count = hist(adfs, bins=edges, color='k')[0]
xlim(0, 1)
ylim(ymax=count.max()) # normalize to max
xticks([0, 0.5, 1])
yticks([0, count.max()])
xlabel('fractional active duration')
ylabel('neuron count')
gcfm().set_window_title('active_duration')
gcf().tight_layout(pad=0.3) # crop figure to contents

figure(figsize=(2.5, 2.92))
count = hist(fsfs, bins=edges, color='k')[0]
xlim(0, 1)
ylim(ymax=count.max()) # normalize to max
xticks([0, 0.5, 1])
yticks([0, count.max()])
xlabel('first spike time fractions')
ylabel('neuron count')
gcfm().set_window_title('first_spike_fractions')
gcf().tight_layout(pad=0.3) # crop figure to contents

figure(figsize=(2.5, 2.92))
count = hist(lsfs, bins=edges, color='k')[0]
xlim(0, 1)
ylim(ymax=count.max()) # normalize to max
xticks([0, 0.5, 1])
yticks([0, count.max()])
xlabel('last spike time fractions')
ylabel('neuron count')
gcfm().set_window_title('last_spike_fractions')
gcf().tight_layout(pad=0.3) # crop figure to contents

show()
