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
spiketypeadfs = {'fast':[], 'slow':[], 'fastasym':[], 'slowasym':[]} # adfs by spike type
rftypeadfs = {'simple':[], 'complex':[], 'LGN':[], None:[]} # adfs by RF type
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
        adf = n.dt / totaltrackdt
        spiketypeadfs[n.spiketype].append(adf)
        rftypeadfs[n.rftype].append(adf)
        adfs.append(adf)
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

# plot active duration distributions according to cell type:
histtype = 'step'
lw = 2

figure(figsize=(3, 3))
hist(spiketypeadfs['fast'], bins=edges, histtype=histtype, lw=lw, color='r')
hist(spiketypeadfs['slow'], bins=edges, histtype=histtype, lw=lw, color='b')
hist(spiketypeadfs['fastasym'], bins=edges, histtype=histtype, lw=lw, color='g')
hist(spiketypeadfs['slowasym'], bins=edges, histtype=histtype, lw=lw, color='e')
y, dy, hw, hl = 40, -10, 0.05, 5
arrow(mean(spiketypeadfs['fast']), y, 0, dy, color='r',
      head_width=hw, head_length=hl, length_includes_head=True)
arrow(mean(spiketypeadfs['slow']), y+dy/2, 0, dy, color='b',
      head_width=hw, head_length=hl, length_includes_head=True)
arrow(mean(spiketypeadfs['fastasym']), y, 0, dy, color='g',
      head_width=hw, head_length=hl, length_includes_head=True)
arrow(mean(spiketypeadfs['slowasym']), y, 0, dy, color='e',
      head_width=hw, head_length=hl, length_includes_head=True)
xlim(0, 1)
#ylim(ymax=count.max()) # normalize to max
xticks([0, 0.5, 1])
#yticks([0, count.max()])
xlabel('fractional active duration')
ylabel('neuron count')
gcfm().set_window_title('spiketype_active_duration')
gcf().tight_layout(pad=0.3) # crop figure to contents

figure(figsize=(3, 3))
hist(rftypeadfs['simple'], bins=edges, histtype=histtype, lw=lw, color='r')
hist(rftypeadfs['complex'], bins=edges, histtype=histtype, lw=lw, color='b')
hist(rftypeadfs[None], bins=edges, histtype=histtype, lw=lw, color='e')
hist(rftypeadfs['LGN'], bins=edges, histtype=histtype, lw=lw, color='g')
y, dy, hw, hl = 38, -9, 0.05, 4
arrow(mean(rftypeadfs['simple']), y, 0, dy, color='r',
      head_width=hw, head_length=hl, length_includes_head=True)
arrow(mean(rftypeadfs['complex']), y+dy/2, 0, dy, color='b',
      head_width=hw, head_length=hl, length_includes_head=True)
arrow(mean(rftypeadfs['LGN']), y, 0, dy, color='g',
      head_width=hw, head_length=hl, length_includes_head=True)
arrow(mean(rftypeadfs[None]), y, 0, dy, color='e',
      head_width=hw, head_length=hl, length_includes_head=True)
xlim(0, 1)
#ylim(ymax=count.max()) # normalize to max
xticks([0, 0.5, 1])
#yticks([0, count.max()])
xlabel('fractional active duration')
ylabel('neuron count')
gcfm().set_window_title('rftype_active_duration')
gcf().tight_layout(pad=0.3) # crop figure to contents

show()
