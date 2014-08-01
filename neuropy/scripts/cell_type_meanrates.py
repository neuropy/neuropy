"""Plot distributions of mean rates across all cells, according to spike and RF type, 
on a log scale. Run with 'run -i scripts/cell_type_meanrates.py' or copy and paste into neuropy
console"""

from __future__ import division
from __future__ import print_function


tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
nbins = 15

spiketyperates = {'fast':[], 'slow':[], 'fastasym':[], 'slowasym':[]}
rftyperates = {'simple':[], 'complex':[], 'LGN':[], None:[]}
allrates = []
for track in tracks:
    for n in track.alln.values():
        spiketyperates[n.spiketype].append(n.meanrate)
        rftyperates[n.rftype].append(n.meanrate)
        allrates.append(n.meanrate)

logallrates = np.log10(allrates)
logmin, logmax = min(logallrates), max(logallrates)
logstart, logend = np.floor(logmin), np.ceil(logmax)
edges = np.logspace(logstart, logend, nbins+1) # nbins+1 points in log space

figure(figsize=(3, 3))
hist(spiketyperates['fast'], bins=edges, color='r')
hist(spiketyperates['slow'], bins=edges, color='b')
hist(spiketyperates['fastasym'], bins=edges, color='g')
hist(spiketyperates['slowasym'], bins=edges, color='e')
ymin, ymax, lw = 28, 30, 2
vlines(10**mean(log10(spiketyperates['fast'])), ymin, ymax, colors='r', lw=lw)
vlines(10**mean(log10(spiketyperates['slow'])), ymin, ymax, colors='b',  lw=lw)
vlines(10**mean(log10(spiketyperates['fastasym'])), ymin, ymax, colors='g', lw=lw)
vlines(10**mean(log10(spiketyperates['slowasym'])), ymin, ymax, colors='e', lw=lw)
xscale('log')
#ylim(0, y.max())
#yticks((0, y.max()))
xlabel('mean firing rate (Hz)')
ylabel('neuron count')
gcfm().window.setWindowTitle('spike_type_meanrates')
tight_layout(pad=0.3) # crop figure to contents

figure(figsize=(3, 3))
hist(rftyperates['simple'], bins=edges, color='r')
hist(rftyperates['complex'], bins=edges, color='b')
hist(rftyperates[None], bins=edges, color='e')
hist(rftyperates['LGN'], bins=edges, color='g')
ymin, ymax, lw = 23, 25, 2
vlines(10**mean(log10(rftyperates['simple'])), ymin, ymax, colors='r', lw=lw)
vlines(10**mean(log10(rftyperates['complex'])), ymin, ymax, colors='b',  lw=lw)
vlines(10**mean(log10(rftyperates['LGN'])), ymin, ymax, colors='g', lw=lw)
vlines(10**mean(log10(rftyperates[None])), ymin, ymax, colors='e', lw=lw)
xscale('log')
#ylim(0, y.max())
#yticks((0, y.max()))
xlabel('mean firing rate (Hz)')
ylabel('neuron count')
gcfm().window.setWindowTitle('RF_type_meanrates')
tight_layout(pad=0.3) # crop figure to contents

show()
