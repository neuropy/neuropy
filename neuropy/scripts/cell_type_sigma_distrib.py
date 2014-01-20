"""Plot distribution of spatial sigmas of designated tracks. Run from within neuropy
using `run -i scripts/sigma_distrib.py`"""

from __future__ import division
from __future__ import print_function

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
bw = 5
bins = np.arange(0, 100+bw, bw)
alltracksbw = 4
alltracksbins = np.arange(0, 100+alltracksbw, alltracksbw)

sigmas = []
for track in tracks:
    tracksigmas = [ n.sigma for n in track.alln.values() ]
    sigmas.append(tracksigmas)
    figure(figsize=(3, 3))
    hist(tracksigmas, fc='k', bins=bins)
    xlim(xmax=110)
    ylim(ymax=20)
    yticks((0, 10, 20))
    xticks((0, 25, 50, 75, 100))
    xlabel('$\sigma$ ($\mu$m)')
    ylabel('neuron count')
    gcfm().window.setWindowTitle(track.absname)
    gcf().tight_layout(pad=0.2) # crop figure to contents

sigmas = np.hstack(sigmas)
figure(figsize=(3, 3))
hist(sigmas, fc='k', bins=alltracksbins)
xlim(xmax=110)
#ylim(ymax=40)
xticks((0, 25, 50, 75, 100))
yticks((0, 20, 40))
xlabel('$\sigma$ ($\mu$m)')
ylabel('neuron count')
gcfm().window.setWindowTitle('all tracks')
gcf().tight_layout(pad=0.2) # crop figure to contents
show()
