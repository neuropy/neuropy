"""Plot distribution of spatial sigmas of designated tracks. Run from within neuropy
using `run -i scripts/sigma_distrib.py`"""

from __future__ import division
from __future__ import print_function

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
bw = 5
bins = np.arange(0, 100+bw, bw)

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
hist(sigmas, fc='k', bins=bins) # bins=40 looks a little more tantalizing
xlim(xmax=110)
ylim(ymax=50)
xticks((0, 25, 50, 75, 100))
yticks((0, 25, 50))
xlabel('$\sigma$ ($\mu$m)')
ylabel('neuron count')
gcfm().window.setWindowTitle('all tracks')
gcf().tight_layout(pad=0.2) # crop figure to contents
show()
