"""Generate a histogram of distance of each neuron to its nearest polytrode site, from all
specified tracks. Run from within neuropy using `run -i scripts/ndistpdf.py`"""

from __future__ import division
from core import dist

figsize = 2.75, 3 # inches
spacing = 65 # um
binw = 2 # um
edges = np.arange(0, spacing/2+binw, binw)
print(edges)

c = 'black'
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
# ptc15 is 2a, ptc22 tracks are both 1a
# calculate recordable areas for each polytrode as half a horizontal or vertical spaceing at
# each edge (top, bottom, left right):
areas = []
for track in tracks:
    xs = np.unique(track.chanpos[:, 0])
    ys = np.unique(track.chanpos[:, 1])
    dx = np.diff(xs)[0]
    dy = np.diff(ys)[0]
    w = max(xs) - min(xs) + dx
    h = max(ys) - min(ys) + dy
    areas.append(w*h) # in square um


N = 0
ds = []
for track in tracks:
    neurons = track.alln.values()
    for neuron in neurons:
        d = [ dist(neuron.pos, cp) for cp in neuron.sort.chanpos ]
        ds.append(min(d))
    N += len(neurons)
ds = np.hstack(ds)

rho = N / sum(areas) # average cell density per unit polytrode area (1/um^2)
theory = rho*2*pi*edges*N # what you'd expect for randomly distributed units

figure(figsize=figsize)
# underplot hlines at electrode site y positions:
n = hist(ds, bins=edges, color=c, ec=c)[0]
plot(edges, theory, 'e--', lw=3)
xlabel('nearest site distance ($\mu$m)')
ylabel('count')
#xticks((0, gca().get_xticks().max())) # get rid of intermediate xicks
tight_layout(0.3)
gcfm().set_window_title("ndistpdf")

show()
