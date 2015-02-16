"""Generate a histogram of distance of each neuron to its nearest polytrode site, from all
specified tracks. Run from within neuropy using `run -i scripts/ndistpdf.py`"""

from __future__ import division
from core import dist

MODEL = False

distfigsize = 3.055, 3 # inches
scatterfigsize = 3, 3 # inches
spacing = 65 # um
binw = 2 # um
edges = np.arange(0, spacing/2+binw, binw)

c = 'black'
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
# ptc15 is 2a, ptc22 tracks are both 1a
# calculate recordable areas for each polytrode as half a horizontal or vertical spacing at
# each edge (top, bottom, left right):
xlims, ylims, areas = [], [], []
for track in tracks:
    xs = np.unique(track.chanpos[:, 0])
    ys = np.unique(track.chanpos[:, 1])
    dx = np.diff(xs)[0]
    dy = np.diff(ys)[0]
    w = max(xs) - min(xs) + dx
    h = max(ys) - min(ys) + dy
    xlims.append([min(xs)-dx/2, max(xs)+dx/2])
    ylims.append([min(ys)-dy/2, max(ys)+dy/2])
    areas.append(w*h) # in square um

N = 0
ds, sigmas = [], []
for track in tracks:
    neurons = track.alln.values()
    for neuron in neurons:
        d = [ dist(neuron.pos, cp) for cp in neuron.sort.chanpos ]
        ds.append(min(d))
        sigmas.append(neuron.sigma)
    N += len(neurons)
ds = np.hstack(ds)
sigmas = np.hstack(sigmas)

# calculate theoretical distribution, derived from concentric annuli:
rho = N / sum(areas) # average cell density per unit polytrode area (1/um^2)
midbins = edges[:-1] + intround(binw/2) # middle of each bin
theory = rho*2*pi*midbins/binw # normed distrib expected for randomly distributed units

# simulate numerically by randomly and uniformly distributing points around a polytrode and
# measuring the distribution of distances between each point and its nearest electrode site:
if MODEL: # can take a long time to run
    umap1a = ptc22.tr1.sort.chanpos
    area1a, xlims1a, ylims1a = areas[1], xlims[1], ylims[1]
    #nseed = intround(rho * area1a)
    nrandom = 100000
    w = xlims1a[1]-xlims1a[0]
    h = ylims1a[1]-ylims1a[0]
    x = np.random.random_sample(nrandom) * w - w/2 # centered around 0
    y = np.random.random_sample(nrandom) * h # all positive
    points = np.vstack([x, y]).T # each row is a point
    modelds = []
    for point in points:
        d = [ dist(point, cp) for cp in umap1a ]
        modelds.append(min(d))
    modelds = np.hstack(modelds)

figure(figsize=distfigsize)
vals = hist(ds, bins=edges, normed=True, color=c, ec=c)[0]
if MODEL:
    model = np.histogram(modelds, bins=edges, density=True)[0]
    plot(midbins, model, 'r--', lw=3)
plot(midbins, theory, 'b--', lw=3)
xlabel('nearest site distance ($\mu$m)')
ylabel('probability density (1/$\mu$m)')
#xticks((0, gca().get_xticks().max())) # get rid of intermediate xicks
xticks(np.arange(0, 30+10, 10))
ylim(ymax=vals.max())
tight_layout(0.3)
gcfm().set_window_title("ndistpdf")

figure(figsize=scatterfigsize)
plot(ds, sigmas, 'k.')
xlabel('nearest site distance ($\mu$m)')
ylabel('spatial extent ($\mu$m)')
xticks(np.arange(0, 30+10, 10))
yticks(np.arange(25, 100+25, 25))
tight_layout(0.3)
gcfm().set_window_title("spatial extent vs nearest site distance")

show()
