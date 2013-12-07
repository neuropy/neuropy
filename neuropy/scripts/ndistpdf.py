"""Generate a histogram of distance of each neuron to its nearest polytrode site, from all
specified tracks"""

from core import dist

figsize = 2.75, 3 # inches
spacing = 65 # um
binw = 2 # um
edges = np.arange(0, spacing/2+binw, binw)

c = 'black'
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]

ds = []
for track in tracks:
    neurons = track.alln.values()
    for neuron in neurons:
        d = [ dist(neuron.pos, cp) for cp in neuron.sort.chanpos ]
        ds.append(min(d))
ds = np.hstack(ds)

figure(figsize=figsize)
# underplot hlines at electrode site y positions:
n = hist(ds, bins=edges, color=c, ec=c)[0]
xlabel('nearest site distance ($\mu$m)')
ylabel('count')
#xticks((0, gca().get_xticks().max())) # get rid of intermediate xicks
tight_layout(0.3)
gcfm().set_window_title("ndistpdf")
