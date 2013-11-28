"""Generate little histograms of x position of 3 tracks, such that they're to scale with plots
from track.npos"""

figsize = 1.4, 1.4

# for 1a 3 column and 2a 2 column probes
binw = 14 # um
binw2 = binw / 2
ledge = -56-binw-binw2
redge = 56+binw-binw2
edges = np.arange(ledge, redge+2*binw, binw)
#edges = np.arange(-56-14, 56+2*14+14, 14) - 14/2

'''
# for 2a 2 column probes:
binw = 7
binw2 = binw / 2
ledge = -28-binw-binw2
redge = 28+binw-binw2
edges = np.arange(ledge, redge+2*binw, binw)
'''
# 2a 2 column probe:
ptc15.tr7c.pospdf(dim='x', edges=edges, figsize=figsize, labels=False)
xticks((-28, 28))
yticks((0, 45))
xlim(-80, 80)
ylim(0, 45)
gcfm().set_window_title("ptc15.tr7c.pospdf(dim='x')")

# 1a 3 column probe:
ptc22.tr1.pospdf(dim='x', edges=edges, figsize=figsize, labels=False)
xticks((-56, 0, 56))
yticks((0, 45))
xlim(-80, 80)
ylim(0, 45)
gcfm().set_window_title("ptc22.tr1.pospdf(dim='x')")

# 1a 3 column probe:
ptc22.tr2.pospdf(dim='x', edges=edges, figsize=figsize, labels=False)
xticks((-56, 0, 56))
yticks((0, 45))
xlim(-80, 80)
ylim(0, 45)
gcfm().set_window_title("ptc22.tr2.pospdf(dim='x')")
