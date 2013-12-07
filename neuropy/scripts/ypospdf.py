"""Generate a histogram of y position of 3 tracks, such that it's to scale with plots
from track.npos"""

figsize = 1.4, 10.57 # inches

c = 'black'
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]

# for 1a 3 column and 2a 2 column probes, with 65 um site spacing (32.5 vertical spacing),
# and aligned sites:
spacing = 65 / 2 # um
binw = spacing / 2 # um
binw2 = binw / 2
ymin, ymax = 32, 1755
ledge = ymin-binw-binw2 # this will become the right edge once axes are inverted
redge = ymax+binw-binw2 # this will become the left edge once axes are inverted
edges = np.arange(ledge, redge+2*binw, binw)
#edges = np.arange(-56-14, 56+2*14+14, 14) - 14/2

ypos = []
chanypos = []
for track in tracks:
    neurons = track.alln.values()
    ypos.append([ n.pos[1] for n in neurons ]) # all y position values
    chanypos.append(track.sort.chanpos[:, 1])
ypos = np.hstack(ypos)
chanypos = np.unique(np.hstack(chanypos))

figure(figsize=figsize)
# underplot hlines at electrode site y positions:
for y in chanypos:
    axhline(y=y, c='e', ls='--', marker=None, zorder=-1)
n = hist(ypos, bins=edges, color=c, ec=c, orientation='horizontal')[0]
gca().invert_yaxis()
xticks((0, gca().get_xticks().max())) # get rid of intermediate xicks
tight_layout(0.3)
gcfm().set_window_title("ypospdf")

print('%d cells are site aligned' % n[1::2].sum())
print('%d cells are inter-site aligned' % n[0::2].sum())
