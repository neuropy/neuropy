"""Plot runspeed"""

from scipy.io import loadmat

rsfname = '/home/mspacek/data/blab/natstate/Ntsr1-Cre_0174/tr2/05/Ntsr1-Cre_0174_s02_e05_runspeed.mat'
rec = nts174.tr2.r05

# load runspeed info:
rsd = loadmat(rsfname, squeeze_me=True) # dict
tspeed, speed = rsd['tspeed'], rsd['speed'] # sec, cm/s, speed can contain nans

# smooth runspeed using overlapping time bins:
width, tres = 10, 0.5 # s
minspeed = 1 # cm/s, otherwise considered at rest
rec.lfp.get_data() # make sure it's loaded so we can access t0 and t1
t0 = rec.lfp.t0 / 1e6 # sec, tspeed starts at t=0, lfp starts a few ms later
t1 = rec.lfp.t1 / 1e6
tranges = core.split_tranges([(t0, t1)], width, tres) # overlapping time bins
tiranges = tspeed.searchsorted(tranges)
nbins = len(tiranges)
tbinspeed = tranges[:, 0]
binspeed = np.zeros(nbins)
for bini, (t0i, t1i) in enumerate(tiranges):
    binspeed[bini] = np.nanmean(speed[t0i:t1i]) # handle nans with nanmean
show()

f = figure(figsize=(1, 3.5))
plot(binspeed, tbinspeed, 'k-')
ylim(t1, t0)
xlim(xlim()[1], xlim()[0])
f.tight_layout(pad=0.3) # crop figure to contents
titlestr = rec.absname + '_runspeed'
gcfm().window.setWindowTitle(titlestr)

# plot runspeed as a color map:
binspeed.shape = -1, 1 # column vector
if minspeed:
    # turn it into a binary rest/run signal:
    binspeed[binspeed < minspeed] = 0 # at rest
    binspeed[binspeed >= minspeed] = 1 # running
f = figure(figsize=(1, 10))
axis('off')
plt.imshow(binspeed, aspect=0.025, cmap='gray_r') # white=rest, black=run
#plt.imshow(binspeed, aspect=0.025, cmap='jet') # blue=rest, red=run
titlestr = rec.absname + '_runspeed_cmap'
gcfm().window.setWindowTitle(titlestr)
f.tight_layout(pad=0.3) # crop figure to contents
show()
