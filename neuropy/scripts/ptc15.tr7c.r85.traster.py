#for k, v in ptc15.tr7c.r85.e0.sweeptranges.items():
    #d = np.diff(v[0])
    #print(k, d)

# sfreq 0.8 cyc/deg, tfreq 2Hz, sweeptime 3s, contrast 1
#sweepis = np.arange(17, 185+24, 24)

# sfreq 0.15-5 cyc/deg, tfreq 2Hz, sweeptime 3s, contrast 0.5
# sweepis are 0-9, 24-33, 48-57, etc.
startis = np.arange(0, 24*8, 24)
l = 10 # len of each range
sweepis = np.hstack([ np.arange(starti, starti+l, 1) for starti in startis ])

ptc15.tr7c.r85.traster(nids='all', sweepis=sweepis, c='k')
