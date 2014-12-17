"""Comparison of STAs and trasters from drifting bars to check if STA polarity is correct
in the code"""

ptc15.tr7c.r70.sta().plot(normed=True, MPL=True)
# ori starts from 0, in steps of 20
ptc15.tr7c.r71.traster(nids=[213], c='k', hlinesweepis=range(1,18), hlinec='k')
ptc15.tr7c.r71.traster(nids=[228], c='k', hlinesweepis=range(1,18), hlinec='k')
ptc15.tr7c.r71.traster(nids=[328], c='k', hlinesweepis=range(1,18), hlinec='k')


ptc22.tr1.r04.sta().plot(normed=True, MPL=True)
# oris are np.arange(198, 360+198, 30) % 360
# [198, 228, 258, 288, 318, 348,  18,  48,  78, 108, 138, 168]
ptc22.tr1.r03.traster(nids=[24], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr1.r03.traster(nids=[131], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr1.r03.traster(nids=[1], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr1.r03.traster(nids=[17], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr1.r03.traster(nids=[22], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr1.r03.traster(nids=[40], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')

ptc22.tr1.r17.sta().plot(normed=True, MPL=True)
ptc22.tr1.r18.traster(nids=[20], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')


ptc22.tr2.r26.sta().plot(normed=True, MPL=True)
# oris are np.arange(270, 360+270, 30) % 360
# [270, 300, 330,   0,  30,  60,  90, 120, 150, 180, 210, 240]
ptc22.tr2.r25.traster(nids=[56], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr2.r25.traster(nids=[246], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
ptc22.tr2.r25.traster(nids=[168], c='bwg', hlinesweepis=range(2,24,2), hlinec='k')
