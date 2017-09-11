"""Plot zoomed in regions of desynched and synched LFP, using same voltage and time scale"""

figsize = 5.5, 2


# cat:
ymin, ymax = -0.75, 0.75 # mV

# desynched:
ptc22.tr1.r08.lfp.plot(t0=500, t1=505, chanis=[-1], gain=1, yunits=None, title=False,
                       relative2t0=False, lim2stim=True, scalebar=False, figsize=figsize)
xlim(xmin=500)
ylim(ymin, ymax)
yticks([-0.5, 0, 0.5])
gcf().tight_layout(pad=0.3) # crop figure to contents
gcfm().window.setWindowTitle('ptc22.tr1.r08_desynched_LFP')
show()
# synched:
ptc22.tr1.r08.lfp.plot(t0=2000, t1=2005, chanis=[-1], gain=1, yunits=None, title=False,
                       relative2t0=False, lim2stim=True, scalebar=True, figsize=figsize)
xlim(xmin=2000)
ylim(ymin, ymax)
yticks([-0.5, 0, 0.5])
tight_layout(pad=0.3) # crop figure to contents
gcfm().window.setWindowTitle('ptc22.tr1.r08_synched_LFP')
show()

# mouse:
ymin, ymax = -0.5, 0.5 # mV

# synched:
nts174.tr2.r05.lfp.plot(t0=500, t1=505, chanis=[-1], gain=1, yunits=None, title=False,
                        relative2t0=False, lim2stim=True, scalebar=False, figsize=figsize)
xlim(xmin=500)
ylim(ymin, ymax)
yticks([-0.3, 0, 0.3])
gcf().tight_layout(pad=0.3) # crop figure to contents
gcfm().window.setWindowTitle('nts174.tr2.r05_synched_LFP')
show()

# desynched:
nts174.tr2.r05.lfp.plot(t0=2000, t1=2005, chanis=[-1], gain=1, yunits=None, title=False,
                        relative2t0=False, lim2stim=True, scalebar=True, figsize=figsize)
xlim(xmin=2000)
ylim(ymin, ymax)
yticks([-0.3, 0, 0.3])
tight_layout(pad=0.3) # crop figure to contents
gcfm().window.setWindowTitle('nts174.tr2.r05_desynched_LFP')
show()
