"""Plot mean firing rates of all cells in each track, in descending order, on a log scale.
Show firing rate threshold as well. Copy and paste into neuropy console"""

fontsize(16)
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
colours = 'r', 'g', 'b'
f = figure()
for track, c in zip(tracks, colours):
    trackrates = sorted([ n.meanrate for n in track.alln.values() ])[::-1] # reverse order
    plot(trackrates, '-', c=c, linewidth=3, label=track.absname)

axhline(y=MINRATE, c='e', ls='--', marker=None)
legend(frameon=False)
yscale('log')
xlabel('neuron count')
ylabel('mean firing rate (Hz)')
f.tight_layout(pad=0.3) # crop figure to contents
