"""Plot mean firing rates of all cells in each track, in descending order, on a log scale.
Show firing rate threshold as well. Also, plot distribution of mean rates across all cells,
again on a log scale. Copy and paste into neuropy console"""

fontsize(18)
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
colours = 'r', 'g', 'b'
nbins = 25
f = figure()
allrates = []
for track, c in zip(tracks, colours):
    trackrates = sorted([ n.meanrate for n in track.alln.values() ])[::-1] # reverse order
    plot(trackrates, '-', c=c, linewidth=3, label=track.absname)
    trackrates = np.asarray(trackrates)
    n = len(trackrates)
    print('%.2f %% neurons < 1 Hz' % ((trackrates < 1).sum() / n * 100))
    print('%.2f %% neurons < %.2f Hz' % ((trackrates < MINRATE).sum() / n * 100, MINRATE))
    allrates.append(trackrates)
allrates = np.hstack(allrates)

axhline(y=MINRATE, c='e', ls='--', marker=None)
legend(frameon=False)
yscale('log')
xlabel('neuron rank')
ylabel('mean firing rate (Hz)')
gcfm().window.setWindowTitle('meanrates_rank')
f.tight_layout(pad=0.3) # crop figure to contents

f = figure()
logallrates = np.log10(allrates)
logmin, logmax = min(logallrates), max(logallrates)
logstart, logend = np.floor(logmin), np.ceil(logmax)
edges = np.logspace(logstart, logend, nbins+1)
# plot meanrates PDF:
y = hist(allrates, bins=edges, color='k')[0]
logx = np.log10(edges)[:-1] # left bin edges
logx = logx + (logx[1] - logx[0]) / 2# middle of bins
# do least-squares LM of fit lognormal distribution to data:
from core import g
from scipy.optimize import leastsq
def cost(p, x, y):
    mu, sigma, A = p
    return A * g(mu, sigma, x) - y
logmean = -1
logsigma = 1
A = max(y)
p = logmean, logsigma, A
result = leastsq(cost, p, args=(logx, y), full_output=True)
p, cov_p, infodict, mesg, ier = result
mu, sigma, A = p
modelx = np.logspace(logstart, logend, 200)
# plot lognormal model:
plot(modelx, A*g(mu, sigma, np.log10(modelx)), 'r-', linewidth=3)
axvline(x=MINRATE, c='e', ls='--', marker=None)
xscale('log')
xlabel('mean firing rate (Hz)')
ylabel('neuron count')
text(0.98, 0.98, 'log($\mu$) = %.2f\n'
                 'log($\sigma$) = %.2f' % (mu, sigma),
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=gca().transAxes,
                 color='red')

gcfm().window.setWindowTitle('meanrates_pdf')
f.tight_layout(pad=0.3) # crop figure to contents
