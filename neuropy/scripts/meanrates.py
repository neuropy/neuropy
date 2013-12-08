"""Plot mean firing rates of all cells in each track, in descending order, on a log scale.
Show firing rate threshold as well. Also, plot distribution of mean rates across all cells,
again on a log scale. Copy and paste into neuropy console"""

fontsize(18)
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2]
colours = 'r', 'g', 'b'
nbins = 25

# Figure 1: plot firing rate rank order:
f1 = figure(1)
allrates = []
for track, c in zip(tracks, colours):
    trackrates = sorted([ n.meanrate for n in track.alln.values() ])[::-1] # reverse order
    plot(trackrates, '--', c=c, linewidth=3, label=track.absname)
    trackrates = np.asarray(trackrates)
    n = len(trackrates)
    print(track.absname)
    print('%.2f %% neurons < 1 Hz' % ((trackrates < 1).sum() / n * 100))
    print('%.2f %% neurons < %.2f Hz' % ((trackrates < MINRATE).sum() / n * 100, MINRATE))
    allrates.append(trackrates)
allrates = np.hstack(allrates)
allrates.sort()
allrates = allrates[::-1] # decreasing order
n = len(allrates)

# plot rank order of all cells across all tracks, in black:
allx = np.arange(n) / len(tracks) # compress x values to fit, effectively the mean
plot(allx, allrates, '--', c='k', linewidth=3, label='mean')

print('%.2f %% total neurons < 1 Hz' % ((allrates < 1).sum() / n * 100))
print('%.2f %% total neurons < %.2f Hz' % ((allrates < MINRATE).sum() / n * 100, MINRATE))

# Figure 2: plot meanrates PDF, distribution as y for modelling
f2 = figure(2)
logallrates = np.log10(allrates)
logmin, logmax = min(logallrates), max(logallrates)
logstart, logend = np.floor(logmin), np.ceil(logmax)
edges = np.logspace(logstart, logend, nbins+1) # nbins+1 points in log space

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
modelx = np.logspace(logstart, logend, 200) # 200 points in log space
modely = A*g(mu, sigma, np.log10(modelx))

# Figure 2: plot lognormal model in magenta:
plot(modelx, modely, 'm-', linewidth=3)

# Figure 2: finish up
axvline(x=MINRATE, c='e', ls='--', marker=None) # plot threshold
xscale('log')
xlabel('mean firing rate (Hz)')
ylabel('neuron count')
text(0.98, 0.98, 'log($\mu$) = %.2f\n'
                 'log($\sigma$) = %.2f' % (mu, sigma),
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=gca().transAxes,
                 color='m')
gcfm().window.setWindowTitle('meanrates_pdf')
f2.tight_layout(pad=0.3) # crop figure to contents

# Figure 1: now plot rank order of lognormal distrib of cells across all tracks, in magenta:
figure(1)
modelx, modely
rank = np.cumsum(modely)
rank = rank / rank.max() * allx.max() # make it range
rank = rank[::-1] # decreasing order
plot(rank, modelx, 'm-', linewidth=3, label='lognormal')

# Figure 1: finish up
axhline(y=MINRATE, c='e', ls='--', marker=None) # plot threshold
legend(frameon=False)
yscale('log')
xlabel('neuron rank')
ylabel('mean firing rate (Hz)')
gcfm().window.setWindowTitle('meanrates_rank')
f1.tight_layout(pad=0.3) # crop figure to contents
