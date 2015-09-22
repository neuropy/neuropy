"""Plot distribution of mean rates of all cells across all 6 natstate recordings. Run with
'run -i scripts/meanrates_natstate.py'"""

from __future__ import division
from __future__ import print_function

from scipy.stats import chisquare, mannwhitneyu
from scipy.optimize import leastsq

import pylab as pl

from core import ceilsigfig, g


figsize = (3, 3) # inches
NULLRATE = 1e-4
ymax = 40 # for meanrates hist
logmin, logmax = -4, 2
ticks = 10**(np.arange(logmin, logmax+1.0, 2.0)) # label every even power of 10
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 in log space
logbinwidth = log10(bins[1])-log10(bins[0])
midbins = np.logspace(logmin+logbinwidth/2, logmax-logbinwidth/2, nbins) # nbins in log space
logmidbins = log10(midbins)

## NOTE: this analysis includes all isolated single units, and is not limited only to
## responsive neurons

def cost(p, x, y):
    """Cost function for LM least-squares fit"""
    mu, sigma, A = p
    return A * g(mu, sigma, x) - y


# sort recordings by their absname:
urecs = [ eval(recname) for recname in sorted(REC2STATETRANGES) ] # unique, no reps, sorted

meanrates = [[], []] # desynched, then synched

for rec in urecs:
    print(rec.absname)
    stranges = REC2STATETRANGES[rec.absname]
    for statei, strange in enumerate(stranges): # desynched, then synched
        dt = (strange[1] - strange[0]) / 1e6 # s
        snids = sorted(rec.alln)
        for nid in snids:
            n = rec.alln[nid]
            si0, si1 = n.spikes.searchsorted(strange)
            nspikes = si1 - si0
            meanrate = nspikes / dt # Hz
            if meanrate == 0:
                meanrate = NULLRATE # replace 0s with NULLRATE
            meanrates[statei].append(meanrate)

meanrates = np.asarray(meanrates)
desynchrates = meanrates[0][meanrates[0] != NULLRATE] # filter out any NULLRATE values
synchrates = meanrates[1][meanrates[1] != NULLRATE]

# plot mean rate distributions in log space:
figure(figsize=figsize)
# plot actual distributions:
shist = hist(synchrates, bins=bins, histtype='step', color='r', zorder=10)[0] # synched
dhist = hist(desynchrates, bins=bins, histtype='step', color='b', zorder=10)[0] # desynched
# do some stats:
u, p = mannwhitneyu(log10(synchrates), log10(desynchrates)) # 1-sided
slogmean = log10(synchrates).mean()
dlogmean = log10(desynchrates).mean()
slogstd = log10(synchrates).std()
dlogstd = log10(desynchrates).std()
smean = 10**(slogmean) # geometric
dmean = 10**(dlogmean)
sstd = 10**(slogstd)
dstd = 10**(dlogstd)
print('smean, dmean:', smean, dmean)
print('slogstd, dlogstd:', slogstd, dlogstd)
'''
# plot approximately equivalent lognormal distribution:
A, mu = 34, slogmean
sigma = log10(synchrates).std()
modelx = np.logspace(logmin, logmax, 200) # 200 points in log space
modely = A*g(mu, sigma, log10(modelx))
plot(modelx, modely, 'e-', lw=1)
'''
# Do least-squares LM fit of a lognormal distribution to the data.
# Everything is done in log space:
# set initial parameters:
sA, dA = max(shist), max(dhist) # Gaussian amplitude
sp = slogmean, slogstd, sA # parameter list
dp = dlogmean, dlogstd, dA
# fit the 3 parameters for both synched and desynched distribs:
sresult = leastsq(cost, sp, args=(logmidbins, shist), full_output=True)
dresult = leastsq(cost, dp, args=(logmidbins, dhist), full_output=True)
sp, cov_p, infodict, mesg, ier = sresult
dp, cov_p, infodict, mesg, ier = dresult
smu, ssigma, sA = sp
dmu, dsigma, dA = dp
# generate fit distributions from fit parameters
modelx = np.logspace(logmin, logmax, 200) # 200 points in log space
smodel = sA*g(smu, ssigma, log10(modelx))
dmodel = dA*g(dmu, dsigma, log10(modelx))
# plot the fit distributions:
plot(modelx, smodel, 'r--', lw=1, alpha=0.5)
plot(modelx, dmodel, 'b--', lw=1, alpha=0.5)
xscale('log')
xlim(xmin=10**logmin, xmax=10**logmax)
ylim(ymax=ymax)
xticks(ticks)
yticks(range(0, 40+10, 10))
xlabel('mean firing rate (Hz)')
ylabel('unit count')
# display geometric means and p value:
text(0.03, 0.98, '$\mu$ = %.2f Hz' % smean, # synched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='r')
text(0.03, 0.90, '$\mu$ = %.2f Hz' % dmean, # desynched
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='b')
text(0.03, 0.82, 'p < %.1g' % ceilsigfig(p, 1),
                 horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
# arrow doesn't display correctly on log axis, use annotate instead:
annotate('', xy=(smean, (6/7)*40), xycoords='data', # synched
             xytext=(smean, 40), textcoords='data',
             arrowprops=dict(fc='r', ec='none', width=1.3, headwidth=7, frac=0.5))
annotate('', xy=(dmean, (6/7)*40), xycoords='data', # desynched
             xytext=(dmean, 40), textcoords='data',
             arrowprops=dict(fc='b', ec='none', width=1.3, headwidth=7, frac=0.5))
titlestr = 'meanrates hist'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

# scatter plot rates in desynched vs synched state
figure(figsize=figsize)
truecols = (meanrates != NULLRATE).all(axis=0) # columns without NULLRATE values
falsecols = (meanrates == NULLRATE).any(axis=0) # columns with NULLRATE values
meanratestrue = meanrates[:, truecols]
meanratesfalse = meanrates[:, falsecols]
# report numbers, fractions and chi2 p values for reliability scatter plot.
nbelowyxline = (meanratestrue[1, :] > meanratestrue[0, :]).sum()
naboveyxline = (meanratestrue[0, :] > meanratestrue[1, :]).sum()
fractionbelowyxline = nbelowyxline / (nbelowyxline + naboveyxline)
chi2, p = chisquare([naboveyxline, nbelowyxline])
pstring = 'p < %g' % ceilsigfig(p)
print('nbelowyxline=%d, naboveyxline=%d, fractionbelowyxline=%.3g, '
      'chi2=%.3g, p=%.3g' % (nbelowyxline, naboveyxline,
                             fractionbelowyxline, chi2, p))
plot([10**logmin, 10**logmax], [10**logmin, 10**logmax], 'e--') # plot y=x line
plot(meanratestrue[1], meanratestrue[0], 'o', mec='k', mfc='None')
plot(meanratesfalse[1], meanratesfalse[0], 'o', mec='e', mfc='None')
xlabel('synchronized rate (Hz)')
ylabel('desynchronized rate (Hz)')
xscale('log')
yscale('log')
xlim(10**(logmin-0.1), 10**logmax)
ylim(10**(logmin-0.1), 10**logmax)
xticks(ticks)
yticks(ticks)
text(0.03, 0.98, '%s' % pstring, horizontalalignment='left', verticalalignment='top',
                 transform=gca().transAxes, color='k')
titlestr = 'meanrates scatter'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

pl.show()
