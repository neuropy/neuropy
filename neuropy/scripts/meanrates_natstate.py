"""Plot distribution of mean rates of all (responsive?) cells across all 6 natstate
recordings. Run with 'run -i scripts/meanrates_natstate.py'"""

from __future__ import division
from __future__ import print_function

import pylab as pl

from core import ceilsigfig

figsize = (3, 3) # inches
NULLRATE = None

## TODO: limit to responsive neurons?
## TODO: limit to only when movie is on screen?
## TODO: scatter plot rates for each neuron?

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
                if NULLRATE:
                    meanrate = NULLRATE # replace 0s with NULLRATE for log scale plotting
                else:
                    continue # skip 0s
            meanrates[statei].append(meanrate)

for statei in range(2):
    meanrates[statei] = np.asarray(meanrates[statei])


# plot mean rate distributions in log space:
logmin, logmax = -4, 2
nbins = 20
bins = np.logspace(logmin, logmax, nbins+1) # nbins+1 points in log space
figure(figsize=figsize)
n1 = hist(meanrates[1], bins=bins, histtype='step', color='r')[0] # synched
n0 = hist(meanrates[0], bins=bins, histtype='step', color='b')[0] # desynched
n = np.hstack([n0, n1])
xlim(xmin=10**logmin, xmax=10**logmax)
ylim(ymax=n.max()+10)
#xticks(ticks)
xscale('log')
xlabel('mean firing rate (Hz)')
ylabel('unit count')
#t, p = ttest_ind(log10(meanrates[0]), log10(meanrates[1]), equal_var=False) # Welch's T-test
u, p = mannwhitneyu(log10(meanrates[0]), log10(meanrates[1])) # 1-sided
smean = 10**(log10(meanrates[1]).mean()) # geometric
dmean = 10**(log10(meanrates[0]).mean())
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
titlestr = 'meanrates'
gcfm().window.setWindowTitle(titlestr)
tight_layout(pad=0.3)

pl.show()
