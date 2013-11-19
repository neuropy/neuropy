"""For each neuron pair, calculate autocorr of each neuron, find number of spike intervals <=
RP, sum this up for both neurons, then temporarily merge the spike trains of the two neurons.
Calculate autocorr of merger, and find its number of spike intervals <= RP. Take difference
between this and the previous sum, normalize by total number of ISIs to get change in fraction
of refractory period violations df. Sort neuron pairs in decreasing order of df. Those with
the biggest df have an autocorr that's most affected by the incorrect merger. This should help
demonstrate that even in the worst case scenario, the autocorr is a poor diagnostic of
contaminated spike trains. Not that in the numerator of the df calcalation, I should really be
counting number of consecutive ISIs <= RP, not number of spike intervals in general in the
autocorr <= RP, but the two are very nearly the same thing. If there's one spike in that
interval, it's very unlikely for there to be another, and so ISI and autocorr for such a short
interval are effectively the same thing.

Meant to be copied and pasted into the neuropy console."""

import util

RP = 2000 # refractory period, us
REXCL = 0 # exclusion radius, um
tracks = ptc15.tr7c, ptc22.tr2, ptc22.tr1 # use this order so ptc22.tr1 comes out on top
colours = 'lightgrey', 'grey', 'black'

fontsize(16)
fig = figure(figsize=(7, 4.4))

fs = [] # collect f arrays
# refractory period violation array dtype, one row per neuron pair,
# each row: [nid0, nid1, f0, f1, df]
fdtype=[('nid0', np.int64), ('nid1', np.int64),
        ('f0', np.float64), ('f1', np.float64), ('df', np.float64)]

for tracki, (track, c) in enumerate(zip(tracks, colours)):
    print(track.absname)
    neurons = track.alln # dict
    nids = sorted(neurons) # sorted keys
    nn = len(nids)
    maxnpairs = nCr(nn, 2) # maximum number of pairs, actual will be less given REXCL
    f = np.zeros(maxnpairs, dtype=fdtype)
    # convert RP to half trange array, in us
    trange = np.array([0, RP]) # only need to work on one half of the autocorrelogram
    pairi = 0
    for nidi0 in range(nn):
        for nidi1 in range(nidi0+1, nn): # for all neuron pairs
            nid0 = nids[nidi0]
            nid1 = nids[nidi1]
            n0, n1 = neurons[nid0], neurons[nid1]
            if core.dist(n0.pos, n1.pos) <= REXCL:
                continue # skip this pair, they're too close together
            s0 = n0.spikes # should already be sorted
            s1 = n1.spikes # should already be sorted
            # taking xcorr here really isn't necessary, only really need to take
            # diff of each spike train and count the number of consecutive intervals
            # (ISIs) that fall below RP, but this is the way I happened to do it at the
            # time, and I don't want to regenerate my plots at the moment, so I'm leaving
            # it be for now:
            ## TODO: this is really inefficient. I should calculate the autocorr of each
            ## cell only once, and then reuse it over and over. Same goes for f0 and f1
            dts0 = util.xcorr(s0, s0, trange)
            dts0 = dts0[dts0 != 0] # remove 0s for autocorr
            dts1 = util.xcorr(s1, s1, trange)
            dts1 = dts1[dts1 != 0] # remove 0s for autocorr
            # now do it for the merged spike trains:
            s = np.concatenate([s0, s1])
            s.sort() # necessary for xcorr
            dts = util.xcorr(s, s, trange)
            dts = dts[dts != 0] # remove 0s for autocorr
            ndts0, ndts1, ndts = len(dts0), len(dts1), len(dts)
            nISI0 = len(s0) - 1
            nISI1 = len(s1) - 1
            nISI = len(s) - 1
            f0 = ndts0 / nISI0
            f1 = ndts1 / nISI1
            df = (ndts - (ndts0 + ndts1)) / nISI
            f[pairi] = nid0, nid1, f0, f1, df
            #f.append(nid0, nid1, f0, f1, df)
            pairi += 1
            print('.', end='')
    npairs = pairi
    print()
    print('maxnpairs = %d' % maxnpairs)
    print('actual npairs = %d' % npairs)
    f = f[:npairs] # trim down to actual number of pairs
    sortis = np.argsort(f['df'])
    sortis = sortis[::-1] # reverse
    f = f[sortis] # sorted by decreasing df values
    print(f)
    plot(f['df'], ls='-', lw=4, marker=None, c=c, label='track %d' % (tracki+1))
    fs.append(f)

legend(frameon=False, loc='center right')
xlim(xmin=-20)
ylim(ymin=-0.0002)
xlabel('pair rank')
ylabel('df')
fig.tight_layout(pad=0.3) # crop figure to contents
