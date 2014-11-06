"""Calculate autocorr of each neuron, find number of spike intervals <= RP (refractory period
violations, RPVs). For each neuron pair, temporarily merge the spike trains of the two
neurons. Calculate autocorr of merger, and find its number of RPVs. Take difference between
this and the sum of the individual cell RPVs, normalize by total number of ISIs to get change
in fraction of refractory period violations df. Sort neuron pairs in decreasing order of df.
Those with the biggest df have an autocorr that's most affected by the incorrect merger. This
should help demonstrate that even in the worst case scenario, the autocorr isn't a great
diagnostic of contaminated spike trains, and in the more typical case, it's an absolutely
terrible diagnostic of contaminated spike trains. Note that in the numerator of the df
calcalation, I should really be counting number of consecutive ISIs <= RP, not number of spike
intervals in general in the autocorr <= RP, but the two are very nearly the same thing. If
there's one spike in that interval, it's very unlikely for there to be another, and so ISI and
autocorr for such a short interval are effectively the same thing.

Meant to be copied and pasted into the neuropy console.

TODO: taking xcorr really isn't necessary, only really need to take diff of each spike train
and count the number of consecutive intervals (ISIs) that fall below RP, but this is the way I
happened to do it at the time, and I don't want to regenerate my plots at the moment, so I'm
leaving it be for now.

TODO: instead of normalizing by number of ISIs (which is typically very large compared to
nRPVs), normalize by the expected coincidence rate of a historyless Poisson spike train with a
mean firing rate that corresponds to each cell.
"""

import util

RP = 750 # refractory period, us
REXCL = 0 # exclusion radius, um
tracks = ptc15.tr7c, ptc22.tr1, ptc22.tr2
colours = 'r', 'g', 'b'

farrs = [] # collect f arrays
fss = [] # collect fs arrays
# refractory period violation array dtype, one row per neuron pair,
# each row: [nid0, nid1, f0, f1, df]
farrdtype=[('nid0', np.int64), ('nid1', np.int64),
           ('f0', np.float64), ('f1', np.float64), ('df', np.float64)]
# convert RP to half trange array, in us
trange = np.array([0, RP]) # only need to work on one half of the autocorrelogram

for track in tracks:
    print(track.absname)
    neurons = track.alln # dict
    nids = sorted(neurons) # sorted keys
    nn = len(nids)
    # precalculate nrpvs and f values for reuse in pair loop:
    nrpvs = np.zeros(nn, dtype=np.int64) # number of refractory period violations, per nid
    fs = np.zeros(nn, dtype=np.float64) # f values, per nid
    for nidi, nid in enumerate(nids):
        s = neurons[nid].spikes # should be sorted
        dts = util.xcorr(s, s, trange)
        dts = dts[dts != 0] # remove 0s for autocorr
        nrpvs[nidi] = len(dts)
        nISI = len(s) - 1
        fs[nidi] = nrpvs[nidi] / nISI
    fss.append(fs)
        
    maxnpairs = nCr(nn, 2) # maximum number of pairs, actual will be less given REXCL
    farr = np.zeros(maxnpairs, dtype=farrdtype)
    pairi = 0
    for nidi0 in range(nn):
        for nidi1 in range(nidi0+1, nn): # for all neuron pairs
            nid0 = nids[nidi0]
            nid1 = nids[nidi1]
            n0, n1 = neurons[nid0], neurons[nid1]
            if core.dist(n0.pos, n1.pos) <= REXCL:
                continue # skip this pair, they're too close together
            s0 = n0.spikes # should be sorted
            s1 = n1.spikes # should be sorted
            # find dts for the merged spike train:
            s = np.concatenate([s0, s1])
            s.sort() # necessary for xcorr
            mdts = util.xcorr(s, s, trange) # dts of merged spike trains
            mdts = mdts[mdts != 0] # remove 0s for autocorr
            nrpv0, nrpv1, nrpvm = nrpvs[nidi0], nrpvs[nidi1], len(mdts)
            nISI = len(s) - 1
            f0, f1 = fs[nidi0], fs[nidi1]
            df = (nrpvm - (nrpv0 + nrpv1)) / nISI
            farr[pairi] = nid0, nid1, f0, f1, df
            #f.append(nid0, nid1, f0, f1, df)
            pairi += 1
            print('.', end='')
    npairs = pairi
    print()
    print('maxnpairs = %d' % maxnpairs)
    print('actual npairs = %d' % npairs)
    print('fs descending: %r' % np.sort(fs)[::-1])
    farr = farr[:npairs] # trim down to actual number of pairs
    sortis = np.argsort(farr['df'])
    sortis = sortis[::-1] # reverse
    farr = farr[sortis] # sorted by decreasing df values
    print(farr)
    farrs.append(farr)

fontsize(16)
fig = figure(figsize=(7, 4.4))
for track, farr, c in zip(tracks, farrs, colours):
    plot(farr['df'], ls='-', lw=4, marker=None, c=c, label=track.absname)
legend(frameon=False, loc='center right')
xlim(xmin=-20)
ylim(ymin=-0.0002)
xlabel('pair rank')
ylabel('df')
fig.tight_layout(pad=0.3) # crop figure to contents
