"""Calculate average ISI histogram for all possible merged pairs of spike trains. Probably not
all that useful. The peaks in these plots are probably dominated by bursty high firing rate
cells"""

nids = np.sort(rec.n.keys())
nn = len(nids)
binw = 1000 # us
bins = np.arange(0, 100000+binw, binw)
hists = []
for nii0 in range(nn):
    for nii1 in range(nii0+1, nn):
        spikes0 = rec.n[nids[nii0]].spikes
        spikes1 = rec.n[nids[nii1]].spikes
        spikes = np.hstack([spikes0, spikes1]) # merge
        spikes.sort()
        isi = np.diff(spikes) # in us
        hist = np.histogram(isi, bins=bins)[0]
        hist = hist / hist.sum() # pmf: normalize so that sum of each hist is 1
        hists.append(hist)
hists = np.vstack(hists).sum(axis=0)
#hists = np.vstack(hists).mean(axis=0)

figure()
bar(bins[:-1]/1000, hists, width=binw/1000)
xlabel('ISI (ms)')
ylabel('mean count across all pairs')
tight_layout(pad=0.3)
