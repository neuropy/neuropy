"""Plot 2D matrix of cell types vs NS response. Run from within neuropy using `run -i
scripts/cell_type_vs_ns_matrix.py`"""

from scipy.stats import chisquare

alpha = 0.01 # for chi squared test

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
tracknames = [track.absname for track in tracks ] + ['all']

spiketype2int = {'fast':0, 'slow':1, 'fastasym':2, 'slowasym':3}
rftype2int = {'simple':4, 'complex':5, 'LGN':6, None: 7}


"""
Lists of cells that respond reliably. This is kind of a hack, should really make an .nsresp
file for every track, analogous to .spiketype and .rftype files:
"""
ns = {}

ns['ptc15.tr7c'] = [4, 33, 50, 63, 64, 85, 87, 91, 93, 116, 138, 145, 148, 161, 180, 185, 190,
206, 213, 221, 222, 228, 278, 328, 332, 340, 345, 414, 415, 427, 434, 479, 510, 514, 522, 540,
546, 558, 563, 569, 586, 642, 662, 674, 677, 691, 699]

ns['ptc22.tr1'] = [1, 2, 3, 5, 10, 17, 18, 20, 23, 24, 26, 27, 29, 32, 40, 43, 46, 50, 52, 53,
67, 74, 90, 92, 94, 96, 99, 100, 101, 102, 109, 111, 127, 128, 133, 135, 141, 162, 175]

ns['ptc22.tr2'] = [1, 7, 8, 11, 19, 21, 22, 32, 35, 39, 41, 46, 52, 68, 142, 172, 174, 177,
189, 206, 214, 219, 224, 237, 253, 262, 272, 274, 285]


m = {}
for trackname in tracknames:
    m[trackname] = np.zeros((2, 8), dtype=np.int64)

for track in tracks:
    nsresponders = ns[track.absname]
    for nid, n in track.alln.items():
        xi = 0 if nid in nsresponders else 1  # matrix rows
        spiketypeyi = spiketype2int[n.spiketype] # matrix spiketype columns
        rftypeyi = rftype2int[n.rftype] # matrix rftype columns
        m[track.absname][xi, spiketypeyi] += 1
        m[track.absname][xi, rftypeyi] += 1
        m['all'][xi, spiketypeyi] += 1
        m['all'][xi, rftypeyi] += 1

for trackname in tracknames:
    figure(figsize=(7, 2))
    imshow(m[trackname], origin='upper', cmap='gray')
    xticks(np.arange(8), ['fast', 'slow', 'fast\nasym', 'slow\nasym',
                          'simple', 'complex', 'LGN\naff', 'unknown'])
    yticks(np.arange(2), ["respond", "don't\nrespond"])
    colorbar(ticks=[m[trackname].min(), m[trackname].max()], label='neuron count')
    gcfm().window.setWindowTitle(trackname + ' cell type vs nsresponse')
    tight_layout(pad=1)

# plot again, splitting spiketype and rftype into subplots:
for trackname in tracknames:
    figure(figsize=(8, 2))
    for spi in range(2):
        a = subplot(1, 2, spi+1) # subplot indices into subplot() are 1-based
        y0, y1 = spi*4, (spi + 1)*4
        imshow(m[trackname][:, y0:y1], origin='upper', cmap='gray')
        xticks(np.arange(4), ['fast', 'slow', 'fast\nasym', 'slow\nasym',
                              'simple', 'complex', 'LGN\naff', 'unknown'][y0:y1])
        yticks(np.arange(2), ["respond", "don't\nrespond"])
        colorbar(ticks=[m[trackname][:, y0:y1].min(), m[trackname][:, y0:y1].max()],
                        label='neuron count')
    gcfm().window.setWindowTitle(trackname + ' cell type vs nsresponse split')
    tight_layout(pad=1)

print('individual chi square tests, alpha=%g' % alpha)
labels = 'fast', 'slow', 'fasym', 'sasym', 'simple', 'complex', 'LGN', 'unknown'
for label, col in zip(labels, m['all'].T):
    chi2, p = chisquare(col)
    print('%s: %r, %.2g, %.2g, %r' % (label, list(col), chi2, p, p<alpha))

print()
print("Seems the following gives the same results:")
print()

mspiketype = m['all'][:, :4]
mrftype = m['all'][:, 4:]

print('spike type group chi square tests, alpha=%g' % alpha)
print(mspiketype)
chi2, p = chisquare(mspiketype)
print('chi2:')
print(chi2)
print('p:')
print(p)
print('p<alpha')
print(p<alpha)
print()

print('RF type group chi square tests, alpha=%g' % alpha)
print(mrftype)
chi2, p = chisquare(mrftype)
print('chi2:')
print(chi2)
print('p:')
print(p)
print('p<alpha')
print(p<alpha)

show()
