"""Plot 2D matrix of rftype classification vs spiketype. Run from within neuropy using `run -i
scripts/cell_type_matrix.py`"""

from __future__ import division

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
tracknames = [track.absname for track in tracks ] + ['all']

spiketype2int = {'fast':0, 'slow':1, 'fastasym':2, 'slowasym':3}
rftype2int = {'simple':0, 'complex':1, 'LGN':2, None: 3}

m = {}
for trackname in tracknames:
    m[trackname] = np.zeros((4, 4), dtype=np.int64)

for track in tracks:
    for nid, n in track.alln.items():
        xi = rftype2int[n.rftype] # matrix rows
        yi = spiketype2int[n.spiketype] # matrix columns
        m[track.absname][xi, yi] += 1
        m['all'][xi, yi] += 1

for trackname in tracknames:
    figure(figsize=(3.25, 3))
    imshow(m[trackname], origin='upper', cmap='gray')
    xticks(np.arange(4), ['fast', 'slow', 'fast asym', 'slow asym'], rotation=90)
    yticks(np.arange(4), ['simple', 'complex', 'LGN aff', 'unknown'])
    colorbar(ticks=[m[trackname].min(), m[trackname].max()], label='neuron count')
    gcfm().window.setWindowTitle(trackname + ' rftype vs spiketype')
    tight_layout(pad=1)

show()
