"""Plot 2D matrix of rftype classification vs spiketype. Run from within neuropy using `run -i
scripts/cell_type_matrix.py`"""

from __future__ import division

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time

spiketype2int = {'fast':0, 'slow':1, 'fastasym':2, 'slowasym':3}
rftype2int = {'simple':0, 'complex':1, 'LGN':2, None: 3}

m = np.zeros((4, 4), dtype=np.int64)

for track in tracks:
    for nid, n in track.alln.items():
        xi = rftype2int[n.rftype] # matrix rows
        yi = spiketype2int[n.spiketype] # matrix columns
        m[xi, yi] += 1

figure(figsize=(3.25, 3))
imshow(m, origin='upper', cmap='gray')
xticks(np.arange(4), ['fast', 'slow', 'fast asym', 'slow asym'], rotation=90)
yticks(np.arange(4), ['simple', 'complex', 'LGN aff', 'unknown'])
colorbar(ticks=[0, m.max()], label='neuron count')
gcfm().window.setWindowTitle('rftype vs spiketype')
tight_layout(pad=1)

show()
