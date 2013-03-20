"""Shows how to use multidimensional scaling to take a N*(N-1)/2 pairwise values,
like spike correlations, and get a 1D projection which ostensibly allows sorting
N units such that those pairs with the greatest similarity are kept as close
to each other as possible, with as little "strain" as possible. This might then
be useful for sorting raster plots to better reveal ensemble activity.

Try messing with kwargs to MDS:

n_components: number of output dimensions
metric: (True, False) compute metric or nonmetric MDS
dissimilarity: ('precomputed', 'euclidean')
"""

import numpy as np
from sklearn.manifold import MDS

cc = ptc22.tr1.r08.cc()
cc.calc()
sim = cc.corrs # 1D array of length npairs
# normalize, might not be necessary, although -ve values might be bad:
sim -= sim.min()
sim /= sim.max()
N = len(cc.nids)
ui = np.triu_indices(N, 1) # upper triangle indices
li = np.tril_indices(N, -1) # lower triangle indices
simm = np.eye(N) # NxN matrix, 1s on diagonal
simm[ui] = sim # fill upper triangle
simm[li] = simm.T[li] # make symmetric by filling lower triangle

mds = MDS(n_components=1, metric=False, dissimilarity="precomputed")
pos = mds.fit_transform(simm)
sortis = pos.ravel().argsort()
print(sortis)
