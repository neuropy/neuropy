"""Plot distribution of spatial sigmas of designated tracks"""

tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time

sigmas = []
for track in tracks:
    for n in track.alln.values():
        sigmas.append(n.sigma)

hist(sigmas, bins=20)
