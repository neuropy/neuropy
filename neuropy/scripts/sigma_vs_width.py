"""Scatter plot spatial sigma vs inter-peak interval, or FWHM of template characteristic peak"""

import scipy
from pylab import get_current_fig_manager as gcfm
from core import intround

def calc_t(nt, tres):
    tend = nt*tres
    t0 = np.arange(0, tend, tres)
    t1 = np.arange(0, tend-tres, newtres)
    aligni = intround(alignf * len(t1))
    return t0, t1, aligni

def argextrema(a):
    """Return indices of local extrema in 1D array a.
    Taken from http://stackoverflow.com/a/9667121/2020363"""
    return np.diff(np.sign(np.diff(a))).nonzero()[0] + 1

def argfwhm(a, exti, fraction=0.5):
    """Find timepoints of full width half max (or whatever fraction is) around extremum exti
    in 1D array a"""
    a = abs(a)
    fm = a[exti] * fraction # fraction of max
    d = a - fm
    lis = np.diff(np.sign(d[:exti])).nonzero()[0]
    ris = np.diff(np.sign(d[exti:])).nonzero()[0] + exti + 1
    assert len(lis) > 0
    assert len(ris) > 0
    # return rightmost of left indices, and leftmost of right indices:
    return lis[-1], ris[0]

# fractional position along waveform to assume characteristic peak is roughly aligned to:
alignf = 0.4
newtres = 1 # tres to interpolate to, in us
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
tracknames = [ track.absname for track in tracks ]

nt = 50
tres = 20
t0, t1, aligni = calc_t(nt, tres) # initial guess, for speed
sigmas = []
dts = []
fwhms = []
figure() # init fig for waveform plots
for track in tracks:
    if tres != track.tres: # recalculate timepoints
        tres = track.tres
        t0, t1, aligni = calc_t(nt, tres)
    nids = sorted(track.alln)
    for nid in nids:
        n = track.alln[nid]
        if nt != n.nt: # recalculate timepoints
            nt = n.nt
            t0, t1, aligni = calc_t(nt, tres)
        sigmas.append(n.sigma)
        maxchani = n.chans.searchsorted(n.maxchan)
        wave = n.wavedata[maxchani]
        # interpolate waveforms for higher rez dt:
        wave = scipy.interpolate.spline(t0, wave, t1)
        plot(wave, '.', ms=2) # overplot all waveforms
        t0i, t1i = wave.argmax(), wave.argmin()
        dt = abs(t0i - t1i) * newtres
        dts.append(dt) # dt between biggest max and min peaks
        extis = argextrema(wave) # indices of all local extrema
        exti = extis[abs(extis - aligni).argmin()] # extremum closest to aligni
        li, ri = argfwhm(wave, exti, fraction=0.5)
        fwhm = abs(li - ri) * newtres
        fwhms.append(fwhm)

# label the wave plots
xlabel('t (us)')
ylabel('voltage (uV)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('waveforms')
tight_layout(pad=0.3)

# scatter plot sigma vs dt
figure()
plot(dts, sigmas, '.')
xlabel('dt (us)')
ylabel('sigma (um)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs dt')
tight_layout(pad=0.3)

# scatter plot sigma vs fwhm
figure()
plot(fwhms, sigmas, '.')
xlabel('FWHM (us)')
ylabel('sigma (um)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs fwhm')
tight_layout(pad=0.3)

# scatter plot fwhm vs dt
figure()
plot(dts, fwhms, '.')
xlabel('dt (us)')
ylabel('FWHM (us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('fwhm vs dt')
tight_layout(pad=0.3)

# plot sigma distribution
figure()
hist(sigmas, bins=30) # bins=40 looks a little more tantalizing
xlabel('sigma (um)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma distrib')
tight_layout(pad=0.3)

# plot dt distribution
figure()
hist(dts, bins=30)
xlabel('dt (us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('dt distrib')
tight_layout(pad=0.3)

# plot fwhm distribution
figure()
hist(fwhms, bins=30)
xlabel('FWHM (us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('fwhm distrib')
tight_layout(pad=0.3)
