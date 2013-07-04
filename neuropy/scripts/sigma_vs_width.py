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
absslopethresh = 0.4 # uV/us
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
tracknames = [ track.absname for track in tracks ]

nt = 50
tres = 20
t0, t1, aligni = calc_t(nt, tres) # initial guess, for speed
sigmas = []
waves = []
dts = []
fwhms = []
slopes = []
durations = []
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
        waves.append(wave)
        plot(t1, wave, '-', lw=1) # overplot all waveforms
        t0i, t1i = wave.argmax(), wave.argmin()
        dt = abs(t0i - t1i) * newtres
        dts.append(dt) # dt between biggest max and min peaks
        extis = argextrema(wave) # indices of all local extrema
        exti = extis[abs(extis - aligni).argmin()] # extremum closest to aligni
        li, ri = argfwhm(wave, exti, fraction=0.5)
        fwhm = abs(li - ri) * newtres
        fwhms.append(fwhm)
        absslope = abs(np.diff(wave)) / newtres # uV/us
        maxabsslope = max(absslope)
        slopes.append(maxabsslope)
        # another way to measure waveform duration is to see over what duration the abs(slope)
        # is greater than something close to 0. Starting from each end, at what timepoint does
        # the slope exceed this minimum threshold? Difference between timepoints is duration
        # of waveform
        flatis = np.where(absslope > absslopethresh)[0]
        if len(flatis) < 2:
            # exclude cells whose slope isn't above threshold for at least two timepoints:
            duration = 0
        else:
            duration = (flatis[-1] - flatis[0]) * newtres
        durations.append(duration)        
        # or as an alternative to using absslopethresh, measure the fwhm of the last
        # extremum in each waveform. Looking at the overplotted waveforms, that, strangely,
        # is where I see the most dichotomy. Not so much in the first peak.

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

# scatter plot sigma vs slope
figure()
plot(slopes, sigmas, '.')
xlabel('slope (uv/us)')
ylabel('sigma (um)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs slope')
tight_layout(pad=0.3)

# scatter plot sigma vs duration
figure()
plot(durations, sigmas, '.')
xlabel('duration (us)')
ylabel('sigma (um)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs duration')
tight_layout(pad=0.3)

# scatter plot fwhm vs dt
figure()
plot(dts, fwhms, '.')
xlabel('dt (us)')
ylabel('FWHM (us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('fwhm vs dt')
tight_layout(pad=0.3)

# scatter plot slope vs dt
figure()
plot(dts, slopes, '.')
xlabel('dt (us)')
ylabel('slope (uv/us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('slope vs dt')
tight_layout(pad=0.3)

# scatter plot duration vs dt
figure()
plot(dts, durations, '.')
xlabel('dt (us)')
ylabel('duration (us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('duration vs dt')
tight_layout(pad=0.3)

# scatter plot duration vs slope
figure()
plot(slopes, durations, '.')
xlabel('slope (uV/us)')
ylabel('duration (us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('duration vs slope')
tight_layout(pad=0.3)

# scatter plot slope vs fwhm
figure()
plot(fwhms, slopes, '.')
xlabel('FWHM (us)')
ylabel('slope (uv/us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('slope vs fwhm')
tight_layout(pad=0.3)

# plot sigma distribution
figure()
hist(sigmas, bins=30)
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

# plot slope distribution
figure()
hist(slopes, bins=30)
xlabel('slope (uV/us)')
title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('slope distrib')
tight_layout(pad=0.3)

# plot duration distribution
figure()
hist(durations, bins=30)
xlabel('duration (us)')
title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('duration distrib')
tight_layout(pad=0.3)

# separate fast vs slow waveforms using trough at 350 us in duration distrib, given
# absslopethresh = 0.4:
fastis = np.asarray(durations) <= 350
slowis = np.asarray(durations) > 350
waves = np.asarray(waves)

figure()
for wave in waves[slowis]:
    plot(t1, wave, 'r-', lw=1)
for wave in waves[fastis]:
    plot(t1, wave, 'g-', lw=1)
xlabel('t (us)')
ylabel('voltage (uV)')
title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('slow and fast waveforms')
tight_layout(pad=0.3)
'''
figure()
for wave in waves[slowis]:
    plot(t1, wave, 'r-', lw=1)
xlabel('t (us)')
ylabel('voltage (uV)')
title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('slow waveforms')
tight_layout(pad=0.3)

figure()
for wave in waves[fastis]:
    plot(t1, wave, 'g-', lw=1)
xlabel('t (us)')
ylabel('voltage (uV)')
title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('fast waveforms')
tight_layout(pad=0.3)
'''
