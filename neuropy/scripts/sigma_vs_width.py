"""Scatter plot spatial sigma vs various measures of spike width. Run from within neuropy
using `run -i scripts/sigma_vs_width.py`"""

## TODO: separate cells according to straight line in duration vs interpeak interval plot,
## should provide cleaner separation than just splitting by duration at 350 us

import scipy
from pylab import get_current_fig_manager as gcfm
from core import intround

# fractional position along waveform to assume characteristic peak is roughly aligned to:
alignf = 0.4
newtres = 1 # tres to interpolate to, in us
absslopethresh = 0.4 # uV/us
tracks = [ptc15.tr7c, ptc22.tr1, ptc22.tr2] # need to be loaded ahead of time
#tracknames = [ track.absname for track in tracks ]

def calc_t(nt, tres):
    """Generate original and desired interpolated timebases t0 and t1, as well as initial
    guess for where the alignment point is in t1"""
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

nt = 50
tres = 20
t0, t1, aligni = calc_t(nt, tres) # initial guess, for speed
sigmas = []
waves = []
dts = [] # interpeak interval
fwhms = [] # full-width half max values
slopes = []
durations = [] # spike duration, measured by time between slope threshold crossings

# collect maxchan waveforms and calculate various measures of waveform duration:
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
        # interpolate waveforms from original t0 timebase to higher rez t1 timebase:
        wave = scipy.interpolate.spline(t0, wave, t1)
        waves.append(wave)
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

# plot unclassified maxchan waveforms:
figure(figsize=(3, 3))
for wave in waves:
    plot(t1, wave, '-', lw=1) # overplot all waveforms
yticks(np.arange(-200, 200+100, 100))
xlabel('t ($\mu$s)')
ylabel('voltage ($\mu$V)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('waveforms')
tight_layout(pad=0.3)
'''
# scatter plot sigma vs dt
figure(figsize=(3, 3))
plot(dts, sigmas, '.')
xlabel('interpeak interval ($\mu$s)')
ylabel('$\sigma$ ($\mu$m)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs dt')
tight_layout(pad=0.3)

# scatter plot sigma vs fwhm
figure(figsize=(3, 3))
plot(fwhms, sigmas, '.')
xlabel('FWHM ($\mu$s)')
ylabel('$\sigma$ ($\mu$m)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs fwhm')
tight_layout(pad=0.3)

# scatter plot sigma vs slope
figure(figsize=(3, 3))
plot(slopes, sigmas, '.')
xlabel('slope ($\mu$V/$\mu$s)')
ylabel('$\sigma$ ($\mu$m)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs slope')
tight_layout(pad=0.3)
'''
# scatter plot sigma vs duration
figure(figsize=(3, 3))
plot(durations, sigmas, 'k.')
ylim(ymin=0)
xticks([0, 200, 400, 600, 800])
xlabel('duration ($\mu$s)')
ylabel('$\sigma$ ($\mu$m)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma vs duration')
tight_layout(pad=0.3)
'''
# scatter plot fwhm vs dt
figure(figsize=(3, 3))
plot(dts, fwhms, '.')
xlabel('interpeak interval ($\mu$s)')
ylabel('FWHM ($\mu$s)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('fwhm vs dt')
tight_layout(pad=0.3)

# scatter plot slope vs dt
figure(figsize=(3, 3))
plot(dts, slopes, '.')
xlabel('interpeak interval ($\mu$s)')
ylabel('slope ($\mu$V/$\mu$s)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('slope vs dt')
tight_layout(pad=0.3)
'''
# scatter plot duration vs dt
figure(figsize=(3, 3))
plot(dts, durations, 'k.')
xlim(xmin=0)
xticks([0, 100, 200, 300, 400])
yticks([0, 200, 400, 600, 800])
xlabel('interpeak interval ($\mu$s)')
ylabel('duration ($\mu$s)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('duration vs dt')
tight_layout(pad=0.3)
'''
# scatter plot duration vs slope
figure(figsize=(3, 3))
plot(slopes, durations, '.')
xlabel('slope ($\mu$V/$\mu$s)')
ylabel('duration ($\mu$s)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('duration vs slope')
tight_layout(pad=0.3)

# scatter plot slope vs fwhm
figure(figsize=(3, 3))
plot(fwhms, slopes, '.')
xlabel('FWHM ($\mu$s)')
ylabel('slope ($\mu$V/$\mu$s)')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('slope vs fwhm')
tight_layout(pad=0.3)
'''
# plot sigma distribution
figure(figsize=(3, 3))
hist(sigmas, bins=15, fc='k')
xlim(xmin=0)
xticks([0, 25, 50, 75, 100])
xlabel('$\sigma$ ($\mu$m)')
ylabel('neuron count')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('sigma distrib')
tight_layout(pad=0.3)

# plot dt distribution
figure(figsize=(3, 3))
hist(dts, bins=30, fc='k')
xlim(xmin=0)
xticks([0, 100, 200, 300, 400])
xlabel('interpeak interval ($\mu$s)')
ylabel('neuron count')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('dt distrib')
tight_layout(pad=0.3)
'''
# plot fwhm distribution
figure(figsize=(3, 3))
hist(fwhms, bins=30)
xlabel('FWHM ($\mu$s)')
ylabel('neuron count')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('fwhm distrib')
tight_layout(pad=0.3)

# plot slope distribution
figure(figsize=(3, 3))
hist(slopes, bins=30)
xlabel('slope ($\mu$V/$\mu$s)')
ylabel('neuron count')
#title('tracks: %r' % tracknames)
gcfm().window.setWindowTitle('slope distrib')
tight_layout(pad=0.3)
'''
# plot duration distribution
figure(figsize=(3, 3))
hist(durations, bins=30,fc='k')
xticks([0, 200, 400, 600, 800])
xlabel('duration ($\mu$s)')
ylabel('neuron count')
#title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('duration distrib')
tight_layout(pad=0.3)

# separate fast vs slow waveforms using trough at 350 us in duration distrib, given
# absslopethresh = 0.4:
fastis = np.asarray(durations) <= 350
slowis = np.asarray(durations) > 350
waves = np.asarray(waves)

# plot classified (fast and slow) waveforms:
figure(figsize=(3, 3))
for wave in waves[slowis]:
    plot(t1, wave, 'r-', lw=1)
for wave in waves[fastis]:
    plot(t1, wave, 'g-', lw=1)
yticks(np.arange(-200, 200+100, 100))
xlabel('t ($\mu$s)')
ylabel('voltage ($\mu$V)')
#title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('slow and fast waveforms')
tight_layout(pad=0.3)
'''
figure(figsize=(3, 3))
for wave in waves[slowis]:
    plot(t1, wave, 'r-', lw=1)
xlabel('t ($\mu$s)')
ylabel('voltage ($\mu$V)')
title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('slow waveforms')
tight_layout(pad=0.3)

figure(figsize=(3, 3))
for wave in waves[fastis]:
    plot(t1, wave, 'g-', lw=1)
xlabel('t ($\mu$s)')
ylabel('voltage ($\mu$V)')
title('tracks: %r, absslopethresh=%.1f' % (tracknames, absslopethresh))
gcfm().window.setWindowTitle('fast waveforms')
tight_layout(pad=0.3)
'''
show()
