"""Filtering functions"""

import numpy as np
import scipy.signal


def notch(data, sampfreq=1000, freq=60, bw=0.25, gpass=0.01, gstop=30, ftype='ellip'):
    """Filter out frequencies in data centered on freq (Hz), of bandwidth +/- bw (Hz) on
    data row indices chanis.

    ftype: 'ellip', 'butter', 'cheby1', 'cheby2', 'bessel'
    """
    w = freq / (sampfreq / 2) # fraction of Nyquist frequency == 1/2 sampling rate
    bw = bw / (sampfreq / 2)
    wp = [w-2*bw, w+2*bw] # outer bandpass
    ws = [w-bw, w+bw] # inner bandstop
    # using more extreme values for gpass or gstop seems to cause IIR filter instability.
    # 'ellip' is the only one that seems to work
    b, a = scipy.signal.iirdesign(wp, ws, gpass=gpass, gstop=gstop, analog=0, ftype=ftype)
    data = scipy.signal.lfilter(b, a, data)
    return data, b, a

def naivenotch(data, sampfreq=1000, freqs=60, bws=1):
    """Filter out frequencies in data centered on freqs (Hz), of bandwidths bws (Hz).
    Filtering out by setting components to 0 is probably naive"""
    nt = data.shape[1]
    tres = 1 / sampfreq
    dt = tres / 1e6 # in sec
    f = np.fft.fftfreq(nt, dt) # fft bin frequencies
    f = f[:nt//2] # grab +ve freqs by splitting f in half
    franges = []
    freqs = tolist(freqs)
    bws = tolist(bws)
    if len(freqs) > 1 and len(bws) == 1:
        bws = bws * len(freqs) # make freqs and bw the same length
    for freq, bw in zip(freqs, bws):
        franges.append(freq-bw)
        franges.append(freq+bw)
    fis = f.searchsorted(franges)
    fis = np.hstack([fis, -fis]) # indices for both +ve and -ve freq ranges
    fis.shape = -1, 2 # reshape to 2 columns
    fdata = np.fft.fft(data)
    for f0i, f1i in fis:
        fdata[:, f0i:f1i] = 0 # replace desired components with 0
        # maybe try using complex average of freq bins just outside of freqs +/- bws
    data = np.fft.ifft(fdata).real # inverse FFT, leave as float
    return data

def bandpass(data, sampfreq=1000, f0=0, f1=7, fr=0.5, gpass=0.01, gstop=30, ftype='ellip'):
    """Bandpass filter data on row indices chanis, between f0 and f1 (Hz), with
    filter rolloff (?) fr (Hz).

    ftype: 'ellip', 'butter', 'cheby1', 'cheby2', 'bessel'
    """
    w0 = f0 / (sampfreq / 2) # fraction of Nyquist frequency == 1/2 sampling rate
    w1 = f1 / (sampfreq / 2)
    wr = fr / (sampfreq / 2)
    if w0 == 0:
        wp = w1
        ws = w1+wr
    elif w1 == 0:
        wp = w0
        ws = w0-wr
    else:
        wp = [w0, w1]
        ws = [w0-wr, w1+wr]
    b, a = scipy.signal.iirdesign(wp, ws, gpass=gpass, gstop=gstop, analog=0, ftype=ftype)
    data = scipy.signal.lfilter(b, a, data)
    return data, b, a

def filter(data, sampfreq=1000, f0=300, f1=None, order=4, rp=None, rs=None,
           btype='highpass', ftype='butter'):
    """Filter data by specifying filter order and btype, instead of gpass and gstop

    ftype: 'ellip', 'butter', 'cheby1', 'cheby2', 'bessel'

    For 'ellip', need to also specify passband and stopband ripple with rp and rs.
    """
    if f1 != None:
        fn = np.array([f0, f1])
    else:
        fn = f0
    wn = fn / (sampfreq / 2)
    b, a = scipy.signal.iirfilter(order, wn, rp=None, rs=None, btype=btype, analog=0,
                                  ftype=ftype, output='ba')
    data = scipy.signal.lfilter(b, a, data)
    return data, b, a

def hilbert(x):
    """Return power (dB wrt 1 mV) and phase (rad) of Hilbert transform of data in x"""
    hx = scipy.signal.hilbert(x) # Hilbert transform of x
    Ex = np.abs(hx) # amplitude == energy?
    Phx = np.angle(hx) # phase
    Px = 10 * np.log(Ex**2) # power in dB wrt 1 mV^2?
    return Px, Phx
