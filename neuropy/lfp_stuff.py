# look at ptc17.44 - spont, 50 Hz and 20Hz are strong, different over cells

def lfpsta(lfp, spikes):
    ts = np.arange(-500000, 500000, lfp.tres)
    dtis = ts[[0, -1]] // lfp.tres
    nts = len(ts)
    nchans, nt = lfp.data.shape
    sta = zeros((nchans, nts), dtype=np.float64)
    tsi = lfp.ts.searchsorted(spikes)
    los = tsi+dtis[0]
    his = tsi+dtis[1] + 1
    loshis = np.column_stack((los, his))
    keeprows = (loshis[:, 0] > 0) * (loshis[:, 1] > 0) # ditch rows with -ve indices
    keeprows *= (loshis[:, 0] < nt) * (loshis[:, 1] < nt) # ditch rows with indices > nt
    loshis = loshis[keeprows]
    for lo, hi in loshis:
        sta += lfp.data[:, lo:hi]
    sta /= len(spikes)
    return ts, sta

def plot_lfpsta(lfp, spikes):
    ts, sta = lfpsta(lfp, spikes)
    figure()
    for chani in range(len(sta)):
        plot(ts//1000, sta[chani], ms=3)

def load_lfp(lfpfname='/data/ptc15/tr7c/92 - track 7c mseq32 0.4deg/92_-_track_7c_mseq32_0.4deg.lfp'):
    from spyke.core import WaveForm
    f = np.load(lfpfname)
    data = f['data']
    chans = f['chans']
    ts = np.arange(f['t0'], f['tend'], f['tres'])
    wave = WaveForm(data=data, ts=ts, chans=chans)
    wave.tres = f['tres']
    return wave

def get_fname():
    dlg = wx.FileDialog(None, message="Open file",
                        defaultDir=os.getcwd(), defaultFile='',
                        wildcard="*.*",
                        style=wx.OPEN)
    if dlg.ShowModal() == wx.ID_OK:
        fname = dlg.GetPath()
        head, tail = os.path.split(fname)
        os.chdir(head) # update cwd
    dlg.Destroy()
    return fname


lfp = load_lfp()
