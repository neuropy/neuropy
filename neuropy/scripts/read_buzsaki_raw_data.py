"""Read raw .dat data file from Buzsaki (from Julia Farms spike sorting meeting) and treat it
as an LFP object. Experiment with different filtering settings"""

from core import LFP

class BZData(LFP):
    def __init__(self, fname='/home/mspacek/work/Buzsaki_raw_data/trace_8Chan_High-Sleep.dat'):
        LFP.__init__(self, Recording(''), fname) # give it a fake recording

    def load(self):
        self.uVperAD = 0.1 # blind guess
        data = np.fromfile(self.fname, dtype=np.int16)
        data = data * self.uVperAD # convert to float uV
        data.shape = -1, 8 # 8 chans, chans are changing fastest in file, 20 kHz sampled data
        data = data.T # transpose to 8 rows
        self.data = data
        nt = data.shape[1]
        self.sampfreq = 20000
        self.tres = 50 # 50 us per sample at 20 kHz
        # fake chanpos and chans, there are correct values, but I dont know what they are:
        self.chanpos = array([[   0,    0],
                              [   0,  100],
                              [   0,  200],
                              [   0,  300],
                              [   0,  400],
                              [   0,  500],
                              [   0,  600],
                              [   0,  700],
                                         ])
        self.chans = np.arange(8)                             
        self.t0 = 0
        self.t1 = (nt-1) * self.tres
        self.PLOTGAIN = 20

    def specgram(self, t0=None, t1=None, f0=None, f1=2000, p0=None, p1=None, chanis=-1,
                 width=2**16, overlap=2**15, cm=None, colorbar=False, figsize=(20, 6.5)):
        LFP.specgram(self, t0, t1, f0, f1, p0, p1, chanis, width, overlap, cm, colorbar,
                     figsize)

    def bandpass(self, chanis=-1, f0=500, f1=0, fr=100, gpass=0.01, gstop=60, ftype='ellip'):
        LFP.bandpass(self, chanis, f0, f1, fr, gpass, gstop, ftype)


bz = BZData()
bz.load()

bz.plot(0, 1)
bz.specgram(0, 500, p0=None, p1=None, f1=2000)
bz.bandpass(chanis=-1, f0=500)
