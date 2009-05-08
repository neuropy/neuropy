"""A (very incomplete) test suite for neuropy"""

from Core import *

def test():

    # check approx()
    t = arange(-1,2)
    ts = arange(-1,2,0.1)
    assert approx(t,ts[[0,10,20]]).all()

    # check sah, make sure it keeps existing data points if some of the new timepoints are the same as the old ones
    r = rand(3)
    rs = sah(t, r, ts, keep=True)
    plot(t, r, 'b.')
    plot(ts, rs, 'g+')
