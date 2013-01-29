"""Create a CodeCorr.shifts() plot for lots of recordings.
Run by calling `%run -i /path/to/shifts.py` within neuropy"""
from pylab import get_current_fig_manager as gcfm

from animal import Animal
try: ptc22; 
except NameError: ptc22 = Animal('/home/mspacek/data/ptc22')
ptc22.load('tr1')

for i in range(4, 23):
    s = 'ptc22.tr1.r%02d.cc().shifts(-12000, 12000, 50)' % i
    exec(s)
    title(s)
    gcfm().window.setWindowTitle(s)
    print(s)
    show()

ptc22.load('tr2')

for i in range(23, 37):
    s = 'ptc22.tr2.r%02d.cc().shifts(-12000, 12000, 50)' % i
    exec(s)
    title(s)
    gcfm().window.setWindowTitle(s)
    print(s)
    show()
