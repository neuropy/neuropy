"""
test_window.py

run 'ipython-qtconsole --pylab=qt',
or 'ipython --pylab=qt',

then type:

from test_window import test, Container
c = Container()
test(c)

"""
from PyQt4 import QtGui

class TestWindow(QtGui.QMainWindow):
    def __init__(self, parent=None, n=0):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle('Test window %d' %n)

class Container(object):
    def __init__(self):
        pass
        
def test(c):
    testwin1 = TestWindow(n=1)
    testwin1.show()
    testwin2 = TestWindow(n=2)
    testwin2.show()
    c.testwin1 = testwin1
    return testwin2
    # both testwin1 and testwin2 will be shown
