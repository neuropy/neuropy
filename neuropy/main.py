"""Main neuropy window"""

from __future__ import division
from __init__ import __version__

__authors__ = ['Martin Spacek']

import os
import sys
import platform

from IPython import embed
from IPython.core import ultratb
from IPython.frontend.terminal.ipapp import load_default_config
# has to come before Qt imports:
from IPython.frontend.qt.console.ipython_widget import IPythonWidget

from PyQt4 import QtCore, QtGui, uic
#from PyQt4.QtCore import Qt
NeuropyUi, NeuropyUiBase = uic.loadUiType('neuropy.ui')


class NeuropyWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self)
        self.ui = NeuropyUi()
        self.ui.setupUi(self) # lay it out

        ipyqtwidget = IPythonWidget()
        self.setCentralWidget(ipyqtwidget)

        #self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('neuropy')

    @QtCore.pyqtSlot()
    def on_actionOpen_triggered(self):
        getExistingDirectory = QtGui.QFileDialog.getExistingDirectory
        path = getExistingDirectory(self, caption="Open recording, track, or animal",
                                    directory=os.getcwd())
        path = str(path)
        print(path)

    @QtCore.pyqtSlot()
    def on_actionQuit_triggered(self):
        self.close() # call close() before destroy() to avoid segfault
        self.destroy()

    @QtCore.pyqtSlot()
    def on_actionAboutNeuropy_triggered(self):
        text = """
        <h2>neuropy %s</h2>
        <p>A tool for neuronal spike data analysis</p>
        <p>Copyright &copy; 2006-2011 Martin Spacek<br>
           University of British Columbia</p>
        <p>Python %s, Qt %s, PyQt %s<br>
        %s</p>""" % (__version__, platform.python_version(),
        QtCore.QT_VERSION_STR, QtCore.PYQT_VERSION_STR, platform.platform())
        QtGui.QMessageBox.about(self, "About neuropy", text)

    @QtCore.pyqtSlot()
    def on_actionAboutQt_triggered(self):
        QtGui.QMessageBox.aboutQt(self)
    '''
    def shell(self):
        embed(display_banner=False, config=load_default_config()) # "self" is accessible
        # embed() seems to override the excepthook, need to reset it:
        set_excepthook()
    '''
    def raise_error(self):
        raise RuntimeError


def set_excepthook():
    """Drops us into IPython's debugger on any error"""
    sys.excepthook = ultratb.FormattedTB(mode='Verbose', call_pdb=1)


if __name__ == "__main__":
    # prevents "The event loop is already running" errors:
    QtCore.pyqtRemoveInputHook()
    set_excepthook()
    app = QtGui.QApplication(sys.argv)
    neuropywindow = NeuropyWindow()
    neuropywindow.show()
    sys.exit(app.exec_())
