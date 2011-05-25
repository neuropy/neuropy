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
from IPython.frontend.qt.kernelmanager import QtKernelManager

from PyQt4 import QtCore, QtGui, uic
#from PyQt4.QtGui import QFont
#from PyQt4.QtCore import Qt
NeuropyUi, NeuropyUiBase = uic.loadUiType('neuropy.ui')


class NeuropyWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        """Note that a lot of default options can be changed by editing
        IPython.frontend.qt.console.console_widget.py and ipython_widget.py"""        
        QtGui.QMainWindow.__init__(self)
        self.ui = NeuropyUi()
        self.ui.setupUi(self) # lay it out

        # might need local_kernel=True kwarg
        # can set paging style, like 'vsplit', 'hsplit', 'none', default is 'inside'
        ipyqtwidget = IPythonWidget(parent=self, local_kernel=True)
        self.ipyqtwidget = ipyqtwidget
        ipyqtwidget.config = load_default_config() # doesn't seem to work
        kernel_manager = QtKernelManager()
        kernel_manager.start_kernel() # might need **kwargs
        kernel_manager.start_channels()
        ipyqtwidget.kernel_manager = kernel_manager
        ipyqtwidget.gui_completion = False
        ipyqtwidget.set_default_style(colors='linux')
        #ipyqtwidget.font = QFont('Lucida Console', 12) # 3rd arg can be e.g. QFont.Bold
        #ipyqtwidget.font.setFixedPitch(True)
        # make "exit" and "quit" typed in ipyqtwidget close self
        ipyqtwidget.exit_requested.connect(self.close)
        self.setCentralWidget(ipyqtwidget)

        # ip.ex() lets you execute strings in user namespace, where ip is the IPython.zmq.zmqshell.ZMQInteractiveShell
        # ip.ev() evaluates string
        # ip.write() writes a string to the shell
        # ip.push() pushes a variables in a dict to user namespace

        # to communicate with the ipy kernel user namespace, probably need to
        # use kernel_manager.kernel.communicate and send a signal somehow
        # this might be too low level though:
        #self.ipyqtwidget.kernel_manager.kernel.communicate('sdfsd')

        self.setGeometry(0, 0, 960, 1080)
        self.setWindowTitle('neuropy')

    @QtCore.pyqtSlot()
    def on_actionOpen_triggered(self):
        getExistingDirectory = QtGui.QFileDialog.getExistingDirectory
        path = getExistingDirectory(self, caption="Open recording, track, or animal",
                                    directory=os.getcwd())
        path = str(path)
        print(path)
        self.ipyqtwidget.execute("from neuropy import *\n"
                                 "r03 = Recording('03')\n"
                                 "r03.load()")

    @QtCore.pyqtSlot()
    def on_actionShell_triggered(self):
        embed(display_banner=False, config=load_default_config()) # "self" is accessible
        # embed() seems to override the excepthook, need to reset it:
        set_excepthook()

    def raise_error(self):
        raise RuntimeError

    @QtCore.pyqtSlot()
    def on_actionQuit_triggered(self):
        self.close() # call close() before destroy() to avoid segfault
        self.destroy()

    def closeEvent(self, event):
        km = self.ipyqtwidget.kernel_manager
        if km and km.channels_running:
            km.shutdown_kernel()
            event.accept()

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


def set_excepthook():
    """Drops us into IPython's debugger on any error"""
    sys.excepthook = ultratb.FormattedTB(mode='Verbose', call_pdb=1)


if __name__ == "__main__":
    # prevents "The event loop is already running" errors when dropping into shell:
    QtCore.pyqtRemoveInputHook()
    # this causes ipython shutdown errors on close. Maybe need to restore InputHook somewhere:
    #set_excepthook()
    app = QtGui.QApplication(sys.argv)
    neuropywindow = NeuropyWindow()
    neuropywindow.show()
    sys.exit(app.exec_())
