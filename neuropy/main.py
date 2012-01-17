"""Main neuropy window"""

from __future__ import division

__authors__ = ['Martin Spacek']

import os
import sys
import platform

from IPython import embed
from IPython.core import ultratb
from IPython.frontend.terminal.ipapp import load_default_config
# has to come before Qt imports:
#from IPython.frontend.qt.console.ipython_widget import IPythonWidget
from IPython.frontend.qt.console.rich_ipython_widget import RichIPythonWidget
from IPython.frontend.qt.kernelmanager import QtKernelManager

from PyQt4 import QtCore, QtGui, uic
#from PyQt4.QtGui import QFont
#from PyQt4.QtCore import Qt
NeuropyUi, NeuropyUiBase = uic.loadUiType('neuropy.ui')

from __init__ import __version__
from core import DATAPATH
from recording import Recording


class NeuropyWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        """Note that a lot of default options can be changed by editing
        IPython.frontend.qt.console.console_widget.py and ipython_widget.py"""        
        QtGui.QMainWindow.__init__(self)
        self.ui = NeuropyUi()
        self.ui.setupUi(self) # lay it out

        # might need local_kernel=True kwarg
        # can set paging: 'inside', 'vsplit', 'hsplit', 'none'
        #ipyqtwidget = IPythonWidget(parent=self, local_kernel=True)
        ipyqtwidget = RichIPythonWidget(parent=self, local_kernel=True) # necessary for inline
        self.ipyqtwidget = ipyqtwidget
        config = load_default_config()
        # this doesn't seem to work as intended:
        config.InteractiveShellApp.exec_lines.append('from recording import Recording')
        self.ipyqtwidget.config = config
        
        kernel_manager = QtKernelManager()
        kernel_manager.start_kernel()
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
        self.path = DATAPATH

        ## TODO: get neuropy imports to work properly without messing up ipython display header
        ## need to be done after window and ipyqtwidget init are finished


    @QtCore.pyqtSlot()
    def on_actionOpen_triggered(self):
        getExistingDirectory = QtGui.QFileDialog.getExistingDirectory
        path = getExistingDirectory(self, caption="Open recording, track, or animal",
                                    directory=self.path)
        path = str(path)
        if path == '':
            return
        self.path = os.path.normpath(os.path.join(path, os.pardir)) # update with path's parent
        self.open_recording(path)

    @QtCore.pyqtSlot()
    def on_action_ptc15_r87_triggered(self):
        path = os.path.join(DATAPATH, 'ptc15/tr7c/87 - track 7c spontaneous craziness')
        self.open_recording(path)

    @QtCore.pyqtSlot()
    def on_action_ptc15_r92_triggered(self):
        path = os.path.join(DATAPATH, 'ptc15/tr7c/92 - track 7c mseq32 0.4deg')
        self.open_recording(path)

    @QtCore.pyqtSlot()
    def on_action_ptc17_r03_triggered(self):
        path = os.path.join(DATAPATH, 'ptc17/tr1/03-tr1-mseq32_40ms')
        self.open_recording(path)

    def open_recording(self, path):
        rec = Recording(path) # init it just to parse its id
        exec_lines = (
        "from recording import Recording\n"
        "r%d = Recording(%r)\n"
        "r%d.load()" % (rec.id, path, rec.id)
        )
        self.ipyqtwidget.execute(exec_lines)

    @QtCore.pyqtSlot()
    def on_actionShell_triggered(self):
        embed(display_banner=False) # "self" is accessible
        # embed() seems to override the excepthook, need to reset it:
        set_excepthook()

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
