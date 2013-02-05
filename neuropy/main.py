"""Main neuropy window"""

from __future__ import division

__authors__ = ['Martin Spacek']

import os
import sys
import platform

from IPython import embed
from IPython.core import ultratb
# has to come before Qt imports:
from IPython.frontend.qt.console.rich_ipython_widget import RichIPythonWidget
from IPython.lib import guisupport
from IPython.utils.path import get_ipython_dir
from IPython.config.loader import PyFileConfigLoader

from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtGui import QFont
#from PyQt4.QtCore import Qt
NeuropyUi, NeuropyUiBase = uic.loadUiType('neuropy.ui')

from __init__ import __version__
from animal import Animal
from track import Track
from recording import Recording

DATAPATH = os.path.expanduser('~/data')
INPROCESS = False # use inprocess kernel? otherwise, use zmq kernel


class NeuropyWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        """Note that a lot of default options can be changed by editing
        IPython.frontend.qt.console.console_widget.py and ipython_widget.py"""
        QtGui.QMainWindow.__init__(self)
        self.ui = NeuropyUi()
        self.ui.setupUi(self) # lay it out

        ipw = RichIPythonWidget(parent=self)
        self.ipw = ipw
        self.setCentralWidget(ipw)

        self.setGeometry(0, 0, 960, 1080)
        self.setWindowTitle('neuropy')
        self.path = DATAPATH

    @QtCore.pyqtSlot()
    def on_actionOpen_triggered(self):
        getExistingDirectory = QtGui.QFileDialog.getExistingDirectory
        path = getExistingDirectory(self, caption="Open recording, track, or animal",
                                    directory=self.path)
        path = str(path)
        if path == '':
            return
        self.path = os.path.normpath(os.path.join(path, os.pardir)) # update with path's parent

        head, tail = os.path.split(path)
        if tail.startswith('pt'): # assume it's an animal
            self.open_animal(path)
        elif tail.startswith('tr'): # assume it's a track
            self.open_track(path)
        else: # asume it's a recording
            self.open_recording(path)

    @QtCore.pyqtSlot()
    def on_action_ptc15_triggered(self):
        path = os.path.join(DATAPATH, 'ptc15')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc17_triggered(self):
        path = os.path.join(DATAPATH, 'ptc17')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc18_triggered(self):
        path = os.path.join(DATAPATH, 'ptc18')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc20_triggered(self):
        path = os.path.join(DATAPATH, 'ptc20')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc21_triggered(self):
        path = os.path.join(DATAPATH, 'ptc21')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc22_triggered(self):
        path = os.path.join(DATAPATH, 'ptc22')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc15_tr7c_triggered(self):
        path = os.path.join(DATAPATH, 'ptc15')
        self.open_animal(path, 'tr7c')

    @QtCore.pyqtSlot()
    def on_action_ptc18_tr1_triggered(self):
        path = os.path.join(DATAPATH, 'ptc18')
        self.open_animal(path, 'tr1')

    @QtCore.pyqtSlot()
    def on_action_ptc22_tr1_triggered(self):
        path = os.path.join(DATAPATH, 'ptc22')
        self.open_animal(path, 'tr1')

    @QtCore.pyqtSlot()
    def on_action_ptc22_tr2_triggered(self):
        path = os.path.join(DATAPATH, 'ptc22')
        self.open_animal(path, 'tr2')

    def open_animal(self, path, tracknames=None):
        a = Animal(path) # init it just to parse its name
        exec_lines = (
        "try: %s; \n"
        "except NameError: %s = Animal(%r)\n"
        "%s.load(%r)" % (a.name, a.name, path, a.name, tracknames)
        )
        self.ipw.execute(exec_lines)

    def open_track(self, path):
        tr = Track(path) # init it just to parse its id
        exec_lines = (
        "tr%d = Track(%r)\n"
        "tr%d.load()" % (tr.id, path, tr.id)
        )
        self.ipw.execute(exec_lines)

    def open_recording(self, path):
        rec = Recording(path) # init it just to parse its id
        exec_lines = (
        "r%d = Recording(%r)\n"
        "r%d.load()" % (rec.id, path, rec.id)
        )
        self.ipw.execute(exec_lines)

    @QtCore.pyqtSlot()
    def on_actionShell_triggered(self):
        ## TODO: this raises an error in IPython 0.14.dev:
        embed()

    @QtCore.pyqtSlot()
    def on_actionQuit_triggered(self):
        self.close()
        #self.destroy() # no longer seems necessary, may cause segfaults?

    @QtCore.pyqtSlot()
    def on_actionAboutNeuropy_triggered(self):
        text = """
        <h2>neuropy %s</h2>
        <p>A tool for neuronal spike data analysis</p>
        <p>Copyright &copy; 2006-2013 Martin Spacek<br>
           University of British Columbia</p>
        <p>Python %s, Qt %s, PyQt %s<br>
        %s</p>""" % (__version__, platform.python_version(),
        QtCore.QT_VERSION_STR, QtCore.PYQT_VERSION_STR, platform.platform())
        QtGui.QMessageBox.about(self, "About neuropy", text)

    @QtCore.pyqtSlot()
    def on_actionAboutQt_triggered(self):
        QtGui.QMessageBox.aboutQt(self)


def config_ipw(ipw):
    """Modify some default settings of IPython Qt widget"""
    ipython_dir = get_ipython_dir()
    profile_dir = os.path.join(ipython_dir, 'profile_default')
    cl = PyFileConfigLoader('ipython_qtconsole_config.py', profile_dir)
    config = cl.load_config()
    ipw.config = config

    ipw.gui_completion = 'ncurses' # 'plain, 'droplist' or 'ncurses'
    ipw.set_default_style(colors='Linux')
    ipw.font = QFont('Lucida Console', 11) # 3rd arg can be e.g. QFont.Bold
    ipw.font.setFixedPitch(True)

def main():
    """Start kernel manager, create window, run app event loop, auto execute some code
    in user namespace. Adapted from IPython example in:
    docs/examples/frontend/inprocess_qtconsole.py"""
    app = guisupport.get_app_qt4()

    if INPROCESS:
        from IPython.kernel.inprocess.ipkernel import InProcessKernel
        from IPython.frontend.qt.inprocess_kernelmanager import QtInProcessKernelManager
        kernel = InProcessKernel(gui='qt4')
        km = QtInProcessKernelManager(kernel=kernel)
        kernel.frontends.append(km)
    else:
        from IPython.frontend.qt.kernelmanager import QtKernelManager
        km = QtKernelManager()
        km.start_kernel()
    km.start_channels()

    neuropywindow = NeuropyWindow()
    ipw = neuropywindow.ipw
    config_ipw(ipw)
    ipw.exit_requested.connect(app.quit)
    ipw.kernel_manager = km
    neuropywindow.show()

    # execute some code directly, note the output appears at the system command line:
    #kernel.shell.run_cell('print "x=%r, y=%r, z=%r" % (x,y,z)')
    # execute some code through the frontend (once the event loop is
    # running). The output appears in the ipw:
    do_later(ipw.execute_file, 'startup.py', hidden=True)
    do_later(ipw.execute_file, 'globals.py', hidden=True)

    guisupport.start_event_loop_qt4(app)

def do_later(func, *args, **kwds):
    from IPython.external.qt import QtCore
    QtCore.QTimer.singleShot(0, lambda: func(*args, **kwds))


if __name__ == "__main__":
    main()
