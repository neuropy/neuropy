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
from PyQt4.QtGui import QFont
#from PyQt4.QtCore import Qt
NeuropyUi, NeuropyUiBase = uic.loadUiType('neuropy.ui')

from __init__ import __version__
from animal import Animal
from track import Track
from recording import Recording

DATAPATH = os.path.expanduser('~/data')


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
        #ipyqtwidget = RichIPythonWidget(parent=self, config=config) # necessary for inline
        ipyqtwidget = RichIPythonWidget(parent=self) # necessary for inline
        self.ipyqtwidget = ipyqtwidget
        #config = load_default_config()
        #import ipdb; ipdb.set_trace()
        from IPython.utils.path import get_ipython_dir
        ipython_dir = get_ipython_dir()
        profile_dir = os.path.join(ipython_dir, 'profile_default')
        from IPython.config.loader import PyFileConfigLoader
        cl = PyFileConfigLoader('ipython_qtconsole_config.py', profile_dir)
        config = cl.load_config()
        #print config
        #import ipdb; ipdb.set_trace()
        self.ipyqtwidget.config = config

        #import ipdb; ipdb.set_trace()
        kernel_manager = QtKernelManager()
        #kernel_manager.config = config
        kernel_manager.start_kernel()
        kernel_manager.start_channels()
        ipyqtwidget.kernel_manager = kernel_manager
        ipyqtwidget.gui_completion = False
        ipyqtwidget.set_default_style(colors='linux')
        ipyqtwidget.font = QFont('Lucida Console', 12) # 3rd arg can be e.g. QFont.Bold
        ipyqtwidget.font.setFixedPitch(True)
        # make "exit" and "quit" typed in ipyqtwidget close self
        ipyqtwidget.exit_requested.connect(self.close)
        self.setCentralWidget(ipyqtwidget)

        # ip.ex() lets you execute strings in user namespace, where ip is the
        # IPython.zmq.zmqshell.ZMQInteractiveShell, as in ip = get_ipython()
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

        # don't know why I need to call "%gui qt4" here, this seems to be required as of
        # ipython 0.13.dev, in order to allow mpl figures to pop up without a .show().
        # This seemed to be automatic (via the config?) in ipython 0.12.dev:
        ipyqtwidget.execute('get_ipython().magic("gui qt4")', hidden=True)
        ipyqtwidget.execute('get_ipython().magic("pylab")', hidden=True)

        # def cf() to close all figures:
        ipyqtwidget.execute("cf = lambda: pylab.close('all')", hidden=True)
        ipyqtwidget.execute('import pylab as pl', hidden=True)
        ipyqtwidget.execute('from recording import Recording', hidden=True)
        ipyqtwidget.execute_file('globals.py', hidden=True)
        # used this instead when for some reason :execute_file wasn't working
        #for line in open('globals.py', 'r'):
        #    ipyqtwidget.execute(line, hidden=True)

        ## TODO: get neuropy imports to work properly without messing up ipython display
        ## header. Need to be done after window and ipyqtwidget init are finished. This may
        ## have been an ipython bug that's since been fixed

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
        "from animal import Animal\n"
        "try: %s; \n"
        "except NameError: %s = Animal(%r)\n"
        "%s.load(%r)" % (a.name, a.name, path, a.name, tracknames)
        )
        self.ipyqtwidget.execute(exec_lines)

    def open_track(self, path):
        tr = Track(path) # init it just to parse its id
        exec_lines = (
        "from track import Track\n"
        "tr%d = Track(%r)\n"
        "tr%d.load()" % (tr.id, path, tr.id)
        )
        self.ipyqtwidget.execute(exec_lines)

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
        <p>Copyright &copy; 2006-2012 Martin Spacek<br>
           University of British Columbia</p>
        <p>Python %s, Qt %s, PyQt %s<br>
        %s</p>""" % (__version__, platform.python_version(),
        QtCore.QT_VERSION_STR, QtCore.PYQT_VERSION_STR, platform.platform())
        QtGui.QMessageBox.about(self, "About neuropy", text)

    @QtCore.pyqtSlot()
    def on_actionAboutQt_triggered(self):
        QtGui.QMessageBox.aboutQt(self)


def set_excepthook():
    """Drops us into IPython's debugger on any error.
    If you hit "C" to continue execution after the error, this can cause the following error
    on ipython shutdown:
    <type 'exceptions.AttributeError'>: 'NoneType' object has no attribute 'ZMQError'
    If instead of you hit CTRL+D, you don't get the above error on shutdown.
    Doesn't seem like a big deal. Maybe need to restore InputHook somewhere?"""
    sys.excepthook = ultratb.FormattedTB(mode='Verbose', call_pdb=1)


if __name__ == "__main__":
    # prevents "The event loop is already running" errors when dropping into shell:
    QtCore.pyqtRemoveInputHook()
    set_excepthook()
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication([])

    neuropywindow = NeuropyWindow()
    neuropywindow.show()

    try:
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        app.exec_()

    '''
    QtCore.pyqtRemoveInputHook()
    set_excepthook()
    app = QtGui.QApplication(sys.argv)
    neuropywindow = NeuropyWindow()
    neuropywindow.show()
    sys.exit(app.exec_())
    '''
