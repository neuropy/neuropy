"""Main neuropy window"""

__authors__ = ['Martin Spacek']

def disable(*modules):
    """ Simulate the non-existance of modules, taken from:
    https://gist.github.com/ChrisBeaumont/5048768
 
    :param modules: One or more modules to mock-uninstall, for the rest of the
                    process
    """
    class ID(object):
        def __init__(self, modules):
            self.forbidden = modules
        def find_module(self, mod_name, pth):
            if mod_name in self.forbidden:
                return self
        def load_module(self, mod_name):
            raise ImportError("forbidden")
 
    id = ID(modules)
    import sys
    sys.meta_path.append(id)

# in case PySide is installed, hide it to prevent conflicts between it and PyQt4 in IPython:
disable('PySide')

import os
import sys
import platform

# instantiate an IPython embedded shell which shows up in the terminal on demand
# and on every exception in the GUI code:
import IPython
from IPython.terminal.ipapp import load_default_config
from IPython.terminal.embed import InteractiveShellEmbed
config = load_default_config()
# automatically call the pdb debugger after every exception, override default config:
config.TerminalInteractiveShell.pdb = True
ipshell = InteractiveShellEmbed(display_banner=False, config=config)

# IPython GUI imports, have to come before Qt imports:
from qtconsole.rich_jupyter_widget import RichJupyterWidget
from IPython.lib import guisupport
from IPython.paths import get_ipython_dir
from traitlets.config.loader import PyFileConfigLoader

from qtconsole.qt import QtCore, QtGui
from PyQt4 import uic
NeuropyUi, NeuropyUiBase = uic.loadUiType('neuropy.ui')

from __init__ import __version__
from animal import Animal
from track import Track
from recording import Recording

from globals import DATAPATH, BLABPATH, SLABPATH
# use inprocess kernel? otherwise, use 2 process kernel. This option doesn't work for now:
INPROCESS = False


class NeuropyWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        """Note that a lot of default options can be changed by editing
        IPython.qtconsole.console_widget.py and ipython_widget.py"""
        QtGui.QMainWindow.__init__(self)
        self.ui = NeuropyUi()
        self.ui.setupUi(self) # lay it out

        # gui_completion: 'plain, 'droplist' or 'ncurses'
        ipw = RichJupyterWidget(parent=self, gui_completion='plain')
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
        if tail[0].isnumeric(): # asume it's a recording
            self.open_recording(path)
        elif tail.startswith('tr'): # assume it's a track
            self.open_track(path)
        else:
            self.open_animal(path) # assume it's an animal

    @QtCore.pyqtSlot()
    def on_action_ptc15_triggered(self):
        path = os.path.join(SLABPATH, 'ptc15')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc17_triggered(self):
        path = os.path.join(SLABPATH, 'ptc17')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc18_triggered(self):
        path = os.path.join(SLABPATH, 'ptc18')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc20_triggered(self):
        path = os.path.join(SLABPATH, 'ptc20')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc21_triggered(self):
        path = os.path.join(SLABPATH, 'ptc21')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc22_triggered(self):
        path = os.path.join(SLABPATH, 'ptc22')
        self.open_animal(path)

    @QtCore.pyqtSlot()
    def on_action_ptc15_tr7c_triggered(self):
        path = os.path.join(SLABPATH, 'ptc15')
        self.open_animal(path, 'tr7c')

    @QtCore.pyqtSlot()
    def on_action_ptc18_tr1_triggered(self):
        path = os.path.join(SLABPATH, 'ptc18')
        self.open_animal(path, 'tr1')

    @QtCore.pyqtSlot()
    def on_action_ptc18_tr2c_triggered(self):
        path = os.path.join(SLABPATH, 'ptc18')
        self.open_animal(path, 'tr2c')


    @QtCore.pyqtSlot()
    def on_action_ptc22_tr1_triggered(self):
        path = os.path.join(SLABPATH, 'ptc22')
        self.open_animal(path, 'tr1')

    @QtCore.pyqtSlot()
    def on_action_ptc22_tr2_triggered(self):
        path = os.path.join(SLABPATH, 'ptc22')
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
        "tr%s = Track(%r)\n"
        "tr%s.load()" % (tr.id, path, tr.id)
        )
        self.ipw.execute(exec_lines)

    def open_recording(self, path):
        rec = Recording(path) # init it just to parse its id
        exec_lines = (
        "r%s = Recording(%r)\n"
        "r%s.load()" % (rec.id, path, rec.id)
        )
        self.ipw.execute(exec_lines)

    @QtCore.pyqtSlot()
    def on_actionShell_triggered(self):
        ipshell()

    @QtCore.pyqtSlot()
    def on_actionQuit_triggered(self):
        self.close() # trigger closeEvent

    def closeEvent(self, event):
        self.stop()

    def stop(self):
        """Taken from ipython/examples/Embedding/inprocess_qtconsole.py. Ensure that
        both the kernel process and the main Python app that launched it are killed."""
        # greatly slows down shutdown, doesn't seem necessary:
        #self.ipw.kernel_client.stop_channels()
        self.ipw.kernel_manager.shutdown_kernel(now=True) # now=True shuts down faster
        guisupport.get_app_qt4().exit()

    @QtCore.pyqtSlot()
    def on_actionAboutNeuropy_triggered(self):
        lf = open('../LICENSE', 'r')
        LICENSE = lf.read()
        lf.close()
        ipverstr = ('%d.%d.%d.%s' % IPython.version_info).rstrip('.')
        system = """<p>Python %s, IPython %s, Qt %s, PyQt %s<br>
                    %s</p>""" % (platform.python_version(), ipverstr,
                                 QtCore.QT_VERSION_STR, QtCore.PYQT_VERSION_STR,
                                 platform.platform())
        text = """
        <h2><a href=http://neuropy.github.io>neuropy</a> %s</h2>
        <p>A tool for neuronal spike, stimulus, and LFP data analysis</p>

        <p>Copyright &copy; 2006-2017 <a href=http://mspacek.github.io>Martin Spacek</a><br>
        
           <a href=http://swindale.ecc.ubc.ca>Swindale</a> Lab,
           University of British Columbia, Vancouver, Canada<br>

           <a href=http://www.neuro.bio.lmu.de/members/system_neuro_busse/busse_l/index.html>
           Busse</a> Lab, Ludwig-Maximilians-University, Munich, Germany</p>

        <p>%s</p>

        %s""" % (__version__, LICENSE, system)
        QtGui.QMessageBox.about(self, "About neuropy", text)

    @QtCore.pyqtSlot()
    def on_actionAboutQt_triggered(self):
        QtGui.QMessageBox.aboutQt(self)


def config_ipw(ipw):
    """Apply and then modify default settings of IPython Qt widget"""
    ipython_dir = get_ipython_dir()
    profile_dir = os.path.join(ipython_dir, 'profile_default')
    cl = PyFileConfigLoader('ipython_config.py', profile_dir)
    config = cl.load_config()
    ipw.config = config

    ipw.set_default_style(colors='Linux')
    ipw.font = QtGui.QFont('Lucida Console', 11) # 3rd arg can be e.g. QFont.Bold
    ipw.font.setFixedPitch(True)

def main():
    """Start kernel manager and client, create window, run app event loop,
    auto execute some code in user namespace. A minimalist example is shown in
    qt_ip_test.py.

    NOTE: Make sure that the Qt v2 API is being used by IPython by running `export
    QT_API=pyqt` at the command line before running neuropy, or by adding it to `.bashrc`"""
    app = guisupport.get_app_qt4()

    if INPROCESS:
        from qtconsole.inprocess import QtInProcessKernelManager
        km = QtInProcessKernelManager()
    else:
        from qtconsole.manager import QtKernelManager
        km = QtKernelManager()
    km.start_kernel()
    km.kernel.gui = 'qt4'
    kc = km.client()
    kc.start_channels()

    nw = NeuropyWindow()
    ipw = nw.ipw
    config_ipw(ipw)
    ipw.kernel_manager = km
    ipw.kernel_client = kc
    ipw.exit_requested.connect(nw.stop)
    nw.show()

    # execute some code through the frontend (once the event loop is running).
    # The output appears in the IPythonWidget (ipw).
    do_later(ipw.execute, 'run -i %s' % 'startup.py', hidden=True)
    do_later(ipw.execute, 'run -i %s' % 'globals.py', hidden=True)
    #guisupport.start_event_loop_qt4(app) # doesn't seem to work in IPy 5.3.0
    app.exec_()

def do_later(func, *args, **kwargs):
    QtCore.QTimer.singleShot(0, lambda: func(*args, **kwargs))


if __name__ == "__main__":
    # prevents "The event loop is already running" messages when calling ipshell(),
    # doesn't seem to be necessary any more in IPy 5.3.0:
    #QtCore.pyqtRemoveInputHook()
    main()
