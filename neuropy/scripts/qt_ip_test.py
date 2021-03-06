"""A comparison of how to create a minimalist qt IPython app using either the inprocess
kernel or the 2 process kernel. See qtconsole.examples.inprocess_qtconsole.py"""

from qtconsole.rich_ipython_widget import RichJupyterWidget
from IPython.lib import guisupport

INPROCESS = False
CLEANSHUTDOWN = False # clean is slow, unclean is fast

def main():
    """Start kernel manager and client, create window, run app event loop"""
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

    widget = RichJupyterWidget()
    widget.kernel_manager = km
    widget.kernel_client = kc
    if CLEANSHUTDOWN: # slow exit on CTRL+D
        def stop():
            kc.stop_channels()
            km.shutdown_kernel()
            app.exit()
        widget.exit_requested.connect(stop)
    else: # fast exit on CTRL+D
        widget.exit_requested.connect(app.quit)
    widget.show()
    guisupport.start_event_loop_qt4(app)

if __name__ == "__main__":
    main()
