"""STA window in Qt. I think what's necessary to make this work in
ipython-qtconsole is to have neuropy (via ipython kernel somehow) send
a zmq message (display_data) back to the frontend telling it to create
an STA window, as well as the results of the STA calc in the message payload
or something. mpl must be doing this every time you make a plot, 
but have to verify"""

from __future__ import division

import sys

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QPixmap, QImage, QPalette, QColor
from PyQt4.QtCore import Qt, QSize

import numpy as np
import scipy as sp
import scipy.ndimage
import matplotlib as mpl
import matplotlib.cm


class STAWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self)
        ts = np.arange(0, 200, 20)
        nids = np.arange(100)
        width, height = 32, 32
        scale = 2 # setting to a float will give uneven sized pixels

        cmap = mpl.cm.jet(np.arange(256), bytes=True) # 8 bit RGBA colormap
        #cmap[:, [0, 1, 2, 3]] = cmap[:, [3, 0, 1, 2]] # 8 bit ARGB colormap
        # from Qt docs, sounds like I should be using ARGB format, but seems like
        # RGBA is the format that works in PyQt4
        colortable = cmap.view(dtype=np.uint32).ravel().tolist() # QVector<QRgb> colors 
        layout = QtGui.QGridLayout() # can set vert and horiz spacing
        #layout.setContentsMargins(0, 0, 0, 0) # doesn't seem to do anything

        for ti, t in enumerate(ts): # time labels along top
            label = QtGui.QLabel(str(t))
            layout.addWidget(label, 0, ti+1)
        for ni, nid in enumerate(nids):
            label = QtGui.QLabel('n'+str(nid)) # nid label on left
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            layout.addWidget(label, ni+1, 0)
            for ti, t in enumerate(ts):
                data = np.uint8(np.random.randint(0, 255, size=(height, width)))
                image = QImage(data.data, width, height, QImage.Format_Indexed8)
                image.ndarray = data # hold a ref, prevent gc
                image.setColorTable(colortable)
                image = image.scaled(QSize(scale*width, scale*height)) # scale it
                pixmap = QPixmap.fromImage(image)
                label = QtGui.QLabel()
                label.setPixmap(pixmap)
                layout.addWidget(label, ni+1, ti+1) # can also control alignment

        mainwidget = QtGui.QWidget(self)
        mainwidget.setLayout(layout)

        scrollarea = QtGui.QScrollArea()
        scrollarea.setWidget(mainwidget)

        self.setCentralWidget(scrollarea)
        #palette = QPalette(QColor(255, 255, 255))
        #self.setPalette(palette) # set white background, or perhaps more


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    stawindow = STAWindow()
    stawindow.show()
    sys.exit(app.exec_())
