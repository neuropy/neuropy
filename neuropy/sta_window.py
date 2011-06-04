"""STA window in Qt"""

from __future__ import division

import sys

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QPixmap, QImage
from PyQt4.QtCore import QSize

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
        nt = len(ts)
        nnids = len(nids)
        pixwidth, pixheight = 32, 32
        scale = 2

        cmap = mpl.cm.jet(np.arange(256), bytes=True) # 8 bit RGBA colormap
        #cmap[:, [0, 1, 2, 3]] = cmap[:, [3, 0, 1, 2]] # 8 bit ARGB colormap
        colortable = cmap.view(dtype=np.uint32).flatten().tolist()# const QVector<QRgb> colors 

        #print(QtGui.QColor(255, 255, 255, 0).rgba())

        layout = QtGui.QGridLayout()
        #layout.setContentsMargins(0, 0, 0, 0)
        for ti, t in enumerate(ts):
            label = QtGui.QLabel(str(t))
            layout.addWidget(label, 0, ti+1)
        for ni, nid in enumerate(nids):
            label = QtGui.QLabel('n'+str(nid))
            layout.addWidget(label, ni+1, 0)
            for ti, t in enumerate(ts):
                #data = np.uint8(np.random.randint(0, 255, size=(pixheight, pixwidth, 4))) # RGBA
                data = np.uint8(np.random.randint(0, 255, size=(pixheight, pixwidth)))
                image = QImage(data.data, pixwidth, pixheight, QImage.Format_Indexed8)
                image.ndarray = data # hold a ref, prevent gc
                image.setColorTable(colortable)
                #print(image.colorTable())
                #for ci in range(256):
                #    argb = int(cmap[ci].view(dtype=np.uint32))
                #    image.setColor(ci, argb)
                pixmap = QPixmap.fromImage(image.scaled(QSize(scale*pixwidth, scale*pixheight)))
                
                #pixmap = QPixmap(32*scale, 32*scale)
                #data = np.uint8(np.random.randint(0, 255, size=(pixheight, pixwidth, 4))) # RGBA
                #data = sp.ndimage.zoom(data, scale, order=0)
                #pixmap.loadFromData(data.data, 'uint')
                #pixmap.ndarray = data # don't let numpy array get gc'd
                # can also call image.setColorTable or setColor, and then use Format_Indexed8 format into it
                label = QtGui.QLabel()
                label.setPixmap(pixmap)
                layout.addWidget(label, ni+1, ti+1) # can also control alignment


        mainwidget = QtGui.QWidget(self)
        mainwidget.setLayout(layout) # or addlayout?

        scrollarea = QtGui.QScrollArea()
        scrollarea.setWidget(mainwidget) # or maybe instead of using setwidget, you just call setlayout or addlayout?
        
        self.setCentralWidget(scrollarea)


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    stawindow = STAWindow()
    stawindow.show()
    sys.exit(app.exec_())
