def rf():
    """Test plotting RFs in a mpl figure.
    Could also do a pure Qt window, either by hijacking the mpl figure,
    or by mimicking the displayhook thing thing that mpl probably does
    when interacting with ipython. Then, in Qt, you can render any widget,
    including all its contained widgets, to an SVG."""
    ts = np.arange(10)
    nids = np.arange(20)
    nt = len(ts)
    nnids = len(nids)
    
    pixwidth, pixheight = 32, 32
    scale = 2
    dpi = 80
    
    # inches
    xoff = 0.5 # makes space for vertical column of nids
    yoff = 0.25 # makes space for horizontal row of ts
    rfwidth = pixwidth * scale / dpi
    rfheight = pixheight * scale / dpi
    xsep = 0.05
    ysep = 0.05
    fs = xoff+nt*(rfwidth + xsep), yoff+nnids*(rfheight + ysep) # figure size
    print(fs)

    import pylab as pl
    import matplotlib as mpl
    import scipy as sp
    import scipy.ndimage
    
    #interactive = mpl.rcParams['interactive']
    #mpl.rcParams['interactive'] = False
    
    #f = mpl.figure.Figure(figsize=fs) # doesn't create a canvas
    f = pl.figure(figsize=fs)#, dpi=80, frameon=True)
    try:
        f.canvas.window().statusBar().hide() # Qt specific
        f.canvas.toolbar.hide() # not Qt specific
        window = f.canvas.window()
    except AttributeError: pass
    #pixmap = QtGui.QPixmap
    #data = np.uint8(np.random.rand(0, 255, size=(32, 32, 4))) # RGBA
    
    #pixmap = QtGui.QPixmap(32, 32)
    #pixmap.loadFromData(data.data, data.nbytes)
    #image = QtGui.QImage(Format_Mono
    #window.setCentralWidget(image)
    
    #a = f.add_axes((0, 0, 1, 1))

    # create ts labels above RFs
    for ti, t in enumerate(ts):
        x = (xoff + ti*(rfwidth + xsep)) / fs[0]
        y = 1 - yoff / fs[1]
        f.text(x, y, str(t), #transform=f.transFigure,
               horizontalalignment='left', verticalalignment='bottom')
    # create nid labels left of RFs
    for ni, nid in enumerate(nids):
        x = xoff / fs[0]
        y = 1 - (yoff + rfheight/2 + ni*(rfheight + ysep)) / fs[1]
        f.text(x, y, 'n'+str(nid), #transform=f.transFigure,
               horizontalalignment='right', verticalalignment='center')
    # plot RFs
    #awidth = rfwidth / fs[0]
    #aheight = rfheight / fs[1]
    #a = mpl.axes.Axes(f, (0.5, 0.5, 0.25, 0.25))#, frameon=False)
    
    for ni, nid in enumerate(nids):
        offsety = (yoff + rfheight + ni*(rfheight + ysep)) * dpi
        for ti, t in enumerate(ts):
            offsetx = (xoff + ti*(rfwidth + xsep)) * dpi
            #a = f.add_axes((x, y, awidth, aheight))
            data = np.random.rand(pixheight, pixwidth)
            data = sp.ndimage.zoom(data, scale, order=0)
            fi = mpl.image.FigureImage(f, offsetx=offsetx, offsety=offsety, data=data)
            f.images.append(fi)
            #a.imshow(data, interpolation='nearest')

    #mpl.rcParams['interactive'] = interactive
    #pl.draw()
