import wx
import numpy as np

class MseqFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwargs)
        self.panel_1 = wx.ScrolledWindow(self, -1, style=wx.TAB_TRAVERSAL)

        self.ns = range(40)
        self.ts = range(0, 200, 20)
        self.bitmaps = {}
        for ni, n in enumerate(self.ns):
            self.bitmaps[n] = {}
            for ti, t in enumerate(self.ts):
                a = np.rand(32*3, 32)*255
                a = a.round().astype(np.uint8)
                im = wx.ImageFromData(32, 32, a.data)
                im = im.Scale(64, 64)
                self.bitmaps[n][t] = wx.StaticBitmap(parent=self.panel_1, bitmap=im.ConvertToBitmap())

        #self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.__set_properties()
        self.__do_layout()
    '''
    def OnPaint(self, event):
        #self.canvas.draw()
        event.Skip()
    '''
    def __set_properties(self):
        self.SetTitle("MseqFrame")
        self.panel_1.SetScrollRate(10, 10)

    def __do_layout(self):
        sizer_1 = wx.GridSizer(1, 1, 0, 0)
        grid_sizer_1 = wx.FlexGridSizer(rows=len(self.ns)+1, cols=len(self.ts)+1, vgap=2, hgap=2)
        grid_sizer_1.Add((1, 1), 0, wx.ADJUST_MINSIZE, 0) # spacer in top left corner
        for t in self.ts:
            grid_sizer_1.Add(wx.StaticText(self.panel_1, -1, "%sms" % t), 0, wx.ADJUST_MINSIZE, 0) # text row along top
        for n in self.ns:
            grid_sizer_1.Add(wx.StaticText(self.panel_1, -1, "n%d" % n), 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL|wx.ADJUST_MINSIZE, 0) # text down left side
            for t in self.ts:
                grid_sizer_1.Add(self.bitmaps[n][t], 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL, 0)
        self.panel_1.SetAutoLayout(True)
        self.panel_1.SetSizer(grid_sizer_1)
        grid_sizer_1.Fit(self.panel_1)
        #grid_sizer_1.SetSizeHints(self.panel_1) # prevents the panel from being resized to something smaller than the above fit size
        '''
        # might be a more direct way to set these:
        for rowi in range(1, len(self.ns)+1):
            print 'rowi:', rowi
            grid_sizer_1.AddGrowableRow(rowi)
        for coli in range(1, len(self.ts)+1):
            print 'coli:', coli
            grid_sizer_1.AddGrowableCol(coli)
        '''
        sizer_1.Add(self.panel_1, 1, wx.ADJUST_MINSIZE|wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        #sizer_1.SetSizeHints(self) # prevents the frame from being resized to something smaller than the above fit size
        self.Layout()


if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    frame_1 = MseqFrame(None, -1, "")
    app.SetTopWindow(frame_1)
    frame_1.Show()
    app.MainLoop()
