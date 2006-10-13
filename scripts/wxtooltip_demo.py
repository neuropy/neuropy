"""Example of how to use wx tooltips on a matplotlib figure window.
A tooltip pops up and tracks the mouse when the mouse is above the
axes, updating its value with the current x and y mouse position in
data coordinates. This has been tested on win32. There seem to be some
display problems under wxGTK."""

import matplotlib as mpl
mpl.use('WXAgg')
mpl.interactive(False)

import pylab as pl
from pylab import get_current_fig_manager as gcfm
import wx


class wxToolTipExample(object):

    def __init__(self):
        self.f = pl.figure()
        self.a = self.f.add_subplot(111)

        # create a long tooltip with newline to get around wx bug (in v2.6.3.3)
        # where newlines aren't recognized on subsequent self.tooltip.SetTip() calls
        self.tooltip = wx.ToolTip(tip='tip with a long %s line and a newline\n' % (' '*100))
        self.tooltip.Enable(False) # leave disabled for now
        self.tooltip.SetDelay(0) # set popup delay in ms
        gcfm().canvas.SetToolTip(self.tooltip) # connect the tooltip to the canvas

        self.f.canvas.mpl_connect('motion_notify_event', self.onmotion)

    def onmotion(self, event):
        """Called during mouse motion over figure"""
        if event.xdata != None and event.ydata != None: # mouse is inside the axes
            tip='x=%f\ny=%f' % (event.xdata, event.ydata)
            self.tooltip.SetTip(tip) # update the tooltip
            self.tooltip.Enable(True) # make sure it's enabled
        else: # mouse is outside the axes
            self.tooltip.Enable(False) # disable the tooltip

example = wxToolTipExample()

pl.show()
