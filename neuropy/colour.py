"""Define colour-related things"""

from __future__ import division
from __future__ import print_function

from copy import copy

import numpy as np

RED = '#ff0000'
ORANGE = '#ff7f00'
YELLOW = '#ffff00'
DARKYELLOW = '#efef00'
GREEN = '#00ff00'
CYAN = '#00ffff'
LIGHTBLUE = '#007fff'
MIDBLUE = '#0040ff'
BLUE = '#0000ff'
VIOLET = '#7f00ff'
MAGENTA = '#ff00ff'
GREY = '#7f7f7f'
WHITE = '#ffffff'
BROWN = '#af5050'
DARKGREY = '#303030'
LIGHTBLACK = '#202020'
BLACK = '#000000'


class ColourDict(dict):
    """Just an easy way to cycle through colours given some index,
    like say a chan id or a neuron id. Better than using a generator,
    cuz you don't need to keep calling .next(). This is like a dict
    of infinite length. Copied from spyke.plot"""
    def __init__(self, colours=None, nocolour=None):
        self.colours = colours
        self.nocolour = nocolour

    def __getitem__(self, key):
        if key < 0: # invalid index into self.colours
            return self.nocolour
        i = key % len(self.colours) - 1 # convert 1-based indices into 0-based
        return self.colours[i]

    def __setitem__(self, key, val):
        raise RuntimeError('ColourDict is unsettable')


def hex2floatrgb(s):
    """Convert RGB hex string s into an RGB float array"""
    s = s[len(s)-6:len(s)] # get last 6 characters
    r, g, b = s[0:2], s[2:4], s[4:6]
    r, g, b = int(r, base=16), int(g, base=16), int(b, base=16)
    return np.float64([r, g, b]) / 255


# cluster colours for plotting on a white background:
CCWHITE = [RED, ORANGE, DARKYELLOW, GREEN, CYAN, LIGHTBLUE, VIOLET, MAGENTA, BLACK, BROWN]
CCWHITEDICT = ColourDict(colours=CCWHITE)
CCWHITERGB = [ hex2floatrgb(c) for c in CCWHITE ]
CCWHITERGBDICT = ColourDict(colours=CCWHITERGB)

# cluster colours for plotting on a black background:
CCBLACK = [RED, ORANGE, YELLOW, GREEN, CYAN, LIGHTBLUE, VIOLET, MAGENTA, WHITE, BROWN]
CCBLACKDICT = ColourDict(colours=CCBLACK)
CCBLACKRGB = [ hex2floatrgb(c) for c in CCBLACK ]
CCBLACKRGBDICT = ColourDict(colours=CCBLACKRGB)
