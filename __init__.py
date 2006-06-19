"""Tools for analyzing neuronal and stimulus data, in python"""

try:
    reload(Core)
except NameError:
    import Core
from Core import *

try:
    reload(Test)
except NameError:
    import Test
from Test import test

__importall__ = [Core] # list of mods imported as "from mod import *"
__all__ = []
import types
for mod in __importall__:
    for key, value in mod.__dict__.iteritems():
        if type(value) is not types.ModuleType and not key.startswith('_'):
            #print key
            __all__.append(key)
'''
try:
    reload(Core)
    #__main__.refreshmod(Core, importall=True)
except:
    print 'couldn''t reload Core'
    import Core
from Core import *

try:
    __main__.refreshmod(Test, importall=False)
except:
    print 'couldn''t reload Test'
    import Test
'''
