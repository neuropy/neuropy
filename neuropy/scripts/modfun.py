def delmod(mod, importall=False):
    import copy
    d = mod.__dict__.copy()
    dk = copy.copy(d.keys())
    for key in dk:
        if '__name__' not in key:# and type(nd[key]) is not types.ModuleType:
            del d[key]
            print 'deleted %s key:' % mod.__name__, key

    if importall:
        import __main__
        md = __main__.__dict__
        print dk
        for key in dk:
            if not key.startswith('__'):# and type(md[key]) is not types.ModuleType:
                try:
                    del md[key]
                    print 'deleted md key:', key
                except:
                    pass

    del md[mod.__name__]

def refreshmod(mod, importall=False):
    """Reloads all neuropy names into namespace"""

    delmod(mod, importall)

    try:
        exec 'reload(%s)' % mod.__name__ in globals()
    except:
        exec 'import %s' % mod.__name__ in globals()
        exec 'reload(%s)' % mod.__name__ in globals()
    if importall:
        exec 'from %s import *' % mod.__name__ in globals() # a direct import * is apparently not allowed in a function


'''
# set neuropy's __all__ attrib for designating which names will be imported when
# you go 'from neuropy import *'
__importall__ = [core] # list of mods imported above as "from mod import *"
__all__ = ['Neuron', 'Experiment', 'test']
import types
for mod in __importall__:
    for key, value in mod.__dict__.iteritems():
        # don't add modules or special names designated with a leading underscore
        if type(value) is not types.ModuleType and not key.startswith('_'):
            #print key
            __all__.append(key)
'''
