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
